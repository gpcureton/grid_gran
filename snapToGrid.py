#!/usr/bin/env python
# encoding: utf-8
"""
snapToGrid.py

This module contains the SnapToGrid class, which allows the regridding 
of a dataset from one grid to another.

This python class implements a gridding procedure developed
by Nick Bearson, and modified by Nadia Smith, both from SSEC.

Created by Geoff Cureton on 2010-04-17.
Copyright (c) 2010 University of Wisconsin SSEC. All rights reserved.
"""

import string, sys
from glob import glob
from os import path,uname
from time import time
import shlex, subprocess

import numpy as np
from matplotlib.mlab import griddata
import scipy.weave as weave
from scipy.weave import converters

import tables as pytables

import ctypes
from numpy.ctypeslib import ndpointer


class SnapToGrid:

    """
    SnapToGrid


    """

    def __init__(self,interp='nn'):
        """
        __init__

        geoFileList: Common list of geolocation files used for every dataset 
                     in this object.
        dataFileList: Dictionary of one or more datasets, each containing files
                      corresponding to those in geoFileList.
        gridLat,gridLon : Equal angle lat and lon grids which the data files are 
                          gridded to.
        grid2GranFileIdx: At each grid point, grid2GranFileIdx gives the index in 
                          geoFileList (and in any dataFileList datafile lists) of 
                          the file which was gridded to that point. 
        grid2GranIdx: At each grid point, grid2GranIdx gives the index of the pixel 
                      which was gridded to that grid point, from the file given by 
                      dataFileList.

        """

        # A default 1 degree grid...
        degInc = 1.
        grids = np.mgrid[-90.:90.+degInc:degInc,-180.:180.:degInc]
        gridLat,gridLon = grids[0],grids[1]

        self.interp=interp
        self.gridLat = gridLat
        self.gridLon = gridLon
        self.gridData = {}
        self.geoFileList = []
        self.dataFileList = {}
        self.grid2GranIdx = np.ones(np.shape(gridLat),dtype=np.int64) * -99999
        self.grid2GranFileIdx = np.ones(np.shape(gridLat),dtype=np.int64) * -99999

    @staticmethod
    def fileToGridObj(gridFile):
        """
        fileToGridObj

        gridFile: Name of file containing gridded data

        """
        print "Reading the file %s" % (gridFile)
        fileObj = pytables.openFile(gridFile,"r")

        # Create a new object
        gridObj = SnapToGrid()

        # Read file related datasets into SnapToGrid object
        gridObj.geoFileList = list(fileObj.getNode("/fileData/geolocationFiles")[:])

        gridObj.grid2GranIdx = fileObj.getNode("/fileData/grid2GranIdx")[:,:]
        gridObj.grid2GranFileIdx = fileObj.getNode("/fileData/grid2GranFileIdx")[:,:]

        for node in fileObj.walkNodes('/fileData/dataFileLists',classname='Array') :
            gridObj.dataFileList[node.name] = list(node[:])

        # Read grid related datasets into SnapToGrid object
        gridObj.gridLat = fileObj.getNode("/gridData/Latitude")[:,:]
        gridObj.gridLon = fileObj.getNode("/gridData/Longitude")[:,:]
        try :
            gridObj.gridDegen = fileObj.getNode("/gridData/gridDegen")[:,:]
        except :
            pass

        for node in fileObj.walkNodes('/gridData/gridDataSets',classname='Array') :
            gridObj.gridData[node.name] = node[:,:]

        fileObj.close()

        return gridObj

    def gridObjToFile(self,gridFile='gridFile.h5'):
        """
        gridObjToFile

        gridFile: Name of file containing gridded data

        """
        print "Creating the file %s" % (gridFile)
        fileObj = pytables.openFile(gridFile,"w")

        # Write the grid latitude and longitude
        print "Creating the Latitude and Longitude datasets"
        fileObj.createArray("/gridData","Latitude",self.gridLat, \
                    createparents=True)
        fileObj.createArray("/gridData","Longitude",self.gridLon, \
                    createparents=True)
        try :
            fileObj.createArray("/gridData","gridDegen",self.gridDegen, \
                        createparents=True)
        except :
            pass

        # Write the various gridded datasets
        for dSet in self.gridData.keys() :
            print "Dataset : ",dSet
            fileObj.createArray("/gridData/gridDataSets",dSet,self.gridData[dSet], \
                    createparents=True)

        # Write the list of geolocation files
        print "Creating the geolocationFiles dataset"
        fileObj.createArray("/fileData","geolocationFiles", \
                np.array(self.geoFileList,dtype=np.str), \
                createparents=True)

        # Write grid2GranFileIdx and grid2GranIdx
        print "Creating the grid2GranFileIdx and grid2GranIdx datasets"
        fileObj.createArray("/fileData","grid2GranFileIdx",self.grid2GranFileIdx, \
                    createparents=True)
        fileObj.createArray("/fileData","grid2GranIdx",self.grid2GranIdx, \
                    createparents=True)

        # Write the various datasets lists
        for dSet in self.dataFileList.keys() :
            print "Dataset : ",dSet
            fileObj.createArray("/fileData/dataFileLists",dSet, \
                    np.array(self.dataFileList[dSet],dtype=np.str), \
                    createparents=True)


        # Flush the file, and close it
        fileObj.flush()
        fileObj.close()


    def __call__(self, lat, lon, data, gridLat, gridLon):
        """
        __call__

        This class special method takes as arguments the latitude 
        and longitude arrays, the data array which we wish to 
        regrid, and the lat and lon grids we are mapping to.

        returns: numpy array of regridded data. Grid cells not
        covered by original data are masked.
        """
        return self.__granToGrid(lat,lon,data,gridLat,gridLon)

    def __granToGrid(self, lat, lon, data, gridLat, gridLon):
        """
        __granToGrid

        This private class method takes as arguments the latitude 
        and longitude arrays, the data array which we wish to 
        regrid, and the lat and lon grids we are mapping to.

        returns: numpy array of regridded data. Grid cells not
        covered by original data are masked.
        """
        newData = griddata(lon,lat,data,gridLon,gridLat,interp=self.interp)
        return newData

    @staticmethod
    def regrid(lat, lon, data, gridLat, gridLon,interp='nn'):
        """
        regrid

        This static class method takes as arguments the latitude 
        and longitude arrays, the data array which we wish to regrid, 
        and the lat and lon grids we are mapping to.

        returns: numpy array of regridded data. Grid cells not
        covered by original data are masked.
        """
        newData = griddata(lon,lat,data,gridLon,gridLat,interp=interp)
        return newData

    @staticmethod
    def snapGrid_numpy(lat, lon, data, gridLat, gridLon, gridData, saveGrid=True):
        """
        snapGrid (numpy)

        This static class method takes as arguments the latitude 
        and longitude arrays, the data array which we wish to regrid, 
        and the lat and lon grids we are mapping to.

        Input...
            lat: 1D array of latitude values
            lon: 1D array of longitude values
            data: 1D array of data corresponding to lat and lon
            gridLat: 2D array defining new latitude grid
            gridLon: 2D array defining new longitude grid
            gridData: 2D array defining output data grid. May already 
                have valid data
            saveGrid: Boolean switch for whether new grid indicies
                  are saved and returned.

        Returns...
            gridData: 2D array containing input data gridded to output data grid 
            dataIdx:  2D array containing indicies of input data which are gridded 
                      to output grid

        """

        print "Shape of data is",np.shape(data)
        print "Shape of gridData is",np.shape(gridData)
        
        dataIdx = np.ones(np.shape(gridData),dtype=np.int64) * -99999

        gridLatInc = np.abs(gridLat[1,0]-gridLat[0,0])
        gridLonInc = np.abs(gridLon[0,1]-gridLon[0,0])

        for idx in range(np.size(lat)):
            latVal,lonVal,dataVal = lat[idx],lon[idx],data[idx]

            # Determine bracketing indicies
            latGridIdxLo = int(np.floor((latVal-gridLat[0,0])/gridLatInc))
            latGridIdxHi = latGridIdxLo + 1
            lonGridIdxLo = int(np.floor((lonVal-gridLon[0,0])/gridLonInc))
            lonGridIdxHi = lonGridIdxLo + 1

            gridPoints = (
                [ latGridIdxLo, latGridIdxLo, latGridIdxHi, latGridIdxHi ],
                [ lonGridIdxLo, lonGridIdxHi, lonGridIdxLo, lonGridIdxHi ]
            )

            minDist = 1000.

            try :
                for points in zip(gridPoints[0],gridPoints[1]) :
                    dist = np.sqrt((latVal-gridLat[points])**2. + (lonVal-gridLon[points])**2.)
                    if (dist < minDist) :
                        snapPoints = points
                        minDist = dist

                gridData[snapPoints] = dataVal
            except IndexError :
                #print "Index error for point: ",points
                pass

            dataIdx[snapPoints] = idx

        return gridData,dataIdx

    @staticmethod
    def snapGrid_weave(lat, lon, data, gridLat, gridLon, gridData, \
        gridDataIdx, fillVal, cellAverage=False):
        """
        snapGrid (weave)

        This static class method takes as arguments the latitude 
        and longitude arrays, the data array which we wish to regrid, 
        and the lat and lon grids we are mapping to.

        Input...
            lat: 1D array of latitude values
            lon: 1D array of longitude values
            data: 1D array of data corresponding to lat and lon
            gridLat: 2D array defining new latitude grid
            gridLon: 2D array defining new longitude grid
            gridData: 2D array defining output data grid. May already 
                have valid data
            gridDegen: 2D array giving the degeneracy of each griddpoint.
            gridDataIdx: 2D array containing indicies of input data which are gridded 
                      to output grid
            fillVal: The desired fill value for this grid type.

        Returns...
            gridData: 2D array containing input data gridded to output data grid.
                      Same as gridData
            dataIdx:  2D array containing indicies of input data which are gridded 
                      to output grid. Same as gridDataIdx

        """

        codeSnapGrid = """

            using namespace std;

            long nData, nGridRows,nGridCols, idx;
            int latGridIdxLo, latGridIdxHi, lonGridIdxLo, lonGridIdxHi;

            double latVal,lonVal,dataVal;
            double gridLatInc,gridLonInc;

            int gridLatPt,gridLonPt;
            int gridLatPoints[4], gridLonPoints[4];
            double minDist,dist,latDist,lonDist;
            double gridLatVal,gridLonVal;
            int crnrPt,snapCrnrPt;
            bool rowInBounds,colInBounds;

            nData = (long) Ndata[0];
            nGridRows = (long) NgridData[0];
            nGridCols = (long) NgridData[1];

            printf("Shape of data is (%ld,)\\n",nData);
            printf("Shape of gridData is (%ld, %ld)\\n", nGridRows,nGridCols);

            // Determine the lat and lon grid spacings
            gridLatInc = fabs(gridLat(1,0)-gridLat(0,0));
            gridLonInc = fabs(gridLon(0,1)-gridLon(0,0));

            // Loop through original data, finding matching gridpoint and assigning
            // data value to that point

            for (idx=0;idx<nData;idx++){
                latVal = lat(idx);
                lonVal = lon(idx);
                dataVal = data(idx);

                latGridIdxLo = (int) floor((latVal-gridLat(0,0))/gridLatInc);
                latGridIdxHi = latGridIdxLo + 1;
                lonGridIdxLo = (int) floor((lonVal-gridLon(0,0))/gridLonInc);
                lonGridIdxHi = lonGridIdxLo + 1;

                rowInBounds = true;
                colInBounds = true;

                if ((latGridIdxLo<0) || (latGridIdxHi>=nGridRows)){
                    rowInBounds = false;
                }
                if ((lonGridIdxLo<0) || (lonGridIdxHi>=nGridCols)){
                    colInBounds = false;
                }

                if (!rowInBounds){
                    continue;
                }else if (!colInBounds){
                    continue;
                }else{
                    gridLatPoints[0] = latGridIdxLo;
                    gridLatPoints[1] = latGridIdxLo;
                    gridLatPoints[2] = latGridIdxHi;
                    gridLatPoints[3] = latGridIdxHi;

                    gridLonPoints[0] = lonGridIdxLo;
                    gridLonPoints[1] = lonGridIdxHi;
                    gridLonPoints[2] = lonGridIdxLo;
                    gridLonPoints[3] = lonGridIdxHi;

                    minDist = 1000.;
                    snapCrnrPt = 0;

                    for (crnrPt=0;crnrPt<4;crnrPt++){

                        gridLatPt = (int) gridLatPoints[crnrPt];
                        gridLonPt = (int) gridLonPoints[crnrPt];

                        gridLatVal = gridLat(gridLatPt,gridLonPt);
                        gridLonVal = gridLon(gridLatPt,gridLonPt);

                        latDist = latVal-gridLatVal;
                        lonDist = lonVal-gridLonVal;

                        dist = sqrt(latDist*latDist + lonDist*lonDist); 

                        if (dist < minDist){
                            snapCrnrPt = crnrPt;
                            minDist = dist;
                        }
                    }

                    gridLatPt = (int) gridLatPoints[snapCrnrPt];
                    gridLonPt = (int) gridLonPoints[snapCrnrPt];

                    // Assign data value to this grid cell
                    if (cellAverage){
                        if (fabs(gridData(gridLatPt,gridLonPt) - fabs(fillVal)) < 0.001){
                            gridData(gridLatPt,gridLonPt) = dataVal;
                        }else{
                            gridData(gridLatPt,gridLonPt) += dataVal;
                        }
                    }else{
                        gridData(gridLatPt,gridLonPt) = dataVal;
                    }
                    gridDegen(gridLatPt,gridLonPt) += 1;

                    // TODO : Save subsequent data indices in the same gridcell
                    gridDataIdx(gridLatPt,gridLonPt) = idx;
                }

            }


        """

        gridDegen = -999 * np.ones(np.shape(gridData),dtype=np.int32)
        cellAverage = int(cellAverage)

        weave.inline(codeSnapGrid,
            arg_names=['lat','lon','data',\
                       'gridLat','gridLon','gridData',\
                       'gridDataIdx','cellAverage','gridDegen','fillVal'],
            type_converters=converters.blitz,
            headers=['<math.h>'],
            libraries=['m'],
            #include_dirs=self.include_dirs,
            force=0)

        if cellAverage :
            return gridData,gridDataIdx,gridDegen
        else :
            return gridData,gridDataIdx

    @staticmethod
    def snapGrid_ctypes(dataLat, dataLon, data, gridLat, gridLon, gridData, gridDataIdx):
        """
        snapGrid (ctypes)

        This static class method takes as arguments the latitude 
        and longitude arrays, the data array which we wish to regrid, 
        and the lat and lon grids we are mapping to.

        Input...
            lat: 1D array of latitude values
            lon: 1D array of longitude values
            data: 1D array of data corresponding to lat and lon
            gridLat: 2D array defining new latitude grid
            gridLon: 2D array defining new longitude grid
            gridData: 2D array defining output data grid. May already 
                have valid data
            gridDataIdx: 2D array containing indicies of input data which are gridded 
                      to output grid

        Returns...
            gridData: 2D array containing input data gridded to output data grid.
                      Same as gridData
            dataIdx:  2D array containing indicies of input data which are gridded 
                      to output grid. Same as gridDataIdx

        """

        nData = np.int64(data.size)
        gridRows = np.int32(gridLat.shape[0])
        gridCols = np.int32(gridLat.shape[1])

        gridData = np.float64(gridData)
        gridDataIdx  = np.int64(gridDataIdx)

        libDir = path.dirname(__file__)
        libFile = 'libgriddingAndGranulation.so.1.0.1'
        libFile = "%s/%s" % (libDir,libFile)
        lib = ctypes.cdll.LoadLibrary(libFile)

        snapGrid_ctypes = lib.gran2grid
        snapGrid_ctypes.restype = None
        snapGrid_ctypes.argtypes = [
                ndpointer(ctypes.c_double,ndim=1,shape=(nData),flags='C_CONTIGUOUS'),
                ndpointer(ctypes.c_double,ndim=1,shape=(nData),flags='C_CONTIGUOUS'),
                ndpointer(ctypes.c_double,ndim=1,shape=(nData),flags='C_CONTIGUOUS'),
                ctypes.c_int64,
                ndpointer(ctypes.c_double,ndim=2,shape=(gridRows,gridCols),flags='C_CONTIGUOUS'),
                ndpointer(ctypes.c_double,ndim=2,shape=(gridRows,gridCols),flags='C_CONTIGUOUS'),
                ndpointer(ctypes.c_double,ndim=2,shape=(gridRows,gridCols),flags='C_CONTIGUOUS'),
                ndpointer(ctypes.c_int64,ndim=2,shape=(gridRows,gridCols),flags='C_CONTIGUOUS'),
                ctypes.c_int32,
                ctypes.c_int32
                ]

        t1 = time()

        '''
        int snapGrid_ctypes(double *lat, 
                        double *lon, 
                        double *data, 
                        long nData, 
                        double *gridLat,
                        double *gridLon,
                        double *gridData,
                        long *gridDataIdx,
                        int nGridRows,
                        int nGridCols
                        )
        '''

        retVal = snapGrid_ctypes(dataLat,
                                 dataLon,
                                 data,
                                 nData,
                                 gridLat,
                                 gridLon,
                                 gridData,
                                 gridDataIdx,
                                 gridRows,
                                 gridCols)

        return gridData,gridDataIdx

