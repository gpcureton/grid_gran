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

import numpy as np
from matplotlib.mlab import griddata
import scipy.weave as weave
from scipy.weave import converters

class SnapToGrid:

    """
    SnapToGrid


    """

    def __init__(self,interp='nn'):
        """
        __init__

        geoFile: filename of a representative geolocation file
        dataFile: filename of a representative data file

        Takes no arguments, and does nothing... yet.
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

    def writeToFile(self,gridFile='gridFile.h5'):
        """
        writeToFile

        gridFile: filename of a representative geolocation file

        """
        fileObj = pytables.openFile(gridFile,"w")
        # do some stuff...
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
                        if (fabs(gridData(gridLatPt,gridLonPt) - fillVal) < 0.001){
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

        gridDegen = np.zeros(np.shape(gridData),dtype=np.int32)
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
