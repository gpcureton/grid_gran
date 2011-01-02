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

        Takes no arguments, and does nothing... yet.
        """
        self.interp=interp
        # TODO : Make grid arrays data for the class instance.
        # TODO :     This way, we only read in the grids once,
        # TODO :     and can keep calling new data arrays for 
        # TODO :     same grid. 


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
    def snapGrid(lat, lon, data, gridLat, gridLon,saveGrid=True):
        """
        snapGrid

        This static class method takes as arguments the latitude 
        and longitude arrays, the data array which we wish to regrid, 
        and the lat and lon grids we are mapping to.

        Input...
            lat: 1D array of latitude values
            lon: 1D array of longitude values
            data: 1D array of data corresponding to lat and lon
            gridLat: 2D array defining new latitude grid
            gridLon: 2D array defining new longitude grid
            saveGrid: Boolean switch for whether new grid indicies
                  are saved and returned.

        Returns...
            newLat: 1D array containing gridded latitude for each value in data
            newLon: 1D array containing gridded longitude for each value in data
            latIdx: 1D array containing indicies of gridLat which data values were 
                    mapped to
            lonIdx: 1D array containing indicies of gridLon which data values were 
                    mapped to
            covered by original data are masked.

        Example code....

            N = 10000
            dataLat = 30.*np.random.rand(N)
            dataLon = 30.*np.random.rand(N)
            data = 0.5*(dataLat + dataLon)

            ppl.figure()
            ppl.scatter(dataLon,dataLat,s=1.,c=data,faceted=False)
            ppl.colorbar()
            ppl.show()

            degInc = 0.2
            grids = np.mgrid[0.:30.+degInc:degInc,0.:30.:degInc]
            gridLat,gridLon = grids[0],grids[1]

            reload(utils.snapToGrid)
            gridData = utils.snapToGrid.SnapToGrid.snapGrid(dataLat, dataLon, data, gridLat, gridLon)

            ppl.figure()
            ppl.scatter(gridLon,gridLat,s=2.5,c=gridData,faceted=False)
            ppl.colorbar()
            ppl.show()

            gridDataMasked = ma.masked_less(gridData,0.)

            ppl.figure()
            ppl.scatter(gridLon,gridLat,s=2.5,c=gridDataMasked,faceted=False)
            ppl.colorbar()
            ppl.show()

        """
        # Array for the gridded data, initialised with -999.9
        gridData = np.ones(np.shape(gridLat))*-999.9

        gridLatInc = np.abs(gridLat[1,0]-gridLat[0,0])
        gridLonInc = np.abs(gridLon[0,1]-gridLon[0,0])

        print "gridLatInc = ",gridLatInc
        print "gridLonInc = ",gridLonInc
        print "shape(gridData) = ",np.shape(gridData)

        #weave.inline(codeSnapGrid,
            #arg_names=['EBBT_indices','L_to_EBBT_tp','L_to_EBBT_rad','bandIndex','npts','radiance','bTemp'],
            #include_dirs=self.include_dirs,
            #headers=self.headers,
            #force=0)

        codeGetBracket = """
            int latGridIdxLo;
            int latGridIdxHi;
            int lonGridIdxLo;
            int lonGridIdxHi;

            latGridIdxLo = 1;
            latGridIdxHi = 2;
            lonGridIdxLo = 3;
            lonGridIdxHi = 4;

            printf("latGridIdxLo = %d\\n",latGridIdxLo);
            printf("latGridIdxHi = %d\\n",latGridIdxHi);
            printf("lonGridIdxLo = %d\\n",lonGridIdxLo);
            printf("lonGridIdxHi = %d\\n",lonGridIdxHi);

            //latGridIdxLo = latVal/gridLatInc;
            //latGridIdxHi = latGridIdxLo + 1;
            //lonGridIdxLo = lonVal/gridLonInc;
            //lonGridIdxHi = lonGridIdxLo + 1;
        """

        for idx in range(np.size(lat)):
            latVal,lonVal,dataVal = lat[idx],lon[idx],data[idx]
            #print "latVal,lonVal,dataVal = ",latVal,lonVal,dataVal

            # Determine bracketing indicies
            latGridIdxLo = int(np.floor(latVal/gridLatInc))
            latGridIdxHi = latGridIdxLo + 1
            lonGridIdxLo = int(np.floor(lonVal/gridLonInc))
            lonGridIdxHi = lonGridIdxLo + 1

            #print "latGridIdxLo,latGridIdxHi = ",latGridIdxLo,latGridIdxHi
            #print "lonGridIdxLo,lonGridIdxHi = ",lonGridIdxLo,lonGridIdxHi

            #weave.inline(codeGetBracket,
                #arg_names=['latVal','lonVal','gridLatInc','gridLonInc'],
                #include_dirs=self.include_dirs,
                #headers=self.headers,
                #force=0)

            #print "(C) latGridIdxLo,latGridIdxHi = ",latGridIdxLo,latGridIdxHi
            #print "(C) lonGridIdxLo,lonGridIdxHi = ",lonGridIdxLo,lonGridIdxHi

            gridPoints = (
                [ latGridIdxLo, latGridIdxLo, latGridIdxHi, latGridIdxHi ],
                [ lonGridIdxLo, lonGridIdxHi, lonGridIdxLo, lonGridIdxHi ]
            )

            minDist = 1000.

            try :
                for points in zip(gridPoints[0],gridPoints[1]) :
                    dist = np.sqrt((latVal-gridLat[points])**2. + (lonVal-gridLon[points])**2.)
                    #print "points,gridLat[points],gridLon[points],minDist,dist = ",points,gridLat[points],gridLon[points],minDist,dist
                    if (dist < minDist) :
                        #print "dist < minDist !!!"
                        #snapLatIdx = latGridIdx
                        #snapLonIdx = lonGridIdx
                        snapPoints = points
                        minDist = dist
                    #gridData[snapLatIdx,snapLonIdx] = data[idx]
                    #print "snapPoints,gridLat[snapPoints],gridLon[snapPoints],minDist,dist = ",snapPoints,gridLat[snapPoints],gridLon[snapPoints],minDist,dist
                    #print ""

                gridData[snapPoints] = data[idx]
            except IndexError :
                pass

            #print ""
            #print ">>>>>>>>>>>>>>>>>>>>>>>>>"
        return gridData
