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
    def snapGrid_numpy(lat, lon, data, gridLat, gridLon,saveGrid=True):
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

        """
        # Array for the gridded data, initialised with -999.9
        gridData = np.ones(np.shape(gridLat))*-999.9

        gridLatInc = np.abs(gridLat[1,0]-gridLat[0,0])
        gridLonInc = np.abs(gridLon[0,1]-gridLon[0,0])

        print "gridLatInc = ",gridLatInc
        print "gridLonInc = ",gridLonInc
        print "shape(gridData) = ",np.shape(gridData)

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

                gridData[snapPoints] = dataVal
            except IndexError :
                pass

            #print ""
            #print ">>>>>>>>>>>>>>>>>>>>>>>>>"
        return gridData

    @staticmethod
    def snapGrid_weave(lat, lon, data, gridLat, gridLon,saveGrid=True):
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

        """
        # Array for the gridded data, initialised with -999.9
        gridData = np.ones(np.shape(gridLat))*-999.9

        gridLatInc = np.double(np.abs(gridLat[1,0]-gridLat[0,0]))
        gridLonInc = np.double(np.abs(gridLon[0,1]-gridLon[0,0]))

        print "gridLatInc = ",gridLatInc
        print "gridLonInc = ",gridLonInc
        print "shape(gridData) = ",np.shape(gridData)

        N_data = np.size(lat)

        codeSnapGrid = """


            printf("exp(2.71) = %f\\n",exp(2.71));
            printf("floor(2.71) = %f\\n",floor(2.71));
            printf("(int) floor(2.71) = %d\\n",(int) floor(2.71));

            long idx;
            long latGridIdxLo, latGridIdxHi, lonGridIdxLo, lonGridIdxHi;
            double minDist,dist,latVal,lonVal,dataVal;
            double gridLatInc,gridLonInc;

            //gridLatInc = abs(gridLat[1][0]-gridLat[0][0]);
            //gridLonInc = abs(gridLon[0][1]-gridLon[0][0]);
            //printf("gridLatInc = %lf\\n",(double) gridLatInc);
            //printf("gridLonInc = %lf\\n",(double) gridLonInc);

            printf("N_data = %d\\n",N_data);

            int rows,cols;
            for (rows=0;rows<5;rows++){
                for (cols=0;cols<5;cols++){
                    printf("gridLat[%d][%d] = %lf\\n",rows,cols,gridLat[rows,cols]);
                }
            }

            for (rows=0;rows<5;rows++){
                for (cols=0;cols<5;cols++){
                    printf("gridLon[%d][%d] = %lf\\n",rows,cols,gridLon[rows,cols]);
                }
            }


            idx = 0;
            for (idx=0;idx<1;idx++){
                latVal = lat[idx];
                lonVal = lon[idx];
                dataVal = data[idx];

                printf("\\nlatVal,lonVal,dataVal = %f, %f, %f\\n",latVal,lonVal,dataVal);
                printf("gridLatInc,gridLonInc = %lf,%lf\\n",(double)gridLatInc,(double)gridLonInc);

                latGridIdxLo = (long) floor(latVal/(double)gridLatInc);
                latGridIdxHi = latGridIdxLo + 1;
                lonGridIdxLo = (long) floor(lonVal/(double)gridLonInc);
                lonGridIdxHi = lonGridIdxLo + 1;

                long gridPoints[2][4] = {
                    {latGridIdxLo, latGridIdxLo, latGridIdxHi, latGridIdxHi},
                    {lonGridIdxLo, lonGridIdxHi, lonGridIdxLo, lonGridIdxHi}
                };

                minDist = 1000.;
                long gridLatPt,gridLonPt;
                double latDist,lonDist;
                int crnrPt,snapCrnrPt;

                for (crnrPt=0;crnrPt<4;crnrPt++){
                    printf("%ld, %ld\\n",gridPoints[0][crnrPt],gridPoints[1][crnrPt]);

                    gridLatPt = gridPoints[0][crnrPt];
                    gridLonPt = gridPoints[1][crnrPt];

                    printf("%ld, %ld\\n",gridLatPt,gridLonPt);

                    //latDist = latVal-gridData[gridLatPt][gridLonPt];
                    //lonDist = lonVal-gridLon[gridLatPt][gridLonPt];

                    //dist = sqrt(latDist*latDist + lonDist*lonDist); 
                    //printf("dist = %lf\\n",dist);

                    //if (dist < minDist){
                    //    snapCrnrPt = crnrPt;
                    //    minDist = dist;
                    //}
                }
                //gridData[gridPoints[0][crnrPt]][] = dataVal


            }


        """

        weave.inline(codeSnapGrid,
            arg_names=['N_data','lat','lon','data','gridLat','gridLon','gridData'],
            headers=['<math.h>'],
            libraries=['m'],
            #include_dirs=self.include_dirs,
            force=1)


            #print ""
            #print ">>>>>>>>>>>>>>>>>>>>>>>>>"
        return gridData
