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
        lat  = np.ravel(lat)
        lon  = np.ravel(lon)
        data = np.ravel(data)
        newData = griddata(lat,lon,data,gridLat,gridLon,interp=self.interp)
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
        lat  = np.ravel(lat)
        lon  = np.ravel(lon)
        data = np.ravel(data)
        newData = griddata(lat,lon,data,gridLat,gridLon,interp=interp)
        return newData

