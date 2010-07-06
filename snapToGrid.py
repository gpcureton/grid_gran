import numpy as np
from matplotlib.mlab import griddata

class SnapToGrid:

    """
    CLASS SnapToGrid

    This python class implements a gridding prodedure developed
    by Nick Bearson, and modified by Nadia Smith, both from SSEC.

    Upon instantiation the class object contains the definition
    of the required grid. Subsequent passage of the 

    Geoff Cureton, SSEC
    April 2010
    """

    def __init__(self):
        """
        __init__

        Takes no arguments, and does nothing... yet.
        """
        
        pass

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
        newData = griddata(lat,lon,data,gridLat,gridLon)
        return newData

    @staticmethod
    def regrid(lat, lon, data, gridLat, gridLon):
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
        newData = griddata(lat,lon,data,gridLat,gridLon)
        return newData

