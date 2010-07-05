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
        Takes no arguments, and does nothing... yet.
        """
        
        pass

    def __call__(self, lat, lon, data, gridLat, gridLon):
        """
        Takes as arguments the latitude and longitude spacings,
        and creates a SnapToGrid object containing numpy arrays 
        of the latitude and longitude grid definitions.
        """

        lat  = np.ravel(lat)
        lon  = np.ravel(lon)
        data = np.ravel(data)
        newData = griddata(lat,lon,data,gridLat,gridLon)
        return newData

    def granToGrid(self,dataGranule):
        """
        granToGrid(dataGranule)

        This class method takes as input a numpy array 
        containing data product that we desire to grid.

        """
        pass
