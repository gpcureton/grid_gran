#!/usr/bin/env python
# encoding: utf-8
"""
netCDF2bin_imager.py

Purpose: Read band arrays from a netcdf4/cf created by pytroll/mpop, and write them to flat 
binary files.

Input:
    * netCDF file containing band arrays

Output:
    * flat binary files for each band

Details:
    * 

Preconditions:
    * netCDF4 python module

Optional:
    * 

Minimum commandline:

    python netCDF2bin_imager.py  --input_files=INPUTFILES

where...

    INPUTFILES: The fully qualified path to the input files. May be a directory or a file glob.


Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2014-03-25.
Copyright (c) 2014 University of Wisconsin Regents. All rights reserved.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

file_Date = '$Date: 2014-01-23 09:02:22 -0600 (Thu, 23 Jan 2014) $'
file_Revision = '$Revision: 1885 $'
file_Author = '$Author: kathys $'
file_HeadURL = '$HeadURL: https://svn.ssec.wisc.edu/repos/jpss_adl/trunk/scripts/edr/adl_viirs_edr.py $'
file_Id = '$Id: adl_viirs_edr.py 1885 2014-01-23 15:02:22Z kathys $'

__author__ = 'Geoff Cureton <geoff.cureton@ssec.wisc.edu>'
__version__ = '$Id: adl_viirs_edr.py 1885 2014-01-23 15:02:22Z kathys $'
__docformat__ = 'Epytext'

import os, sys, logging, traceback
from os import path,uname,environ
import string
import re
import uuid
from shutil import rmtree,copyfile
from glob import glob
from time import time
from datetime import datetime,timedelta

import numpy as np
from numpy import ma
import copy

import matplotlib
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure

matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# This must come *after* the backend is specified.
import matplotlib.pyplot as ppl

from mpl_toolkits.basemap import Basemap

from netCDF4 import Dataset
from netCDF4 import num2date


def listData(file_attrs,file_dicts):
    '''
    List various file and data attributes
    '''

    print "\n### File Attributes ###\n"
    for keys in file_attrs.keys():
        print "{} : {}".format(keys,file_attrs[keys])

    print "\n### File Dimensions ###\n"
    print "Keys = {}".format(file_dicts['dimensions'].keys())
    for keys in file_dicts['dimensions'].keys():
        print "{} : {}\n".format(keys,file_dicts['dimensions'][keys])

    print "\n### File Variables ###\n"
    var_attrs = ['_FillValue','scale_factor','add_offset','coordinates','long_name','resolution','units','ndim','size','dtype','shape','dimensions']
    print "Keys = {}".format(file_dicts['variables'].keys())
    for var_key in file_dicts['variables'].keys():
        print "\nvar = {}".format(var_key)
        var = file_dicts['variables'][var_key]
        for attr in var_attrs:
            try:
                var_attr = getattr(var,attr)
                print "{}.{} = {}".format(var_key,attr,var_attr)
            except:
                pass



def plotData(data,data_mask,units,pngName,stride=1,dpi=200):
    '''
    Plot the input dataset in swath projection
    '''

    # General Setup
    dataShape = data.shape
    aspectRatio = float(dataShape[0])/float(dataShape[1])

    figWidth,figHeight = 4.,aspectRatio*4.
    ax_rect = [0.05, 0.15, 0.90, 0.80  ] # [left,bottom,width,height]

    fig = Figure(figsize=(figWidth,figHeight))
    canvas = FigureCanvas(fig)

    ax = fig.add_axes(ax_rect)

    data = ma.array(data,mask=data_mask)

    im = ax.imshow(data,interpolation='nearest')

    if units == '%':
        txt = ax.set_title('Reflectance')
    elif units == 'K':
        txt = ax.set_title('Brightness Temperature')
    else:
        txt = ax.set_title('scaled data')


    ppl.setp(ax.get_xticklines(),visible=False)
    ppl.setp(ax.get_yticklines(),visible=False)
    ppl.setp(ax.get_xticklabels(), visible=False)
    ppl.setp(ax.get_yticklabels(), visible=False)

    ax.set_aspect('equal')

    cax_rect = [0.05 , 0.05, 0.90 , 0.05 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')

    if units == '%':
        txt = cax.set_title('Reflectance (%)')
    elif units == 'K':
        txt = cax.set_title('Brightness Temperature (K)')
    else:
        txt = cax.set_title('scaled units')

    # Redraw the figure
    canvas.draw()
    canvas.print_figure(pngName,dpi=dpi)


def plotMapData(lat,lon,data,data_mask,units,pngName,stride=1,pointSize=1.,scatterPlot=False, \
        lat_0=None,lon_0=None,latMin=None,lonMin=None,latMax=None,lonMax=None,dpi=200):
    '''
    Plot the input dataset in mapped to particular projection
    '''

    # General Setup
    figWidth,figHeight = 8.,8.
    ax_rect = [0.05, 0.15, 0.90, 0.80  ] # [left,bottom,width,height]

    fig = Figure(figsize=(figWidth,figHeight))
    canvas = FigureCanvas(fig)

    ax = fig.add_axes(ax_rect)

    lat_0,lon_0 = 0.,0.

    #print "lonMin = {} degrees".format(lonMin)
    #print "lonMax = {} degrees".format(lonMax)
    #print "latMin = {} degrees".format(latMin)
    #print "latMax = {} degrees".format(latMax)

    m = Basemap(projection='cyl',lon_0=lon_0,lat_0=lat_0,lat_ts=lat_0,ax=ax,
            llcrnrlon=lonMin,
            llcrnrlat=latMin,
            urcrnrlon=lonMax,
            urcrnrlat=latMax
            )

    x,y=m(lon[::stride,::stride],lat[::stride,::stride])

    m.drawcoastlines()
    #m.drawstates()
    m.drawcountries()
    m.fillcontinents(color='0.85',zorder=0)
    m.drawparallels(np.arange( -90, 91,30), color = '0.25', linewidth = 0.5,labels=[1,0,1,0])
    m.drawmeridians(np.arange(-180,180,30), color = '0.25', linewidth = 0.5,labels=[1,0,1,0])

    data = ma.array(data[::stride,::stride],mask=data_mask[::stride,::stride])

    if scatterPlot:
        cs = m.scatter(x,y,s=pointSize,c=data,axes=ax,edgecolors='none')
    else:
        cs = m.pcolor(x,y,data,axes=ax,edgecolors='none',antialiased=False)

    if units == '%':
        txt = ax.set_title('Reflectance')
    elif units == 'K':
        txt = ax.set_title('Brightness Temperature')
    else:
        txt = ax.set_title('scaled data')

    ppl.setp(ax.get_xticklines(),visible=False)
    ppl.setp(ax.get_yticklines(),visible=False)
    ppl.setp(ax.get_xticklabels(), visible=False)
    ppl.setp(ax.get_yticklabels(), visible=False)

    #ax.set_aspect('equal')

    cax_rect = [0.05 , 0.05, 0.90 , 0.05 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes
    #cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    cb = fig.colorbar(cs, cax=cax, orientation='horizontal')

    if units == '%':
        txt = cax.set_title('Reflectance (%)')
    elif units == 'K':
        txt = cax.set_title('Brightness Temperature (K)')
    else:
        txt = cax.set_title('scaled units')

    # Redraw the figure
    canvas.draw()
    canvas.print_figure(pngName,dpi=dpi)


class dataset():
    '''
    Define a "dataset" datatpe, which contains useful attributes for a particular
    dataset.
    '''
    def __init__(self):
        self.bandNum
        self.units = None
        self.dpi = 200
        self.arrayName = None
        self.mapRes = None


def _argparse():
    '''
    Method to encapsulate the option parsing and various setup tasks.
    '''

    import argparse

    endianChoices = ['little','big']

    dataChoices=['REF1','REF2','REF3','BT1','BT2','REF','BT','ALL']

    mapRes = ['c','l','i']

    defaults = {
                'input_file':None,
                'dataset':'ALL',
                'stride':1,
                'unscaled':False,
                'scatter_plot':False,
                'dpi':200,
                'list_attrs':False,
                'pointSize':1,
                'make_plots':False,
                'map_plot':False,
                'no_binaries':False
                }

    description = '''Read band arrays from a netcdf4/cf created by pytroll/mpop, and write them to flat 
binary files.'''
    usage = "usage: %prog [mandatory args] [options]"
    version = __version__

    parser = argparse.ArgumentParser()

    # Mandatory arguments

    parser.add_argument('-i','--input_file',
                      action='store',
                      dest='input_file',
                      type=str,
                      required=True,
                      help="The fully qualified path to the input file."
                      )

    # Optional arguments 

    parser.add_argument('--dset',
                      action="store",
                      dest="dataset",
                      default=defaults["dataset"],
                      type=str,
                      choices=dataChoices,
                      help='''The imager dataset(s) toVIIRS algorithm to run.\n\n
                              Possible values are...
                              {}.
                           '''.format(dataChoices.__str__()[1:-1])
                      )

    parser.add_argument('--unscaled',
                      action="store_true",
                      dest="doUnscale",
                      default=defaults["unscaled"],
                      help="Attempt to unscale the band arrays before output to flat binary files. [default: {}]".format(defaults["unscaled"])
                      )
    
    parser.add_argument('--make_plots',
                      action="store_true",
                      dest="doPlots",
                      default=defaults["make_plots"],
                      help="Plot the band arrays to png files. [default: {}]".format(defaults["make_plots"])
                      )

    parser.add_argument('--map_plot',
                      action="store_true",
                      dest="doMapPlots",
                      default=defaults["map_plot"],
                      help="Plot the band arrays to png files, using a map projection (has no effect if --make_plots is not set). [default: {}]".format(defaults["map_plot"])
                      )

    parser.add_argument('--no_binaries',
                      action="store_true",
                      dest="noBinaries",
                      default=defaults["no_binaries"],
                      help="Plot the band arrays to png files. [default: {}]".format(defaults["make_plots"])
                      )

    parser.add_argument('--plotMin',
                      action="store",
                      dest="plotMin",
                      type=float,
                      help="Minimum value to plot."
                      )

    parser.add_argument('--plotMax',
                      action="store",
                      dest="plotMax",
                      type=float,
                      help="Maximum value to plot."
                      )

    parser.add_argument('-d','--dpi',
                      action="store",
                      dest="dpi",
                      default='200.',
                      type=float,
                      help="The resolution in dots per inch of the output png file. [default: {}]".format(defaults["dpi"])
                      )

    parser.add_argument('-S','--stride',
                      action="store",
                      dest="stride",
                      default=defaults["stride"],
                      type=int,
                      help="Sample every STRIDE pixels in the band data. [default: {}]".format(defaults["stride"])
                      )

    parser.add_argument('--scatter_plot',
                      action="store_true",
                      dest="doScatterPlot",
                      default=defaults["scatter_plot"],
                      help="Generate the plot using a scatterplot approach."
                      )

    parser.add_argument('--list_attrs',
                      action="store_true",
                      dest="doListAttrs",
                      default=defaults["list_attrs"],
                      help="List some file and dataset attributes."
                      )

    parser.add_argument('-P','--pointSize',
                      action="store",
                      dest="pointSize",
                      default=defaults["pointSize"],
                      type=float,
                      help="Size of the plot point used to represent each pixel. [default: {}]".format(defaults["pointSize"])
                      )

    parser.add_argument('--lat_0',
                      action="store",
                      dest="lat_0",
                      type=float,
                      help="Center latitude of plot."
                      )

    parser.add_argument('--lon_0',
                      action="store",
                      dest="lon_0",
                      type=float,
                      help="Center longitude of plot."
                      )

    parser.add_argument('--latMin',
                      action="store",
                      dest="latMin",
                      type=float,
                      help="Minimum latitude to plot."
                      )

    parser.add_argument('--latMax',
                      action="store",
                      dest="latMax",
                      type=float,
                      help="Maximum latitude to plot."
                      )

    parser.add_argument('--lonMin',
                      action="store",
                      dest="lonMin",
                      type=float,
                      help="Minimum longitude to plot."
                      )

    parser.add_argument('--lonMax',
                      action="store",
                      dest="lonMax",
                      type=float,
                      help="Maximum longitude to plot."
                      )

    args = parser.parse_args()

    return args

def main():

    options = _argparse()

    stride = options.stride
    lat_0  = options.lat_0
    lon_0  = options.lon_0
    latMin = options.latMin
    lonMin = options.lonMin
    latMax = options.latMax
    lonMax = options.lonMax
    doScatterPlot = options.doScatterPlot
    pointSize = options.pointSize
    dpi = options.dpi

    # Create the dataset dictionaries
    #dataDict = {}
    #if options.dataset == "REF":

    #if options.dataset == "BT":
    #if options.dataset == "ALL":

    dataDict = {
            'REF1',
            'REF2',
            'REF3',
            'BT1',
            'BT2',
            'REF',
            'BT',
            'ALL'
            }

    # Open the file object
    print "Opening {} ...".format(options.input_file)
    fileObj = Dataset(options.input_file)

    file_attrs = {}
    file_attrs['history'] = fileObj.history
    file_attrs['orbit'] = fileObj.orbit
    file_attrs['platform'] = fileObj.platform

    file_dicts = {}
    file_dicts['dimensions'] = fileObj.dimensions
    file_dicts['variables'] = fileObj.variables

    # List various file and data attributes
    if options.doListAttrs:
        listData(file_attrs,file_dicts)

    if not (options.noBinaries and not options.doPlots):

        print "We are reading the file Datasets"

        # Output the lats and lons

        lats = file_dicts['variables']['lat0'][:,:].astype(np.float32)
        lons = file_dicts['variables']['lon0'][:,:].astype(np.float32)

        if not options.noBinaries:
            geo_fileName = string.replace(options.input_file,'.nc','_latitude.bin')
            print "Writing latitude to {}...".format(geo_fileName)
            lats.tofile(geo_fileName)

            geo_fileName = string.replace(options.input_file,'.nc','_longitude.bin')
            print "Writing longitude to {}...".format(geo_fileName)
            lons.tofile(geo_fileName)

        # Output the bands

        for band_var in ['Image0','Image1']:
        #for band_var in ['Image0']:
            nBands = file_dicts['variables'][band_var].shape[2]
            dtype = file_dicts['variables'][band_var].dtype.name
            units = ''

            for band in range(nBands):
            #for band in [0]:
                var = file_dicts['variables'][band_var][:,:,band]
                fillValue = file_dicts['variables'][band_var]._FillValue
                print "fill value is {}".format(fillValue)
                if band_var == 'Image1':
                    var_mask = ma.masked_inside(var,32720,32768).mask
                else:
                    var_mask = ma.masked_equal(var,fillValue).mask

                scale_factor = file_dicts['variables'][band_var].scale_factor[band]
                offset =  file_dicts['variables'][band_var].add_offset[band]
                print "{} band {} (scale,offset) = ({},{})".format(band_var,band,scale_factor,offset)

                if options.doUnscale :
                    dtype = 'float32'
                    units = file_dicts['variables'][band_var].units
                    var = var * scale_factor + offset
                    var = ma.array(var,mask=var_mask,fill_value=-999.9).filled()

                print "{} band {} fill value = {}".format(band_var,band,np.min(var))
                print "{} band {} min value = {}".format(band_var,band,np.min(ma.array(var,mask=var_mask)))
                print "{} band {} max value = {}".format(band_var,band,np.max(ma.array(var,mask=var_mask)))

                bin_fileName = string.replace(options.input_file,'.nc','')
                bin_fileName = "{}_{}_band{}_{}.bin".format(bin_fileName,band_var,band,dtype)

                if not options.noBinaries:
                    print "Writing {} band {} to {}...".format(band_var,band,bin_fileName)
                    var.astype(dtype).tofile(bin_fileName)

                if options.doPlots:
                    png_fileName = string.replace(bin_fileName,'_{}.bin'.format(dtype),'.png')
                    print "Plotting {} band {} to {}...".format(band_var,band,png_fileName)
                    if options.doMapPlots:
                        plotMapData(lats,lons,var,var_mask,units,png_fileName,stride=stride,lat_0=lat_0,lon_0=lon_0,latMin=latMin,lonMin=lonMin,latMax=latMax,lonMax=lonMax,scatterPlot=doScatterPlot,pointSize=pointSize,dpi=dpi)
                    else:
                        plotData(var,var_mask,units,png_fileName,stride=stride,dpi=dpi)


    print "Closing {} ...".format(options.input_file)
    fileObj.close()

    sys.exit(0)


if __name__=='__main__':
    sys.exit(main())  
