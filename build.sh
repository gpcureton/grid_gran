#!/bin/bash
# Commands to build the grid_gran dynamically linked library

#
# buildbucket_gcc64 using either docker or singularity
#

# gitlab.ssec.wisc.edu:5555/cspp/buildbucket/buildbucket_gcc64:latest
curdir=$PWD
echo $(which gcc)
export HDF4_DIR=/opt/hdf4
export HDF5_DIR=/opt/hdf5
export NETCDF_DIR=/opt/netcdf

rm -f *.o
gcc -c -fPIC -O3 griddingAndGranulation.c -o griddingAndGranulation.o
gcc -shared -O3 -Wl,-soname,libgriddingAndGranulation.so.1 -lm -o libgriddingAndGranulation.so.1.0.1 griddingAndGranulation.o

find ./ -name "libgriddingAndGranulation.so.1.0.1" -print
