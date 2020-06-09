#!/bin/bash

#      This file is a component of the volcanic ash transport and dispersion model Ash3d,
#      written at the U.S. Geological Survey by Hans F. Schwaiger (hschwaiger@usgs.gov),
#      Larry G. Mastin (lgmastin@usgs.gov), and Roger P. Denlinger (roger@usgs.gov).

#      The model and its source code are products of the U.S. Federal Government and therefore
#      bear no copyright.  They may be copied, redistributed and freely incorporated 
#      into derivative products.  However as a matter of scientific courtesy we ask that
#      you credit the authors and cite published documentation of this model (below) when
#      publishing or distributing derivative products.

#      Schwaiger, H.F., Denlinger, R.P., and Mastin, L.G., 2012, Ash3d, a finite-
#         volume, conservative numerical model for ash transport and tephra deposition,
#         Journal of Geophysical Research, 117, B04204, doi:10.1029/2011JB008968. 

#      We make no guarantees, expressed or implied, as to the usefulness of the software
#      and its documentation for any purpose.  We assume no responsibility to provide
#      technical support to users of this software.

# Shell script that converts a grib2 file to netcdf using netcdf-java
# This script takes one command-line argument, the name of the grib file

# Please edit these variables to match your system and location of netcdf-java
JAVAHOME="/usr/local/bin/"
NCJv="${HOME}/ncj/netcdfAll-4.5.jar"
#JAVAHOME="/usr/bin/"
#NCJv="${HOME}/ncj/netcdfAll-4.6.14.jar"

echo "------------------------------------------------------------"
echo "running grib2nc.sh $1"
echo `date`
echo "------------------------------------------------------------"

FILEROOT=$1
GribFile="${FILEROOT}"
NCFile="${FILEROOT}.nc"

${JAVAHOME}java -Xmx2048m -classpath ${NCJv} ucar.nc2.dataset.NetcdfDataset \
     -in ${GribFile} -out ${NCFile} -IsLargeFile

rm ${FILEROOT}.ncx2 ${FILEROOT}.gbx9

