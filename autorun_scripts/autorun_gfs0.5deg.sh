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

# Shell script that manages the download of the gfs 0.5 degree data files for the
# current date, and converts the file to NetCDF.
# This script expects a command line argument indicating which forecast package to download.
#   autorun_gfs0.5deg.sh 0   for the 00 foracast package

# Please edit the line below to be consistant with the install directory specified in
# the makefile
INSTALLDIR="/opt/USGS/"

FC=$1

if [[ "$FC" -eq 0 ]] ; then
 FChour="00"
 FChourR="0.0"
fi
if [[ "$FC" -eq 6 ]] ; then
 FChour="06"
 FChourR="6.0"
fi
if [[ "$FC" -eq 12 ]] ; then
 FChour="12"
 FChourR="12.0"
fi
if [[ "$FC" -eq 18 ]] ; then
 FChour="18"
 FChourR="18.0"
fi
if [[ "$FC" -eq 24 ]] ; then
 FChour="24"
 FChourR="24.0"
fi

yearmonthday=`date -u +%Y%m%d`

echo "------------------------------------------------------------"
echo "running autorun_gfs0.5deg ${yearmonthday} ${FChour} script"
echo "------------------------------------------------------------"

SCRIPTDIR="${INSTALLDIR}bin/autorun_scripts"

#script that gets the wind files
echo "  Calling ${SCRIPTDIR}/get_gfs0.5deg.sh ${yearmonthday} ${FChour}"
${SCRIPTDIR}/get_gfs0.5deg.sh ${yearmonthday} ${FChour}

#script that converts grib2 to netcdf
echo "  Calling ${SCRIPTDIR}/convert_gfs0.5deg.sh ${yearmonthday} ${FChour}"
${SCRIPTDIR}/convert_gfs0.5deg.sh ${yearmonthday} ${FChour}

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "finished autorun_gfs0.5deg_00 script"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

