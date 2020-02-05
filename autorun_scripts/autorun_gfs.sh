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
#   autorun_gfs.sh 0p25 0   for the 0.25 degree 00 forecast package

# Please edit the line below to be consistant with the install directory specified in
# the makefile
INSTALLDIR="/opt/USGS"

if [ $# -eq 0 ]
  then
  echo "No arguments supplied"
  echo "Usage: autorun_gfs.sh Resolution FCpackage"
  echo "       where Resolution = 1p00, 0p50, or 0p25"
  echo "             FCpackage  = 0, 6, 12, 18 or 24"
  exit
fi

GFS=$1
FC=$2

case ${GFS} in
 0p25)
  echo "GFS 0.25 degree"
  ;;
 0p50)
  echo "GFS 0.50 degree"
  ;;
 1p00)
  echo "GFS 1.00 degree"
  ;;
 *)
  echo "GFS product not recognized"
  echo "Valid values: 0p25, 0p50, 1p00"
  exit
esac

case ${FC} in
 0)
  FChour="00"
  FChourR="0.0"
  ;;
 6)
  FChour="06"
  FChourR="6.0"
  ;;
 12)
  FChour="12"
  FChourR="12.0"
  ;;
 18)
  FChour="18"
  FChourR="18.0"
  ;;
 24)
  FChour="24"
  FChourR="24.0"
  ;;
 *)
  echo "GFS forecast package not recognized"
  echo "Valid values: 0, 6, 12, 18, 24"
  exit
esac

yearmonthday=`date -u +%Y%m%d`

echo "------------------------------------------------------------"
echo "running autorun_gfs ${GFS} ${yearmonthday} ${FChour} script"
echo "------------------------------------------------------------"

SCRIPTDIR="${INSTALLDIR}/bin/autorun_scripts"

#script that gets the wind files
echo "  Calling ${SCRIPTDIR}/get_gfs.sh ${GFS} ${yearmonthday} ${FChour}"
${SCRIPTDIR}/get_gfs.sh ${GFS} ${yearmonthday} ${FChour}

script that converts grib2 to netcdf
echo "  Calling ${SCRIPTDIR}/convert_gfs.sh ${GFS} ${yearmonthday} ${FChour}"
${SCRIPTDIR}/convert_gfs.sh ${GFS} ${yearmonthday} ${FChour}

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "finished autorun_gfs script"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

