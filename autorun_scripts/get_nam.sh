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

NAM=$1
yearmonthday=$2
FChour=$3
#SERVER="https://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod"
SERVER="ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod"

echo "------------------------------------------------------------"
echo "running get_nam.sh script for ${NAM} $yearmonthday ${FChour}"
echo `date`
echo "------------------------------------------------------------"
t0=`date`

case ${NAM} in
 196)
  # HI 2.5 km
  HourMax=36
  HourStep=1
  #        nam.t00z.hawaiinest.hiresf00.tm00.grib2
  FilePre="nam.t${FChour}z.hawaiinest.hiresf"
  FilePost=".tm00.grib2"
  ;;
 091)
  # AK 2.95 km
  HourMax=36
  HourStep=1
  #        nam.t06z.alaskanest.hiresf00.tm00.grib2
  FilePre="nam.t${FChour}z.alaskanest.hiresf"
  FilePost=".tm00.grib2"
  ;;
 *)
  echo "NAM product not recognized"
  echo "Valid values: 091, 196"
  exit
esac

rc=0
WINDROOT="/bsuscratch/${USER}/WindFiles"
NAMDATAHOME="${WINDROOT}/nam/${NAM}"
install -d ${NAMDATAHOME}
if [[ $? -ne 0 ]] ; then
   echo "Error:  Download directory ${NAMDATAHOME} cannot be"
   echo "        created or has insufficient write permissions."
   rc=$((rc + 1))
   exit $rc
fi

#name of directory containing current files
FC_day=${yearmonthday}_${FChour}

#******************************************************************************
#START EXECUTING

#go to correct directory
cd $NAMDATAHOME
mkdir -p $FC_day
cd $FC_day

t=0
while [ "$t" -le ${HourMax} ]; do
  if [ "$t" -le 9 ]; then
      hour="0$t"
   else
      hour="$t"
  fi
  INFILE=${FilePre}${hour}${FilePost}
  fileURL=${SERVER}/nam.${yearmonthday}/$INFILE
  time wget ${fileURL}
  /opt/USGS/bin/gen_GRIB2_index $INFILE
  t=$(($t+${HourStep}))
done

mkdir -p $NAMDATAHOME/latest
cd $NAMDATAHOME/latest
rm nam.*
ln -s ../$FC_day/* .

t1=`date`
echo "download start: $t0"
echo "download   end: $t1"
