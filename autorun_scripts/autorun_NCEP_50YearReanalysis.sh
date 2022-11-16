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

# Shell script that manages the download of the NCEP 2.5 degree Reanalysis data files.

# Please edit the line below to be consistant with the install directory specified in
# the makefile
INSTALLDIR="/cm/shared/apps/metreader/gcc/12.1.0/5f23d0b"
# This is the location where the downloaded windfiles will be placed.
# Please edit this to suit your system.
WINDROOT="/bsuscratch/${USER}/WindFiles"

echo "------------------------------------------------------------"
echo "running NCEP_50year_Reanalysis.sh"
echo `date`
echo "------------------------------------------------------------"

SCRIPTDIR="${INSTALLDIR}/bin/autorun_scripts"
echo "SCRIPTDIR=$SCRIPTDIR"

starttime=`date`                  #record when we're starting the download
echo "starting autorun_NCEP_50YearReanalysis.sh at $starttime"

rc=0
NCEPDATAHOME="${WINDROOT}/NCEP"
install -d ${NCEPDATAHOME}
if [[ $? -ne 0 ]]; then
   echo "Error:  Download directory ${NCEPDATAHOME} cannot be"
   echo "        created or has insufficient write permissions."
   rc=$((rc + 1))
   exit $rc
fi
echo "NCEPDATAHOME=$NCEPDATAHOME"

y=`date +%Y`
monthnow=`date +%m`
daynow=`date +%d`
echo "year=$y"

echo "monthnow=$monthnow, daynow=$daynow"
#If this is the beginning of the year, make sure we get the completed files from last year.
if [ "$monthnow" -eq 1 ] && [ "$daynow" -lt 8 ] ; then
   oldyear=`expr $y - 1`
   echo "getting the remainder of $oldyear"
   ${SCRIPTDIR}/get_NCEP_50YearReanalysis.sh $oldyear
   echo "done getting remainder of $oldyear"
fi

#if the directory for this year doesn't exist (e.g. it's Jan. 1), create it
echo "making sure the directory for year ${y} exists"
if [ ! -r "${NCEPDATAHOME}/${y}" ]         
then
   echo "It doesnt.  Creating directory for year ${y}"
   mkdir ${NCEPDATAHOME}/${y}
else
   echo "Good.  It does."
   echo "Copying contents of ${NCEPDATAHOME}/${y} to ${NCEPDATAHOME}/backup"
   cp -v ${NCEPDATAHOME}/${y}/* ${NCEPDATAHOME}/backup
   echo "all done copying files"
fi

if test -r ${NCEPDATAHOME}/dbuffer
then
    echo "moving to ${NCEPDATAHOME}/dbuffer"
    cd ${NCEPDATAHOME}/dbuffer
  else
    mkdir ${NCEPDATAHOME}/dbuffer
    echo "moving to ${NCEPDATAHOME}/dbuffer"
    cd ${NCEPDATAHOME}/dbuffer
fi

echo "Calling ${SCRIPTDIR}/get_NCEP_50YearReanalysis.sh $y"
${SCRIPTDIR}/get_NCEP_50YearReanalysis.sh $y

yearmonthday=`date -u +%Y%m%d`       #current year, month & day (e.g. 20110119)
echo ${yearmonthday} > ${NCEPDATAHOME}/last_downloaded.txt   #write date of last download to text file

endtime=`date`

echo "download started at $starttime"
echo "download ended at $endtime"

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "finished NCEP_50year_Reanalysis.sh"
echo `date`
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

