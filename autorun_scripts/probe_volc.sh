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
  echo "Usage: probe_volc.sh product FCpackage"
  echo "       where Resolution = gfs1p00, gfs0p50, gfs0p25, or nam091"
  echo "             FCpackage  = 0, 6, 12, 18 or 24"
  exit
fi

PROD=$1
FC=$2

case ${PROD} in
 gfs0p25)
  echo "GFS 0.25 degree"
  ;;
 gfs0p50)
  echo "GFS 0.50 degree"
  HourMax=198
  HourStep=3
  ;;
 gfs1p00)
  echo "GFS 1.00 degree"
  ;;
 nam091)
  echo "NAM 2.95 km"
  HourMax=34
  #HourMax=1
  HourStep=1
  ;;
 *)
  echo "Product not recognized"
  echo "Valid values: gfs0p25, gfs0p50, gfs1p00, nam092"
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
  echo "Forecast package not recognized"
  echo "Valid values: 0, 6, 12, 18, 24"
  exit
esac

yearmonthday=`date -u +%Y%m%d`
# Here you can over-ride the date if need be
#yearmonthday="20200626"

echo "------------------------------------------------------------"
echo "running probe_volc ${PROD} ${yearmonthday} ${FChour} script"
echo "------------------------------------------------------------"

SCRIPTDIR="${INSTALLDIR}/bin/autorun_scripts"

SONDEDIR="/data/WindFiles/sonde"
volc=`cat ${SONDEDIR}/volc.dat | cut -d' ' -f1`
lon=`cat  ${SONDEDIR}/volc.dat | cut -d' ' -f2`
lat=`cat  ${SONDEDIR}/volc.dat | cut -d' ' -f3`
echo "$volc $lon $lat"
#mkdir -p ${SONDEDIR}/${volc}/${yearmonthday}
#cd ${SONDEDIR}/${volc}/${yearmonthday}
mkdir -p ${SONDEDIR}/${volc}/temp
cd ${SONDEDIR}/${volc}/temp

rm -f probe.log
touch probe.log
t=0
while [ "$t" -le ${HourMax} ]
do
  if [ "$t" -le 9 ]; then
      hour="00$t"
   elif [ "$t" -le 99 ]; then
      hour="0$t"
   else
      hour="$t"
  fi

  case ${PROD} in
   gfs0p50)
    ARGS="1 0 $lon $lat T 4 1 2 3 5 4 20 2"
    WINDPATH=/data/WindFiles/gfs/gfs.${yearmonthday}${FChour}
    WINDFILE=${yearmonthday}${FChour}.f${hour}.nc
    CLOUD="ncks -C -d lon,${lon} -d lat,${lat} -H -Q -s '%f ' -v Total_cloud_cover_isobaric ${WINDFILE} > Cloud.dat"
    ;;
   nam091)
    hour=`echo "${hour}" | cut -c2-3`
    ARGS="1 1 $lon $lat F 4 1 2 3 5 4 13 3"
    WINDPATH=/data/WindFiles/nam/091/${yearmonthday}_${FChour}
    WINDFILE=nam.t00z.alaskanest.hiresf${hour}.tm00.loc.grib2
    CLOUD="touch Cloud.dat"
    ;;
   *)
    echo "Product not recognized"
    echo "Valid values: gfs0p25, gfs0p50, gfs1p00, nam092"
    exit
  esac


  ln -s ${WINDPATH}/${WINDFILE} .
  echo "${INSTALLDIR}/bin/probe_Met ${WINDFILE} ${ARGS}"
  ${INSTALLDIR}/bin/probe_Met ${WINDFILE} ${ARGS}
  HourOffset=`echo "${FChour} + ${t}"  | bc`
  NewYYYYMMDD=`date -d"${yearmonthday} +${HourOffset} hour" -u +%Y%m%d`
  Newhour=`date -d"${yearmonthday} +${HourOffset} hour" -u +%H`
  mkdir -p ${SONDEDIR}/${volc}/${NewYYYYMMDD}
  mv NWP_prof.dat ${SONDEDIR}/${volc}/${NewYYYYMMDD}/${volc}_${PROD}_phuvt_${NewYYYYMMDD}_${Newhour}.dat
  #ncks -C -d lon,${lon} -d lat,${lat} -H -Q -s '%f ' -v Total_cloud_cover_isobaric ${WINDFILE} > Cloud.dat
  ${CLOUD}
  mv Cloud.dat ${SONDEDIR}/${volc}/${NewYYYYMMDD}/${volc}_${PROD}_cloud_${NewYYYYMMDD}_${Newhour}.dat
  rm ${WINDFILE}
  #echo "$t : Mapping ${WINDFILE} to ${volc}_gfs_phuvt_${NewYYYYMMDD}_${Newhour}.dat" >> probe.log
  t=$((t+${HourStep}))
done

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "finished probe_volc script"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

