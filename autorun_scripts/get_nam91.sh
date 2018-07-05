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

echo "get_nam91.sh:  checking input arguments"
if [ -z $5 ]
then
  echo "Error: Insufficient command-line arguments"
  echo "Usage:  get_nam91.sh YYYY MM DD FChour HH"
  echo "        where FChour = 00, 06, 12, or 18"
  echo "        and HH = 0,1,2,...36"
  exit 1
else
  YYYY=$1
  MM=$2
  DD=$3
  FChour=$4
  HH=$5
fi


yearmonthday=${YYYY}${MM}${DD}

echo "------------------------------------------------------------"
echo "running get_nam91.sh script for $yearmonthday ${FChour}"
echo `date`
echo "------------------------------------------------------------"
t0=`date`

WINDROOT="/data/WindFiles"
NAMDATAHOME="${WINDROOT}/nam/ak03km"

#name of directory containing current files
FC_day=${yearmonthday}_${FChour}
#FC_day="temp"

#******************************************************************************
#START EXECUTING

# Note: Since grid 91 files are so big (~800 Mb/hourly-file), we will follow the
#       instructions on
#       http://nomads.ncep.noaa.gov/txt_descriptions/fast_downloading_grib.shtml
#       and first download the idx file, process it, and grab only the grib
#       layers needed.

#go to correct directory
cd $NAMDATAHOME
mkdir $FC_day
cd $FC_day

t=0
SERVER=http://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod

# Variables as a function of isobaric2
#  Geopotential_height_isobaric
#  u-component_of_wind_isobaric
#  v-component_of_wind_isobaric
#  Vertical_velocity_pressure_isobaric
#  Temperature_isobaric
#  Relative_humidity_isobaric
#  Specific_humidity_isobaric
niso2=42
iso2=("10 mb" "20 mb" "30 mb" "50 mb" "75 mb" "100 mb" "125 mb" "150 mb" "175 mb" "200 mb" "225 mb" "250 mb" "275 mb" "300 mb" "325 mb" "350 mb" "375 mb" "400 mb" "425 mb" "450 mb" "475 mb" "500 mb" "525 mb" "550 mb" "575 mb" "600 mb" "625 mb" "650 mb" "675 mb" "700 mb" "725 mb" "750 mb" "775 mb" "800 mb" "825 mb" "850 mb" "875 mb" "900 mb" "925 mb" "950 mb" "975 mb" "1000 mb")
nvar_iso2=7
var_iso2=("HGT" "UGRD" "VGRD" "VVEL" "TMP" "RH" "SPFH")
# Variables as a function of isobaric3
#  Cloud_mixing_ratio_isobaric
#  Snow_mixing_ratio_isobaric
niso3=40
iso3=("30 mb" "50 mb" "75 mb" "100 mb" "125 mb" "150 mb" "175 mb" "200 mb" "225 mb" "250 mb" "275 mb" "300 mb" "325 mb" "350 mb" "375 mb" "400 mb" "425 mb" "450 mb" "475 mb" "500 mb" "525 mb" "550 mb" "575 mb" "600 mb" "625 mb" "650 mb" "675 mb" "700 mb" "725 mb" "750 mb" "775 mb" "800 mb" "825 mb" "850 mb" "875 mb" "900 mb" "925 mb" "950 mb" "975 mb" "1000 mb")
nvar_iso3=2
var_iso3=("CLWMR" "SNMR")
# Variables as a function of height_above_ground7
#  u-component_of_wind_height_above_ground
#  v-component_of_wind_height_above_ground
nhag7=2
hag7=("10 m above ground" "80 m above ground")
nvar_hag7=2
var_hag7=("UGRD" "VGRD")
# Variables as a function of depth_below_surface_layer
#  Volumetric_Soil_Moisture_Content_depth_below_surface_layer
ndbsl=4
dbsl=("0-0.1 m below ground" "0.1-0.4 m below ground" "0.4-1 m below ground" "1-2 m below ground")
nvar_dbsl=1
var_dbsl=("SOILW")
# Variables as a function of surface
#  Planetary_Boundary_Layer_Height_surface
#  Frictional_Velocity_surface
#  Snow_depth_surface
#  Surface_roughness_surface
#  Wind_speed_gust_surface
#  Categorical_Rain_surface
#  Categorical_Snow_surface
#  Categorical_Freezing_Rain_surface
#  Categorical_Ice_Pellets_surface
#  Precipitation_rate_surface
nsurf=1
surf=("surface")
nvar_surf=10
var_surf=("HPBL" "FRICV" "SNOD" "SFCR" "GUST" "CRAIN" "CSNOW" "CFRZR" "CICEP" "PRATE")
# Variables as a function of some special 2d variable
#  Pressure_cloud_base
#  Pressure_cloud_tops
#  Total_cloud_cover_entire_atmosphere
nmisc2d=3
misc2d=("cloud base" "cloud top" "entire atmosphere (considered as a single layer)")
nvar_misc2d=("PRES" "PRES" "TCDC")

# To download the whole forcast package, uncomment the while loop
#while [ "$t" -le 36 ]; do
t=${HH}

  if [ "$t" -le 9 ]; then
      hour="0$t"
   else
      hour="$t"
  fi
  INFILE=nam.t${FChour}z.alaskanest.hiresf${hour}.tm00.grib2
  INFILEx=${INFILE}.idx
  MyINFILE=nam.t${FChour}z.alaskanest.hiresf${hour}.tm00.avo.grib2
  URL=${SERVER}/nam.${yearmonthday}/${INFILE}
  echo "$URL"
  ~/bin/get_inv.pl $URL.idx > my_inv

  # Get all variables that are a function of isobaric2
  rm iso2.grib2
  touch iso2.grib2
  for (( iv=0;iv<$nvar_iso2;iv++))
  do
    for (( id=0;id<$niso2;id++))
    do
      echo "${iv} ${var_iso2[iv]} ${iso2[id]}"
      grep ":${var_iso2[iv]}:" my_inv | grep ":${iso2[id]}:" > rec.tmp
      cat rec.tmp | ~/bin/get_grib.pl $URL tmp.grib2
      cat iso2.grib2 tmp.grib2 >> iso2.grib2
    done
  done
  # Get all variables that are a function of isobaric3
  rm iso3.grib2
  touch iso3.grib2
  for (( iv=0;iv<$nvar_iso3;iv++))
  do
    for (( id=0;id<$niso3;id++))
    do
      echo "${iv} ${var_iso3[iv]} ${iso3[id]}"
      grep ":${var_iso3[iv]}:" my_inv | grep ":${iso3[id]}:" > rec.tmp
      cat rec.tmp | ~/bin/get_grib.pl $URL tmp.grib2
      cat iso3.grib2 tmp.grib2 >> iso3.grib2
    done
  done
  # Get all variables that are function of height_above_ground7
  rm hag7.grib2
  touch hag7.grib2
  for (( iv=0;iv<$nvar_hag7;iv++))
  do
    for (( id=0;id<$nhag7;id++))
    do
      echo "${iv} ${var_hag7[iv]} ${hag7[id]}"
      grep ":${var_hag7[iv]}:" my_inv | grep ":${hag7[id]}:" > rec.tmp
      cat rec.tmp | ~/bin/get_grib.pl $URL tmp.grib2
      cat hag7.grib2 tmp.grib2 >> hag7.grib2
    done
  done
  # Get all variables that are function of depth_below_surface_layer
  rm dbsl.grib2
  touch dbsl.grib2
  for (( iv=0;iv<$nvar_dbsl;iv++))
  do
    for (( id=0;id<$ndbsl;id++))
    do
      echo "${iv} ${var_dbsl[iv]} ${dbsl[id]}"
      grep ":${var_dbsl[iv]}:" my_inv | grep ":${dbsl[id]}:" > rec.tmp
      cat rec.tmp | ~/bin/get_grib.pl $URL tmp.grib2
      cat dbsl.grib2 tmp.grib2 >> dbsl.grib2
    done
  done
  # Get all variables that are function of surface
  rm surf.grib2
  touch surf.grib2
  for (( iv=0;iv<$nvar_surf;iv++))
  do
    for (( id=0;id<$nsurf;id++))
    do
      echo "${iv} ${var_surf[iv]} ${surf[id]}"
      grep ":${var_surf[iv]}:" my_inv | grep ":${surf[id]}:" > rec.tmp
      cat rec.tmp | ~/bin/get_grib.pl $URL tmp.grib2
      cat surf.grib2 tmp.grib2 >> surf.grib2
    done
  done
  # Get all variables that are function of some special 2d variable
  rm misc2d.grib2
  touch misc2d.grib2
  for (( iv=0;iv<$nvar_misc2d;iv++))
  do
    echo "${iv} ${var_misc2d[iv]} ${misc2d[iv]}"
    grep ":${var_misc2d[iv]}:" my_inv | grep ":${misc2d[iv]}:" > rec.tmp
    cat rec.tmp | ~/bin/get_grib.pl $URL tmp.grib2
    cat misc2d.grib2 tmp.grib2 >> misc2d.grib2
  done

  # Now bundle all these grib2 files into a grib2 file for this timestep
  cat iso2.grib2 iso3.grib2 hag7.grib2 dbsl.grib2 surf.grib2 misc2d.grib2 > ${MyINFILE}
  /opt/USGS/bin/gen_GRIB2_index ${MyINFILE}

  rm iso2.grib2 iso3.grib2 hag7.grib2 dbsl.grib2 surf.grib2 misc2d.grib2 rec.tmp tmp.grib2

#  t=$(($t+1))
#done

