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

# Shell script that looks for windfiles that need to be purged

echo "------------------------------------------------------------"
echo "running prune_windfiles.sh"
echo `date`
echo "------------------------------------------------------------"

GFS_retain=17
PUFF_retain=1
NAM196_retain=0
NAM091_retain=0
MetProf_retain=30
Hysplit_retain=7

WINDROOT="/data/WindFiles"
GFSDATAHOME="${WINDROOT}/gfs"
PUFFDATAHOME="${WINDROOT}/puff/gfs"
NAM196DATAHOME="${WINDROOT}/nam/196"
NAM091DATAHOME="${WINDROOT}/nam/091"
METPROFDATAHOME="${WINDROOT}/MetProfiles"
HYSPLITDATAHOME="${WINDROOT}/Hysplit_traj"

#******************************************************************************
#FIRST, DELETE OLD FILES
echo "Deleting all GFS windfiles older than ${GFS_retain} days"
find ${GFSDATAHOME} -type f -mtime +${GFS_retain} -exec rm '{}' \;
find ${GFSDATAHOME} -type d -empty -exec rmdir '{}' \;
echo "Deleting all Puff windfiles older than ${PUFF_retain} days"
find ${PUFFDATAHOME} -type f -mtime +${PUFF_retain} -exec rm '{}' \;
echo "Deleting all nam-HI windfiles older than ${NAM196_retain} days"
find ${NAM196DATAHOME} -type f -mtime +${NAM196_retain} -exec rm '{}' \;
find ${NAM196DATAHOME} -type d -empty -exec rmdir '{}' \;
echo "Deleting all nam-AK windfiles older than ${NAM091_retain} days"
find ${NAM091DATAHOME} -type f -mtime +${NAM091_retain} -exec rm '{}' \;
find ${NAM091DATAHOME} -type d -empty -exec rmdir '{}' \;

echo "Deleting all Met Sonde files older than ${MetProf_retain} days"
find ${METPROFDATAHOME} -type f -mtime +${MetProf_retain} -exec rm '{}' \;
#******************************************************************************

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo ""
echo `date`
echo "GFS directories remaining:"
echo `ls -l ${GFSDATAHOME}`
echo "------------------------------------------------------------"
echo "Puff files remaining:"
echo `ls -l ${PUFFDATAHOME}`
echo "------------------------------------------------------------"
echo "NAM196 directories remaining:"
echo `ls -l ${NAM196DATAHOME}`
echo "------------------------------------------------------------"
echo "NAM1091 directories remaining:"
echo `ls -l ${NAM091DATAHOME}`
echo "------------------------------------------------------------"

echo "------------------------------------------------------------"
echo "MetProf files remaining:"
echo `ls -l ${METPROFDATAHOME}`

echo "------------------------------------------------------------"
echo "finished purge_windfiles.sh"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
