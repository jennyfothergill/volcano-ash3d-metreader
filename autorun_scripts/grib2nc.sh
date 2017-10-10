#!/bin/bash

FILEROOT=$1
#GribFile="${FILEROOT}.grib2"
GribFile="${FILEROOT}"
NCFile="${FILEROOT}.nc"

java -Xmx2048m -classpath ~/ncj/netcdfAll-4.5.jar ucar.nc2.dataset.NetcdfDataset \
     -in ${GribFile} -out ${NCFile} -IsLargeFile
rm ${FILEROOT}.ncx2 ${FILEROOT}.gbx9

