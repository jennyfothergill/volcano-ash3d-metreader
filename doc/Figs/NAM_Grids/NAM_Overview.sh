#!/bin/bash
${GMTpre[GMTv]} gmtset PS_MEDIA letter

# We need to know if we must prefix all gmt commands with 'gmt', as required by version 5
GMTv=5
type gmt >/dev/null 2>&1 || { echo >&2 "Command 'gmt' not found.  Assuming GMTv4."; GMTv=4;}
GMTpre=("-" "-" "-" "-" " "   "gmt ")
GMTelp=("-" "-" "-" "-" "ELLIPSOID" "PROJ_ELLIPSOID")
GMTnan=("-" "-" "-" "-" "-Ts" "-Q")
GMTrgr=("-" "-" "-" "-" "grdreformat" "grdconvert")


FILE="NCEP-Grids"
CLON="250.0"
CLAT="50.0"

mapscale="6i"

PROJg=-JG${CLON}/${CLAT}/${mapscale}

AREA="-Rg"
DETAIL=-Dc
COAST="-G220/220/220 -W -Ggrey -Slightblue"

# Just coast
${GMTpre[GMTv]}  pscoast $AREA $PROJg $DETAIL $COAST -K > temp.ps

# Global data
#awk '{print $1, $2, 1.0}' n002_boundary.dat | ${GMTpre[GMTv]}  psxy $AREA $PROJg -P -K -O -W2 -V  >> temp.ps
#awk '{print $1, $2, 1.0}' n004_boundary.dat | ${GMTpre[GMTv]}  psxy $AREA $PROJg -P -K -O -W2 -V  >> temp.ps
#awk '{print $1, $2, 1.0}' nGCp_boundary.dat | ${GMTpre[GMTv]}  psxy $AREA $PROJg -P -K -O -W4 -V  >> temp.ps
#awk '{print $1, $2, 1.0}' nGNp_boundary.dat | ${GMTpre[GMTv]}  psxy $AREA $PROJg -P -K -O -W4 -V  >> temp.ps
#awk '{print $1, $2, 1.0}' nMCp_boundary.dat | ${GMTpre[GMTv]}  psxy $AREA $PROJg -P -K -O -W4 -V  >> temp.ps
#awk '{print $1, $2, 1.0}' n193_boundary.dat | ${GMTpre[GMTv]}  psxy $AREA $PROJg -P -K -O -W4 -V  >> temp.ps


awk '{print $1, $2, 1.0}' Bound/n181_boundary.dat | ${GMTpre[GMTv]}  psxy $AREA $PROJg -P -K -O -W4 -V  >> temp.ps
awk '{print $1, $2, 1.0}' Bound/n182_boundary.dat | ${GMTpre[GMTv]}  psxy $AREA $PROJg -P -K -O -W4 -V  >> temp.ps
awk '{print $1, $2, 1.0}' Bound/n194_boundary.dat | ${GMTpre[GMTv]}  psxy $AREA $PROJg -P -K -O -W4 -V  >> temp.ps
awk '{print $1, $2, 1.0}' Bound/n196_boundary.dat | ${GMTpre[GMTv]}  psxy $AREA $PROJg -P -K -O -W4 -V  >> temp.ps
awk '{print $1, $2, 1.0}' Bound/n198_boundary.dat | ${GMTpre[GMTv]}  psxy $AREA $PROJg -P -K -O -W4 -V  >> temp.ps
awk '{print $1, $2, 1.0}' Bound/n211_boundary.dat | ${GMTpre[GMTv]}  psxy $AREA $PROJg -P -K -O -W4 -V  >> temp.ps
awk '{print $1, $2, 1.0}' Bound/n221_boundary.dat | ${GMTpre[GMTv]}  psxy $AREA $PROJg -P -K -O -W4 -V  >> temp.ps
awk '{print $1, $2, 1.0}' Bound/n216_boundary.dat | ${GMTpre[GMTv]}  psxy $AREA $PROJg -P -K -O -W4 -V  >> temp.ps
awk '{print $1, $2, 1.0}' Bound/n243_boundary.dat | ${GMTpre[GMTv]}  psxy $AREA $PROJg -P -K -O -W4 -V  >> temp.ps


   # gridlines and close postscript file
${GMTpre[GMTv]}  psbasemap -B360g90:."${FILE}": $AREA $PROJg -O >> temp.ps

ps2epsi temp.ps
epstopdf temp.epsi

mv temp.pdf Overview_${FILE}.pdf
convert temp.ps -rotate 90 Overview_${FILE}.png
rm temp.ps temp.epsi
