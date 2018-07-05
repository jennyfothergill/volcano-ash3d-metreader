#!/bin/bash
FILE=$1
# We need to know if we must prefix all gmt commands with 'gmt', as required by version 5
GMTv=5
type gmt >/dev/null 2>&1 || { echo >&2 "Command 'gmt' not found.  Assuming GMTv4."; GMTv=4;}
GMTpre=("-" "-" "-" "-" " "   "gmt ")
GMTelp=("-" "-" "-" "-" "ELLIPSOID" "PROJ_ELLIPSOID")
GMTnan=("-" "-" "-" "-" "-Ts" "-Q")
GMTrgr=("-" "-" "-" "-" "grdreformat" "grdconvert")

#${GMTpre[GMTv]} gmtset PS_MEDIA letter

mapscale="6i"
PROJp=-JJ180.0/${mapscale} 

pltgrid="F"
#if [ "$FILE" = "n002" ]; then
#  DETAILp=-Dl
#  pltgrid="T"
#  lonw=0.0
#  lone=360.0
#el
if [ "$FILE" = "n003" ]; then
  DETAILp=-Dl
  pltgrid="T"
  lonw=0.0
  lone=360.0
#elif [ "$FILE" = "n004" ]; then
#  DETAILp=-Dl
#  pltgrid="F"
#  lonw=0.0
#  lone=360.0
#elif [ "$FILE" = "n193" ]; then
#  DETAILp=-Dl
#  pltgrid="F"
#  lonw=0.0
#  lone=360.0
elif [ "$FILE" = "nGCp" ]; then
  DETAILp=-Dl
  pltgrid="F"
  lonw=-180.0
  lone=180.0
#elif [ "$FILE" = "nGNp" ]; then
#  DETAILp=-Dl
#  pltgrid="F"
#  lonw=-180.0
#  lone=180.0
fi

lats=-90.0
latn=90.0
#echo "-R${lonw}/${lats}/${lone}/${latn}r"
AREAp="-R${lonw}/${lats}/${lone}/${latn}r"

DETAILp=-Dc
COASTp="-G220/220/220 -W"

ALLp="$AREAp $PROJp"

# Just coast
${GMTpre[GMTv]} pscoast $AREAp $PROJp $DETAILp $COASTp -K > temp.ps
if [ "$pltgrid" = "T" ]; then
  awk '{print $1, $2, 1.0}' Grid/${FILE}_grid.dat | ${GMTpre[GMTv]} psxy $AREAp $PROJp -K -O  -Sc0.1p -Gblack -W0.1  >> temp.ps
fi

   # gridlines and close postscript file
${GMTpre[GMTv]} psbasemap -B0g0 $ALLp -O >> temp.ps

ps2epsi temp.ps 
epstopdf temp.epsi
mv temp.pdf ${FILE}.pdf
convert ${FILE}.pdf -rotate 90 -background white -flatten ${FILE}.png
rm temp.ps temp.epsi
