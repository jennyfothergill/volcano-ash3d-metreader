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

pltgrid="F"
if [ "$FILE" = "n181" ]; then
  PROJp=-JM-80.0/20.0/${mapscale} 
  DETAILp=-Dl
  pltgrid="F"
  lonw=-100.0
  lone=-60.148
  lats=0.1380001
  latn=30.054
elif [ "$FILE" = "n182" ]; then
  PROJp=-JM-155.0/20.0/${mapscale}
  DETAILp=-Dh
  pltgrid="F"
  lonw=-170.0
  lone=-140.084
  lats=8.133
  latn=32.973
elif [ "$FILE" = "n243" ]; then
  PROJp=-JM215.0/20.0/${mapscale}
  DETAILp=-Dl
  pltgrid="T"
  lonw=190.0
  lone=240.0
  lats=10.0
  latn=50.0
fi

#${GMTpre[GMTv]} 
#echo "-R${lonw}/${lats}/${lone}/${latn}r"
AREAp="-R${lonw}/${lats}/${lone}/${latn}r"

DETAILp=-Dh
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
