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
#if [ "$FILE" = "n005" ]; then
#  PROJp=-JS255.0/90/${mapscale}
#  DETAILp=-Dl
#  pltgrid="T"
#  xmin=-4953.03
#  xmax=4952.97
#  ymin=-9144.025
#  ymax=1523.974
#  projline="invproj +proj=stere  +lon_0=255  +lat_0=90 +k_0=0.933 +R=6371.229"
#elif [ "$FILE" = "n104" ]; then
#  PROJp=-JS255.0/90/${mapscale}
#  PROJp=-Js255.0/90/90/0.933
#  DETAILp=-Dl
#  xmin=-6761.21
#  xmax=6489.02
#  ymin=-9846.821
#  ymax=45.47379
#  projline="invproj +proj=stere  +lon_0=255  +lat_0=90 +k_0=0.933 +R=6371.229"
#elif [ "$FILE" = "n216" ]; then
if [ "$FILE" = "n216" ]; then
  PROJp=-JS225.0/90/${mapscale}
  DETAILp=-Dl
  pltgrid="T"
  xmin=-4225.928
  xmax=1984.072
  ymin=-5408.941
  ymax=-638.9415
  projline="invproj +proj=stere  +lon_0=225  +lat_0=90 +k_0=0.933 +R=6371.229"
#elif [ "$FILE" = "n242" ]; then  
#  PROJp=-JS225.0/90/${mapscale}
#  DETAILp=-Dl
#  pltgrid="F"
#  xmin=-4225.87071154953
#  xmax=1984.12928845047
#  ymin=-5408.86785597764
#  ymax=-638.867855977640
#  projline="invproj +proj=stere  +lon_0=225  +lat_0=90 +k_0=0.933 +R=6371.229"
elif [ "$FILE" = "n194" ]; then
  PROJp=-JM-75.5/15.0/${mapscale}
  DETAILp=-Dh
  pltgrid="F"
  xmin=0.0
  xmax=1357.5
  ymin=1585.61
  ymax=2358.11
  projline="invproj +proj=merc  +lat_ts=15.0 +lon_0=-75.5  +R=6371.229"
elif [ "$FILE" = "n196" ]; then
  PROJp=-JM198.475/20.0/${mapscale}
  DETAILp=-Df
  pltgrid="F"
  xmin=0.0
  xmax=800.0
  ymin=1920.618
  ymax=2480.618
  projline="invproj +proj=merc  +lat_ts=20.0 +lon_0=198.475  +R=6371.229"
#elif [ "$FILE" = "n198" ]; then  PROJp=-JS210.0/90/${mapscale}
#  PROJp=-JS210.0/90/${mapscale}
#  DETAILp=-Dl
#  pltgrid="F"
#  xmin=-2619.397
#  xmax=2285.875
#  ymin=-4810.103
#  ymax=-1524.047
#  projline="invproj +proj=stere  +lon_0=210  +lat_0=90 +k_0=0.933 +R=6371.229"
elif [ "$FILE" = "n091" ]; then  PROJp=-JS210.0/90/${mapscale}
  PROJp=-JS210.0/90/${mapscale}
  DETAILp=-Dl
  pltgrid="F"
  xmin=-2619.397
  xmax=2285.875
  ymin=-4810.103
  ymax=-1524.047
  projline="invproj +proj=stere  +lon_0=210  +lat_0=90 +k_0=0.933 +R=6371.229"
elif [ "$FILE" = "n211" ]; then
  PROJp=-JL265/25.0/25.0/25.0/${mapscale}
  DETAILp=-Dl
  pltgrid="T"
  xmin=-4226.108
  xmax=3250.823
  ymin=-832.6978
  ymax=4368.646
  projline="invproj +proj=lcc  +lon_0=265.0 +lat_0=25.0 +lat_1=25.0 +lat_2=25.0 +R=6371.229"
#elif [ "$FILE" = "n212" ]; then
#  PROJp=-JL265/25.0/25.0/25.0/${mapscale}
#  DETAILp=-Dl
#  pltgrid="T"
#  xmin=-4226.108
#  xmax=3250.731
#  ymin=-832.6978
#  ymax=4368.582
#  projline="invproj +proj=lcc  +lon_0=265.0 +lat_0=25.0 +lat_1=25.0 +lat_2=25.0 +R=6371.229"
#elif [ "$FILE" = "n218" ]; then
#  PROJp=-JL265/25.0/25.0/25.0/${mapscale}
#  DETAILp=-Dl
#  pltgrid="F"
#  xmin=-4226.108
#  xmax=3246.974
#  ymin=-832.6978
#  ymax=4372.859
#  projline="invproj +proj=lcc  +lon_0=265.0 +lat_0=25.0 +lat_1=25.0 +lat_2=25.0 +R=6371.229"
elif [ "$FILE" = "n221" ]; then
  PROJp=-JL-107/50.0/50.0/50.0/${mapscale}
  DETAILp=-Dl
  pltgrid="F"
  xmin=-5632.668
  xmax=5664.457
  ymin=-4612.566
  ymax=4347.222
  projline="invproj +proj=lcc  +lon_0=-107.0 +lat_0=50.0 +lat_1=50.0 +lat_2=50.0 +R=6371.229"
#elif [ "$FILE" = "n227" ]; then
#  PROJp=-JL265/25.0/25.0/25.0/${mapscale}
#  DETAILp=-Dl
#  pltgrid="F"
#  xmin=-4226.108
#  xmax=3250.179
#  ymin=-832.6978
#  ymax=4368.198
#  projline="invproj +proj=lcc  +lon_0=265.0 +lat_0=25.0 +lat_1=25.0 +lat_2=25.0 +R=6371.229"
fi

lonw=`echo "${xmin} ${ymin}" | ${projline} -f %.10f | cut -f1`
lats=`echo "${xmin} ${ymin}" | ${projline} -f %.10f | cut -f2`
lone=`echo "${xmax} ${ymax}" | ${projline} -f %.10f | cut -f1`
latn=`echo "${xmax} ${ymax}" | ${projline} -f %.10f | cut -f2`

#echo "-R${lonw}/${lats}/${lone}/${latn}r"
AREAp="-R${lonw}/${lats}/${lone}/${latn}r"

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

