#!/bin/bash

#####
#       USAGE STATEMENT
#####
function usage() {
    echo "Usage: trg_schem.sh -i STR,DIP,RAK [-a PSFILE] [-x OPT]" 1>&2
    echo 1>&2
    echo "-i STR,DIP,RAK   Input strike, dip, and rake" 1>&2
    echo "-a PSFILE        Add schematic to an existing PostScript file" 1>&2
    echo "-x OPT           Move/resize schematic (x0,y0,wid); origin at bottom left" 1>&2
    echo 1>&2
    exit 1
}
if [ $# -eq 0 ]; then usage; fi

#####
#	PARSE COMMAND LINE
#####
TSTR=""
TDIP=""
TRAK=""
KINEMATICS=""
APPEND=""
ALTER=""
while [ "$1" != "" ]
do
    case $1 in
        -i)shift;KINEMATICS="$1";;
        -a)shift;APPEND="$1";;
        -x)shift;ALTER="$1";;
         *)echo No option $1;usage;;
    esac
    shift
done

# Fault kinematics: str,dip,rak
if [ "$KINEMATICS" == "" ]
then
    echo "$0: define fault kinematics" 1>&2
    usage
fi
TSTR=`echo $KINEMATICS | awk -F, '{print $1}'`
TDIP=`echo $KINEMATICS | awk -F, '{print $2}'`
TRAK=`echo $KINEMATICS | awk -F, '{print $3}'`
if [ "$TSTR" == "" -o "$TDIP" == "" -o "$TRAK" == "" ]
then
    echo "$0: failed to parse $KINEMATICS into strike,dip,rake" 1>&2
    usage
fi

# Dimensions and location of figure
if [ -z $ALTER ]
then
    X0="0.5"
    Y0="0.5"
    WID="1.0"
else
    X0=`echo $ALTER | awk -F, '{print $1}'`
    Y0=`echo $ALTER | awk -F, '{print $2}'`
    WID=`echo $ALTER | awk -F, '{print $3}'`
fi
if [ "$X0" == "" -o "$Y0" == "" -o "$WID" == "" ]
then
    echo "$0: failed to parse $ALTER into x0,y0,wid" 1>&2
    usage
fi
SHIFT="-Xa${X0}i -Ya${Y0}i"
HEIGHT=`echo $WID | awk '{print $1*1.5}'`
PROJ="-JX${WID}i/${HEIGHT}i"
LIMS="-R0/${WID}/0/${HEIGHT}"

if [ -z $APPEND ]
then
    PSFILE="trg_schem.ps"
    echo 0 0 | gmt psxy $PROJ $LIMS $SHIFT -K > $PSFILE
else
    PSFILE="$APPEND"
fi

#####
#	PLOT SCHEMATIC
#####
gmt psxy $PROJ $LIMS -W1p -Gwhite $SHIFT -K -O >> $PSFILE << EOF
0.0 0.0
$WID 0.0
$WID $HEIGHT
0.0 $HEIGHT
0.0 0.0
EOF

WID2=`echo $WID 0.5 | awk '{print $1*$2}'`
X1=`echo $WID 0.1 | awk '{print $1*$2}'`
H1=`echo $HEIGHT 1.78 2 | awk '{print $1*$2/$3}'`
H2=`echo $HEIGHT 1.57 2 | awk '{print $1*$2/$3}'`
HSTR=`echo $HEIGHT 1.35 2 | awk '{print $1*$2/$3}'`
HDIP=`echo $HEIGHT 1.15 2 | awk '{print $1*$2/$3}'`
HRAK=`echo $HEIGHT 0.95 2 | awk '{print $1*$2/$3}'`
FONT1=`echo 11 $WID 1 | awk '{print $1*$2/$3}'`
FONT2=`echo 10 $WID 1 | awk '{print $1*$2/$3}'`
gmt pstext $PROJ $LIMS -F+f+j $SHIFT -N -K -O >> $PSFILE << EOF
$WID2 $H1 $FONT1,0 CB Target Fault
$WID2 $H2 $FONT1,0 CB Kinematics
$X1 $HSTR $FONT2,2 LB Strike: $TSTR\260
$X1 $HDIP $FONT2,2 LB Dip: $TDIP\260
$X1 $HRAK $FONT2,2 LB Rake: $TRAK\260
EOF

XF="$WID2"  # center of fault in x direction
YF=`echo $HEIGHT 0.45 2 | awk '{print $1*$2/$3}'`  # center of fault in y direction
LEN=`echo $WID 0.8 | awk '{print $1*$2}'` # length of slip vector arrows
DW=`echo $WID 0.15 | awk '{print $1*$2}'`   # offset of slip vector arrows
FLT=`echo $WID 0.35 | awk '{print $1*$2}'` # side length of square fault
HEAD=`echo $WID 6 | awk '{print $1*$2}'`
TAIL=`echo $WID 1 | awk '{print $1*$2}'`
UPDIP=`echo $WID 2 | awk '{print $1*$2}'`
FAULT=`echo $WID 1 | awk '{print $1*$2}'`
SV=`echo $TSTR $TDIP $TRAK |\
    awk '{
      pi = 4*atan2(1,1)
      d2r = pi/180
      coss = cos((90-$1)*d2r)
      sins = sin((90-$1)*d2r)
      cosd = cos($2*d2r)
      sind = sin($2*d2r)
      cosr = cos($3*d2r)
      sinr = sin($3*d2r)
      x =  cosr*coss - sinr*cosd*sins
      y =  cosr*sins + sinr*cosd*coss
      z = sinr*sind
      print atan2(y,x)/d2r
    }'`
DX=`echo $TSTR $DW | awk '{print $2*sin(($1+90)*0.01745)}'`
DY=`echo $TSTR $DW | awk '{print $2*cos(($1+90)*0.01745)}'`
# Footwall slip vector
echo $DX $DY $SV $XF $YF $LEN |\
    awk '{print (-1)*$1+$4,(-1)*$2+$5,$3+180,$6}' |\
    gmt psxy $PROJ $LIMS -Sv${HEAD}p+e+jc+a40 -W${TAIL}p -Gblack $SHIFT -K -O >> $PSFILE
# Fault square projected onto horizontal surface
echo $XF $YF $TSTR $FLT `echo $TDIP $FLT | awk '{print $2*cos($1*0.01745)}'` |\
    gmt psxy $PROJ $LIMS -SJ -W${FAULT}p -Gwhite@20 $SHIFT -K -O >> $PSFILE
# highlight updip edge of fault
echo $TSTR $TDIP $FLT $XF $YF |\
    awk '{
      d = $3*0.5
      print d*sin(($1-90)*0.01745)*cos($2*0.01745)+$4,
              d*cos(($1-90)*0.01745)*cos($2*0.01745)+$5,
                      $1,$3"i"}' |\
    gmt psxy $PROJ $LIMS -SV10p+jc -W${UPDIP}p,0/155/0 $SHIFT -K -O >> $PSFILE
# Hanging wall slip vector
echo $DX $DY $SV $XF $YF $LEN |\
    awk '{print $1+$4,$2+$5,$3,$6}' |\
    gmt psxy $PROJ $LIMS -Sv${HEAD}p+e+jc+a40 -W${TAIL}p -Gblack $SHIFT -K -O >> $PSFILE

if [ -z $APPEND ]
then
    echo 0 0 | gmt psxy $PROJ $LIMS -O >> $PSFILE
    ps2pdf $PSFILE
    rm $PSFILE
fi

exit
