#!/bin/sh

#####
#	USAGE STATEMENT
#####
function usage() {
    echo "$0 -f IFILE [-c CPT] [-s SCALE] [-j DX] [-a PSFILE] [-x X0,Y0,WID] [-t FONTSZ] [--shift_ss X,Y]"
    echo "    -f IFILE             Input file (fth fss fno [mag [val]])"
    echo "    -c CPT               Color palette file (color by val)"
    echo "    -s SCALE             Symbol scale factor (default: 0.002)"
    echo "    -j DX                Jitter symbols randomly by up to DX"
    echo "    -a PSFILE            Add ternary diagram to an existing PostScript file"
    echo "    -x X0,Y0,WID         Move or resize triangle(x0,y0,wid); origin at bottom left"
    echo "    -t FONTSZ            Define font size"
    echo "    --shift_ss X,Y       Move labels on ternary diagram (alternatively, _th, _no)"
    echo
    echo " NOTE: all units are inches"
    exit
}

#####
#	PARSE COMMAND LINE
#####
SCALE="0.002"
JITTER=0
APPEND=""
ALTER=""
FONTSZ=""
SS_SHFT="0,0"
TH_SHFT="0,0"
NO_SHFT="0,0"
if [ $# -eq 0 ]; then usage; fi
while [ "$1" != "" ]; do
    case $1 in
        -f)shift;IFILE=$1;;
        -s)shift;SCALE=$1;;
        -c)shift;CPT=$1;;
        -j)shift;JITTER=$1;;
        -a)shift;APPEND=$1;;
        -x)shift;ALTER=$1;;
        -t)shift;FONTSZ=$1;;
        --shift_ss)shift;SS_SHFT="$1";;
        --shift_th)shift;TH_SHFT="$1";;
        --shift_no)shift;NO_SHFT="$1";;
        -h)usage;;
         *)echo No option $1;usage
    esac
    shift
done
if [ -z $SCALE ];then echo Did not define scale, setting to 0.002;SCALE="0.002";fi

#####
#	CHECK FOR INPUT FILE AND FORMAT
#####
if [ -z "$IFILE" ]
then
    echo No input file defined
    exit
elif [ ! -f "$IFILE" ]
then
    echo No input file $IFILE found
    exit
fi
MODE=`sed -n "1p" $IFILE | awk '{print NF}'`
case $MODE in
    3)OPT="-Sc0.2i -Ggrey";;
    4)OPT="-Sci -Ggrey";;
    5)OPT="-Sci -C$CPT";;
    *)echo No option for this many fields in input file; exit;;
esac

#####
#	SET UP TRIANGLE COORDINATES
#####
D2R=`echo 3.14159265 180 | awk '{print $1/$2}'`
HEIGHT="1"
AREF=`echo 35.26 $D2R | awk '{print $1*$2}'`
XMIN=`echo $HEIGHT | awk '{print -$1/sqrt(3)}'`
XMAX=`echo $HEIGHT | awk '{print  $1/sqrt(3)}'`
YMIN=`echo $HEIGHT | awk '{print -$1/3}'`
YMAX=`echo $HEIGHT | awk '{print $1*2/3}'`

#####
#	PLOT TERNARY FIGURE
#####
if [ ! -z $ALTER ]
then
    X0=`echo $ALTER | awk -F, '{print $1}'`
    Y0=`echo $ALTER | awk -F, '{print $2}'`
    WID=`echo $ALTER | awk -F, '{print $3}'`
else
    X0="1.0"
    Y0="1.0"
    WID="5"
fi
SHFT="-Xa${X0}i -Ya${Y0}i"
WID=`echo $XMIN $XMAX $WID | awk '{print $3/($2-$1)}'`
PROJ="-Jx${WID}i $SHFT"
LIMS="-R$XMIN/$XMAX/$YMIN/$YMAX"
if [ -z "$APPEND" ]
then
    PSFILE="ternary.ps"
    echo 0 0 | gmt psxy $PROJ $LIMS -K > $PSFILE
else
    PSFILE="$APPEND"
fi

# Triangle background
gmt psxy $PROJ $LIMS -Gwhite -K -O >> $PSFILE << EOF
$XMIN $YMIN
$XMAX $YMIN
0 $YMAX
$XMIN $YMIN
EOF

# Gridlines
for TDIP in 10 20 30 40 50 60 70 80
do
    case $TDIP in
        50)PEN="1p,55/55/55,4_2:0";;
         *)PEN="0.5p,225/225/225,4_2:0";;
    esac
    TDIPR=`echo $TDIP $D2R | awk '{print $1*$2}'`
    FTH=`echo $TDIPR | awk '{print sin($1)*sin($1)}'`
    LFTOVR=`echo $FTH | awk '{x=sqrt(1-$1);print atan2(x,sqrt(1-x*x))}'`
    grid -x 0 $LFTOVR -dx 0.005 -p |\
        awk '{print $0}END{print '"$LFTOVR"'}' |\
        awk 'BEGIN{
            tdipr = '"$TDIPR"'
            fth = '"$FTH"'
        }{
            pdipr = $1
            fno = sin(pdipr)*sin(pdipr)
            fss = 1-fth-fno
            if (fss<0) {fss = -fss}
            bdipr = atan2(sqrt(fss),sqrt(1-fss))
            print tdipr,bdipr,pdipr
        }' |\
        awk 'BEGIN{
            aref = '"$AREF"'
            hgt = '"$HEIGHT"'
        }{
            tdipr = $1
            bdipr = $2
            pdipr = $3
            psi = atan2(sin(tdipr),sin(pdipr)) - 3.14159265/4
            denom = sin(aref)*sin(bdipr) + cos(aref)*cos(bdipr)*cos(psi)
            h = hgt*sqrt(2)*cos(bdipr)*sin(psi)/(3*denom)
            v = hgt*sqrt(2)*(cos(aref)*sin(bdipr)-sin(aref)*cos(bdipr)*cos(psi))/(3*denom)
            print h,v
        }' |\
        gmt psxy $PROJ $LIMS -W$PEN -K -O >> $PSFILE
done
for BDIP in 10 20 30 40 50 60 70 80
do
    case $BDIP in
        60)PEN="1p,55/55/55,4_2:0";;
         *)PEN="0.5p,225/225/225,4_2:0";;
    esac
    BDIPR=`echo $BDIP $D2R | awk '{print $1*$2}'`
    FSS=`echo $BDIPR | awk '{print sin($1)*sin($1)}'`
    LFTOVR=`echo $FSS | awk '{x=sqrt(1-$1);print atan2(x,sqrt(1-x*x))}'`
    grid -x 0 $LFTOVR -dx 0.005 -p |\
        awk '{print $0}END{print '"$LFTOVR"'}' |\
        awk 'BEGIN{
            bdipr = '"$BDIPR"'
            fss = '"$FSS"'
        }{
            pdipr = $1
            fno = sin(pdipr)*sin(pdipr)
            fth = 1-fss-fno
            if (fth<0) {fth = -fth}
            tdipr = atan2(sqrt(fth),sqrt(1-fth))
            print tdipr,bdipr,pdipr
        }' |\
        awk 'BEGIN{
            aref = '"$AREF"'
            hgt = '"$HEIGHT"'
        }{
            tdipr = $1
            bdipr = $2
            pdipr = $3
            psi = atan2(sin(tdipr),sin(pdipr)) - 3.14159265/4
            denom = sin(aref)*sin(bdipr) + cos(aref)*cos(bdipr)*cos(psi)
            h = hgt*sqrt(2)*cos(bdipr)*sin(psi)/(3*denom)
            v = hgt*sqrt(2)*(cos(aref)*sin(bdipr)-sin(aref)*cos(bdipr)*cos(psi))/(3*denom)
            print h,v
        }' |\
        gmt psxy $PROJ $LIMS -W$PEN -K -O >> $PSFILE
done
for PDIP in 10 20 30 40 50 60 70 80
do
    case $PDIP in
        60)PEN="1p,55/55/55,4_2:0";;
         *)PEN="0.5p,225/225/225,4_2:0";;
    esac
    PDIPR=`echo $PDIP $D2R | awk '{print $1*$2}'`
    FNO=`echo $PDIPR | awk '{print sin($1)*sin($1)}'`
    LFTOVR=`echo $FNO | awk '{x=sqrt(1-$1);print atan2(x,sqrt(1-x*x))}'`
    grid -x 0 $LFTOVR -dx 0.005 -p |\
        awk '{print $0}END{print '"$LFTOVR"'}' |\
        awk 'BEGIN{
            pdipr = '"$PDIPR"'
            fno = '"$FNO"'
        }{
            bdipr = $1
            fss = sin(bdipr)*sin(bdipr)
            fth = 1-fss-fno
            if (fth<0) {fth = -fth}
            tdipr = atan2(sqrt(fth),sqrt(1-fth))
            print tdipr,bdipr,pdipr
        }' |\
        awk 'BEGIN{
            aref = '"$AREF"'
            hgt = '"$HEIGHT"'
        }{
            tdipr = $1
            bdipr = $2
            pdipr = $3
            psi = atan2(sin(tdipr),sin(pdipr)) - 3.14159265/4
            denom = sin(aref)*sin(bdipr) + cos(aref)*cos(bdipr)*cos(psi)
            h = hgt*sqrt(2)*cos(bdipr)*sin(psi)/(3*denom)
            v = hgt*sqrt(2)*(cos(aref)*sin(bdipr)-sin(aref)*cos(bdipr)*cos(psi))/(3*denom)
            print h,v
        }' |\
        gmt psxy $PROJ $LIMS -W$PEN -K -O >> $PSFILE
done

# Triangle outline and labels
gmt psxy $PROJ $LIMS -W1p -K -O >> $PSFILE << EOF
$XMIN $YMIN
$XMAX $YMIN
0 $YMAX
$XMIN $YMIN
EOF
if [ -z "$FONTSZ" ]
then
    FONT_SIZE=`echo $WID 5 24 | awk '{print $1/$2*$3}'`
else
    FONT_SIZE="$FONTSZ"
fi
OFFSET=`echo $WID 5 0.2 | awk '{print $1/$2*$3}'`
SS_DX=`echo $SS_SHFT | awk -F, '{print $1}'`
SS_DY=`echo $SS_SHFT | awk -F, '{print $2+'"$OFFSET"'}'`
echo 0 $YMAX $FONT_SIZE,0 CB Slip | gmt pstext $PROJ $LIMS -F+f+j -D${SS_DX}i/${SS_DY}i -N -K -O >> $PSFILE
SS_DY=`echo $SS_DY $FONT_SIZE | awk '{print $1+$2/72}'`
echo 0 $YMAX $FONT_SIZE,0 CB Strike | gmt pstext $PROJ $LIMS -F+f+j -D${SS_DX}i/${SS_DY}i -N -K -O >> $PSFILE
NO_DX=`echo $NO_SHFT | awk -F, '{print $1}'`
NO_DY=`echo $NO_SHFT | awk -F, '{print $2-'"$OFFSET"'}'`
echo $XMIN $YMIN $FONT_SIZE,0 LT Normal | gmt pstext $PROJ $LIMS -F+f+j -D${NO_DX}i/${NO_DY}i -N -K -O >> $PSFILE
TH_DX=`echo $TH_SHFT | awk -F, '{print $1}'`
TH_DY=`echo $TH_SHFT | awk -F, '{print $2-'"$OFFSET"'}'`
echo $XMAX $YMIN $FONT_SIZE,0 RT Thrust | gmt pstext $PROJ $LIMS -F+f+j -D${TH_DX}i/${TH_DY}i -N -K -O >> $PSFILE

# Plot data
TIME=`date "+%s"`
awk 'BEGIN{
    aref = '"$AREF"'
    hgt  = '"$HEIGHT"'
    scl  = '"$SCALE"'
    jitter = '"$JITTER"'
    srand('"$TIME"')
}{
    fth = $1
    fss = $2
    fno = $3
    tdipr = atan2(sqrt(fth),sqrt(1-fth))
    bdipr = atan2(sqrt(fss),sqrt(1-fss))
    pdipr = atan2(sqrt(fno),sqrt(1-fno))
    psi = atan2(sin(tdipr),sin(pdipr)) - 3.14159265/4
    denom = sin(aref)*sin(bdipr) + cos(aref)*cos(bdipr)*cos(psi)
    h = hgt*sqrt(2)*cos(bdipr)*sin(psi)/(3*denom)
    v = hgt*sqrt(2)*(cos(aref)*sin(bdipr)-sin(aref)*cos(bdipr)*cos(psi))/(3*denom)
    if (jitter>0) {
        h = h + (rand()-0.5)*0.02
        v = v + (rand()-0.5)*0.02
    }
    if (NF==3) {
        print h,v
    } else if (NF==4) {
        print h,v,$4*$4*$4*scl
    } else if (NF==5) {
        print h,v,$5,$4*$4*$4*scl
    }
}' $IFILE |\
    gmt psxy $PROJ $LIMS $OPT -W0.5p -N -K -O >> $PSFILE

if [ -z $APPEND ]
then
    echo 0 0 | gmt psxy $PROJ $LIMS -O >> $PSFILE
    ps2pdf $PSFILE
    rm $PSFILE
fi



