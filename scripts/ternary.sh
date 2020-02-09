#!/bin/bash

#####
#	USAGE STATEMENT
#####
function usage() {
    echo "Usage: ternary.sh -f IFILE [...options...]" 1>&2
    echo 1>&2
    echo "-f IFILE             Input file (fth fss fno [mag [val]])" 1>&2
    echo "-c CPT               Color palette file (color by val in fifth column)" 1>&2
    echo "-s SCALE             Symbol scale factor (default: 0.2)" 1>&2
    echo "-j DX                Jitter symbols randomly by up to DX in (default: none)" 1>&2
    echo "-a PSFILE            Add ternary diagram to an existing PostScript file" 1>&2
    echo "-x X0,Y0,WID         Move/resize triangle; origin at bottom left (default: 1,1,5 in)" 1>&2
    echo "-t FONTSZ            Define font size (default: scaled with WID)" 1>&2
    echo "--shift_ss X,Y       Move labels on ternary diagram (alternatively, _th, _no)" 1>&2
    echo "--frac_simple FRAC   Fraction simple projection (default:0.33)" 1>&2
    echo 1>&2
    exit 1
}

#####
#	PARSE COMMAND LINE
#####
SCALE="0.2"
JITTER="0"
APPEND=""
ALTER=""
FONTSZ=""
SS_SHFT="0,0"
TH_SHFT="0,0"
NO_SHFT="0,0"
FRAC_SIMPLE="0.33"
if [ $# -eq 0 ]; then usage; fi
while [ "$1" != "" ]; do
    case $1 in
        -f) shift;IFILE=$1;;
        -s) shift;SCALE=$1;;
        -c) shift;CPT=$1;;
        -j) shift;JITTER=$1;;
        -a) shift;APPEND=$1;;
        -x) shift;ALTER=$1;;
        -t) shift;FONTSZ=$1;;
        --shift_ss) shift;SS_SHFT="$1";;
        --shift_th) shift;TH_SHFT="$1";;
        --shift_no) shift;NO_SHFT="$1";;
        --frac_simple) shift;FRAC_SIMPLE="$1";;
        -h) usage;;
         *) echo No option $1;usage
    esac
    shift
done
if [ -z $SCALE ]
then
    echo "ternary.sh: no symbol scale defined, setting to 0.2" 1>&2
    SCALE="0.2"
fi
FRAC_SIMPLE=`echo $FRAC_SIMPLE | awk '{if($1<0){print 0}else if($1>1){print 1}else{print $1}}'`

#####
#	CHECK FOR INPUT FILE AND FORMAT
#####
if [ -z "$IFILE" ]
then
    echo "ternary.sh: no input file defined" 1>&2
    usage
elif [ ! -f "$IFILE" ]
then
    echo "ternary.sh: no input file $IFILE found" 1>&2
    usage
fi
MODE=`sed -n "1p" $IFILE | awk '{print NF}'`
case $MODE in
    3) OPT="-Sc${SCALE}i -Ggrey";;
    4) OPT="-Sci -Ggrey";;
    5) OPT="-Sci -C$CPT"
       if [ "$CPT" == "" ]
       then
           echo "ternary.sh: no color palette defined; plotting without colors" 1>&2
           MODE=4
           OPT="-Sci -Ggrey"
       elif [ ! -f "$CPT" ]
       then
           echo "ternary.sh: no color palette file named $CPT found; plotting without colors" 1>&2
           MODE=4
           OPT="-Sci -Ggrey"
       fi
       ;;
    *) echo "ternary.sh: No option for this many fields in input file" 1>&2; usage;;
esac

BIN_DIR="BIN_DIR_CHANGEME"

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
    echo 0 0 | gmt psxy $PROJ $LIMS -K -P --PS_MEDIA=8.5ix11i > $PSFILE
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
    ${BIN_DIR}/grid -x 0 $LFTOVR -dx 0.005 |\
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
            frac = '"$FRAC_SIMPLE"'
        }{
            tdipr = $1
            bdipr = $2
            pdipr = $3
            psi = atan2(sin(tdipr),sin(pdipr)) - 3.14159265/4
            denom = sin(aref)*sin(bdipr) + cos(aref)*cos(bdipr)*cos(psi)
            h = hgt*sqrt(2)*cos(bdipr)*sin(psi)/(3*denom)
            v = hgt*sqrt(2)*(cos(aref)*sin(bdipr)-sin(aref)*cos(bdipr)*cos(psi))/(3*denom)
            hs = hgt*(sin(tdipr)*sin(tdipr)-sin(pdipr)*sin(pdipr))/sqrt(3)
            vs = hgt*(sin(bdipr)*sin(bdipr)-1/3)
            h = frac*hs + (1-frac)*h
            v = frac*vs + (1-frac)*v
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
    ${BIN_DIR}/grid -x 0 $LFTOVR -dx 0.005 |\
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
            frac = '"$FRAC_SIMPLE"'
        }{
            tdipr = $1
            bdipr = $2
            pdipr = $3
            psi = atan2(sin(tdipr),sin(pdipr)) - 3.14159265/4
            denom = sin(aref)*sin(bdipr) + cos(aref)*cos(bdipr)*cos(psi)
            h = hgt*sqrt(2)*cos(bdipr)*sin(psi)/(3*denom)
            v = hgt*sqrt(2)*(cos(aref)*sin(bdipr)-sin(aref)*cos(bdipr)*cos(psi))/(3*denom)
            hs = hgt*(sin(tdipr)*sin(tdipr)-sin(pdipr)*sin(pdipr))/sqrt(3)
            vs = hgt*(sin(bdipr)*sin(bdipr)-1/3)
            h = frac*hs + (1-frac)*h
            v = frac*vs + (1-frac)*v
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
    ${BIN_DIR}/grid -x 0 $LFTOVR -dx 0.005 |\
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
            frac = '"$FRAC_SIMPLE"'
        }{
            tdipr = $1
            bdipr = $2
            pdipr = $3
            psi = atan2(sin(tdipr),sin(pdipr)) - 3.14159265/4
            denom = sin(aref)*sin(bdipr) + cos(aref)*cos(bdipr)*cos(psi)
            h = hgt*sqrt(2)*cos(bdipr)*sin(psi)/(3*denom)
            v = hgt*sqrt(2)*(cos(aref)*sin(bdipr)-sin(aref)*cos(bdipr)*cos(psi))/(3*denom)
            hs = hgt*(sin(tdipr)*sin(tdipr)-sin(pdipr)*sin(pdipr))/sqrt(3)
            vs = hgt*(sin(bdipr)*sin(bdipr)-1/3)
            h = frac*hs + (1-frac)*h
            v = frac*vs + (1-frac)*v
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
echo 0 $YMAX $FONT_SIZE,0 CB Slip |\
    gmt pstext $PROJ $LIMS -F+f+j -D${SS_DX}i/${SS_DY}i -N -K -O >> $PSFILE
SS_DY=`echo $SS_DY $FONT_SIZE | awk '{print $1+$2/72}'`
echo 0 $YMAX $FONT_SIZE,0 CB Strike |\
    gmt pstext $PROJ $LIMS -F+f+j -D${SS_DX}i/${SS_DY}i -N -K -O >> $PSFILE
NO_DX=`echo $NO_SHFT | awk -F, '{print $1}'`
NO_DY=`echo $NO_SHFT | awk -F, '{print $2-'"$OFFSET"'}'`
echo $XMIN $YMIN $FONT_SIZE,0 LT Normal |\
    gmt pstext $PROJ $LIMS -F+f+j -D${NO_DX}i/${NO_DY}i -N -K -O >> $PSFILE
TH_DX=`echo $TH_SHFT | awk -F, '{print $1}'`
TH_DY=`echo $TH_SHFT | awk -F, '{print $2-'"$OFFSET"'}'`
echo $XMAX $YMIN $FONT_SIZE,0 RT Thrust |\
    gmt pstext $PROJ $LIMS -F+f+j -D${TH_DX}i/${TH_DY}i -N -K -O >> $PSFILE

# Plot data
TIME=`date "+%s"`
awk 'BEGIN{
    aref = '"$AREF"'
    hgt  = '"$HEIGHT"'
    frac = '"$FRAC_SIMPLE"'
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
    hs = hgt*(sin(tdipr)*sin(tdipr)-sin(pdipr)*sin(pdipr))/sqrt(3)
    vs = hgt*(sin(bdipr)*sin(bdipr)-1/3)
    h = frac*hs + (1-frac)*h
    v = frac*vs + (1-frac)*v
    if (jitter>0) {
        h = h + (rand()-0.5)*jitter
        v = v + (rand()-0.5)*jitter
    }
    if ('"$MODE"'==3) {
        print h,v
    } else if ('"$MODE"'==4) {
        print h,v,$4*$4*$4*scl/100
    } else if ('"$MODE"'==5) {
        print h,v,$5,$4*$4*$4*scl/100
    }
}' $IFILE |\
    gmt psxy $PROJ $LIMS $OPT -W0.5p -N -K -O >> $PSFILE

if [ -z $APPEND ]
then
    echo 0 0 | gmt psxy $PROJ $LIMS -O >> $PSFILE
    ps2pdf $PSFILE
    rm $PSFILE
fi
