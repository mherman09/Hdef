#!/bin/bash

###############################################################################
# Script for automatically computing and plotting surface displacements
# generated by an earthquake in an elastic half-space.
###############################################################################

###############################################################################
#	PARSE COMMAND LINE
###############################################################################
function usage() {
    echo "Usage: insar.sh SRC_TYPE SRC_FILE AZ INC WVL [...options...]" 1>&2
    echo 1>&2
    echo "Required arguments" 1>&2
    echo "SRC_TYPE            MT, FLT, FFM, or FSP" 1>&2
    echo "SRC_FILE            Name of input fault file" 1>&2
    echo "                      MT:  evlo evla evdp str dip rak mag" 1>&2
    echo "                      FLT: evlo evla evdp str dip rak slip wid len" 1>&2
    echo "                      FFM: finite fault model in USGS .param format" 1>&2
    echo "                      FSP: finite fault model in SRCMOD FSP format" 1>&2
    echo "AZ                  Look azimuth (satellite to ground)" 1>&2
    echo "INC                 Look inclination (0=horizontal, satellite to ground)" 1>&2
    echo "WVL                 Radar wavelength (e.g., Sentinel 1: 0.056 m)" 1>&2
    echo 1>&2
    echo "Optional arguments (many of these are defined automatically)" 1>&2
    echo "-Rw/e/s/n           Map limits" 1>&2
    echo "-autothr  THR       Auto displacement threshold (default: 1 mm)" 1>&2
    echo "-ngrid NN           Number of vertical contour grid points (default: 100/dimension)" 1>&2
    echo "-seg                Plot segmented finite faults" 1>&2
    echo "-emprel EMPREL      Empirical relation for rect source" 1>&2
    echo "-hdefbin DIR        Directory with Hdef executables" 1>&2
    echo "-o FILENAME         Basename for output file" 1>&2
    echo "-noclean            Keep all temporary files (useful for debugging)" 1>&2
    echo 1>&2
    exit 1
}

# Source type and source file are required
if [ $# -eq 0 ]
then
    usage
elif [ $# -lt 5 ]
then
    echo "insar.sh: SRC_TYPE, SRC_FILE, and satellite look arguments required" 1>&2
    usage
fi
SRC_TYPE="$1"
SRC_FILE="$2"
AZ="$3"
INC="$4"
WVL="$5"
shift; shift; shift; shift; shift


# Check that the source type is an available option
if [ $SRC_TYPE != "FFM" -a $SRC_TYPE != "MT" -a $SRC_TYPE != "FSP" -a $SRC_TYPE != "FLT" ]
then
    echo "insar.sh: source type must be FFM, FSP, MT, or FLT" 1>&2
    usage
fi

# Check that input file exists
if [ ! -f $SRC_FILE ]
then
    echo "insar.sh: no source file $SRC_FILE found" 1>&2
    usage
fi

# Check that azimuth, inclination, and wavelength are numbers
CHECK_LOOK=$(echo $AZ $INC $WVL |\
    awk 'BEGIN{e=0}{
        if ($1~/[A-Za-z]/) {print "insar.sh: error in azimuth input:",$1 >> "/dev/stderr";e=1}
        if ($2~/[A-Za-z]/) {print "insar.sh: error in inclination input:",$2 >> "/dev/stderr";e=1}
        if ($3~/[A-Za-z]/) {print "insar.sh: error in wavelength input:",$3 >> "/dev/stderr";e=1}
    }END{print e}')
if [ "$CHECK_LOOK" == "1" ]
then
    usage
fi


# Parse optional arguments
LIMS=""
NN="100" # Background vertical displacement grid is (NN x NN) points
SEG="0"
EMPREL="WC"
AUTO_THR="0.001" # Map limits defined based on this displacement threshold
HDEF_BIN_DIR=""
OFILE="insar"
CLEAN="Y"
while [ "$1" != "" ]
do
    case $1 in
        -R*) LIMS="$1";;
        -autothr) shift;AUTO_THR="$1";;
        -ngrid) shift;NN="$1";NN=$(echo $NN | awk '{printf("%d"),$1}');;
        -seg) SEG="1" ;;
        -emprel) shift;EMPREL="$1";;
        -hdefbin) shift;HDEF_BIN_DIR="$1";;
        -o) shift;OFILE="$1" ;;
        -noclean) CLEAN="N";;
        *) echo "insar.sh: no option \"$1\"" 1>&2; usage;;
    esac
    shift
done

PSFILE="$OFILE.ps"


###############################################################################
#	CHECK FOR REQUIRED EXECUTABLES
###############################################################################

# Check for executables in user-specified directory, if defined
if [ "$HDEF_BIN_DIR" != "" ]
then
    BIN_DIR=$(which $HDEF_BIN_DIR/o92util | xargs dirname)
    if [ "$BIN_DIR" == "" ]
    then
        echo "insar.sh: Hdef executables not found in user-specified HDEF_BIN_DIR=$HDEF_BIN_DIR" 1>&2
        usage
    fi
else
    BIN_DIR=$(which o92util | xargs dirname)
    echo "insar.sh: using Hdef executables found in BIN_DIR=$BIN_DIR" 1>&2
    if [ "$BIN_DIR" == "" ]
    then
        echo "insar.sh: unable to find Hdef executables in PATH; exiting" 1>&2
        usage
    fi
fi


# GMT executables are required for this script!
GMT_DIR=$(which gmt | xargs dirname)
if [ "$GMT_DIR" == "" ]
then
    echo "insar.sh: unable to find GMT executables; exiting" 1>&2
    exit 1
fi


###############################################################################
#	CLEAN UP FUNCTION
###############################################################################
function cleanup () {
    rm -f o92util_auto_lims.dat
    rm -f gmt.*
    rm -f *.tmp
}
if [ "$CLEAN" == "Y" ]
then
    trap "cleanup" 0 1 2 3 8 9
fi


###############################################################################
#	DEFINE APPEARANCE OF DISPLACEMENT CONTOURS
###############################################################################
# LOS phase color palette
gmt makecpt -T-3.15/3.15/0.02 -Crainbow -D > phase.cpt
# if [ -f /usr/local/cpt/cet-c2.cpt ]
# then
#     gmt makecpt -T-3.15/3.15/0.02 -C/usr/local/cpt/cet-c2.cpt -D > phase.cpt
# fi
if [ -f ~/Downloads/ScientificColourMaps6/romaO/romaO.cpt ]
then
    gmt makecpt -T-3.15/3.15/0.02 -C/Users/mwh5316/Downloads/ScientificColourMaps6/romaO/romaO.cpt -D > phase.cpt
fi


# Define the value at which the color bar for LOS displacements
# will saturate, based on maximum LOS displacements.
# IF (MAXIMUM LOS DISPLACEMENT >= THRESHOLD) {USE THIS SATURATION VALUE}
cat > los_scale_max.awk << EOF
{
  if (\$1>=20) {print 10}
  else if (\$1>=10) {print 5}
  else if (\$1>=5) {print 2}
  else if (\$1>=2) {print 1}
  else if (\$1>=1) {print 0.5}
  else if (\$1>=0.5) {print 0.2}
  else if (\$1>=0.2) {print 0.1}
  else if (\$1>=0.1) {print 0.05}
  else if (\$1>=0.05) {print 0.02}
  else if (\$1>=0.02) {print 0.01}
  else if (\$1>=0.01) {print 0.005}
  else if (\$1>=0.005) {print 0.002}
  else {print 0.001}
}
EOF

# Define the annotation increment for the LOS displacement scale bar,
# based on the saturation value above.
# IF (MAXIMUM LOS DISPLACEMENT >= THRESHOLD) {USE THIS ANNOTATION INCREMENT}
cat > los_scale_lbl.awk << EOF
{
  if (\$1>=10) {print 5}
  else if (\$1>=5) {print 2}
  else if (\$1>=2) {print 1}
  else if (\$1>=1) {print 0.5}
  else if (\$1>=0.5) {print 0.2}
  else if (\$1>=0.2) {print 0.1}
  else if (\$1>=0.1) {print 0.05}
  else if (\$1>=0.05) {print 0.02}
  else if (\$1>=0.02) {print 0.01}
  else if (\$1>=0.01) {print 0.005}
  else if (\$1>=0.005) {print 0.002}
  else {print 0.001}
}
EOF


###############################################################################
###############################################################################
# Everything below this point *should* be automated. This script requires the
# tools O92UTIL, GRID, VEC2LOS, WRAPLOS, and FF2GMT from Matt's codes, and
# creates the figure using GMT 5/6 commands. All of the work is performed in
# the same directory that the script is run from.
###############################################################################
###############################################################################

#####
#	INPUT FILES FOR DISPLACEMENT CALCULATION
#####
# Copy source file to temporary file
cp $SRC_FILE ./source.tmp || { echo "insar.sh: error copying EQ source file" 1>&2; exit 1; }

# Elastic half-space properties
LAMDA="4e10"   # Lame parameter
MU="4e10"      # Shear modulus
echo "lame $LAMDA shearmod $MU" > haf.tmp


#####
#	SET UP COMPUTATION GRID
#####
Z="0.0" # Depth is zero on the surface
if [ -z $LIMS ]
then
    # Use "-auto" option in O92UTIL to get rough map limits
    NINIT="10"  # Large initial increment, to get map limits without taking much time
    if [ $SRC_TYPE == "FFM" ]
    then
        ${BIN_DIR}/o92util -ffm source.tmp -auto h $Z $NINIT -auto:thr $AUTO_THR -haf haf.tmp -disp disp.tmp || \
            { echo "insar.sh: error running o92util with FFM source" 1>&2; exit 1; }
    elif [ $SRC_TYPE == "FSP" ]
    then
        ${BIN_DIR}/o92util -fsp source.tmp -auto h $Z $NINIT -auto:thr $AUTO_THR -haf haf.tmp -disp disp.tmp  || \
            { echo "insar.sh: error running o92util with FSP source" 1>&2; exit 1; }
    elif [ $SRC_TYPE == "MT" ]
    then
        ${BIN_DIR}/o92util -mag source.tmp -auto h $Z $NINIT -auto:thr $AUTO_THR -haf haf.tmp -disp disp.tmp  || \
            { echo "insar.sh: error running o92util with MT source" 1>&2; exit 1; }
    elif [ $SRC_TYPE == "FLT" ]
    then
        ${BIN_DIR}/o92util -flt source.tmp -auto h $Z $NINIT -auto:thr $AUTO_THR -haf haf.tmp -disp disp.tmp  || \
            { echo "insar.sh: error running o92util with FLT source" 1>&2; exit 1; }
    else
        echo "insar.sh: no source type named \"$SRC_TYPE\"" 1>&2
        usage
    fi

    gmt gmtinfo -C disp.tmp > lims.tmp || \
        { echo "insar.sh: error determining disp.tmp limits" 1>&2; exit 1; }
    W=`awk '{printf("%.3f"), $1}' lims.tmp`
    E=`awk '{printf("%.3f"), $2}' lims.tmp`
    S=`awk '{printf("%.3f"), $3}' lims.tmp`
    N=`awk '{printf("%.3f"), $4}' lims.tmp`
    echo "Starting map limits: $W $E $S $N"

    # Determine if map has decent aspect ratio and correct as necessary
    # Mercator projection x and y lengths
    X=`echo $W $E | awk '{print $2-$1}'`
    Y=`echo $S $N |\
       awk '{
         v2 = log(sin(3.14159/4+$2/2*0.01745)/cos(3.14159/4+$2/2*0.01745))
         v1 = log(sin(3.14159/4+$1/2*0.01745)/cos(3.14159/4+$1/2*0.01745))
         print v2-v1
       }' |\
       awk '{print $1/0.017}'`

    # Check map aspect ratio (no skinnier than 1.4:1)
    FIX=`echo $X $Y |\
         awk '{
           if ($1>1.4*$2) {print "fixx"}
           else if ($2>1.4*$1) {print "fixy"}
           else {print 1}
         }'`

    # Reduce map limits in long dimension
    if [ $FIX == "fixx" ]
    then
        NEW=`echo $W $E $Y | awk '{print 0.5*($1+$2)-$3*0.70,0.5*($1+$2)+$3*0.70}'`
        W=`echo $NEW | awk '{print $1}'`
        E=`echo $NEW | awk '{print $2}'`
    elif [ $FIX == "fixy" ]
    then
        NEW=`echo $S $N $X $Y |\
             awk '{print 0.5*($1+$2)-0.7*$3/$4*($2-$1),0.5*($1+$2)+0.7*$3/$4*($2-$1)}'`
        S=`echo $NEW | awk '{print $1}'`
        N=`echo $NEW | awk '{print $2}'`
    fi
    # Round map limits to nearest 0.1 degree
    W=`echo "$W $E" | awk '{printf("%.1f"),$1}'`
    E=`echo "$W $E" | awk '{printf("%.1f"),$2}'`
    S=`echo "$S $N" | awk '{printf("%.1f"),$1}'`
    N=`echo "$S $N" | awk '{printf("%.1f"),$2}'`
    echo "Final map limits:    $W $E $S $N"

else
    # Use map limits specified on command line
    W=`echo $LIMS | sed -e "s/\// /g" -e "s/-R//" | awk '{print $1}'`
    E=`echo $LIMS | sed -e "s/\// /g" -e "s/-R//" | awk '{print $2}'`
    S=`echo $LIMS | sed -e "s/\// /g" -e "s/-R//" | awk '{print $3}'`
    N=`echo $LIMS | sed -e "s/\// /g" -e "s/-R//" | awk '{print $4}'`
    echo "Using map limits from command line: $W $E $S $N"
fi


# Locations of displacement computations
if [ $NN -le 1 ]; then echo "NN must be greater than 1" 1>&2; usage; fi
${BIN_DIR}/grid -x $W $E -nx $NN -y $S $N -ny $NN -z $Z -o sta.tmp || \
    { echo "insar.sh: error in program grid" 1>&2; exit 1; }


#####
#	COMPUTE SURFACE DISPLACEMENTS
#####
if [ $SRC_TYPE == "FFM" ]
then
    ${BIN_DIR}/o92util -ffm source.tmp -sta sta.tmp -haf haf.tmp -disp disp.tmp -prog || \
        { echo "insar.sh: error running o92util with FFM source" 1>&2; exit 1; }
elif [ $SRC_TYPE == "FSP" ]
then
    ${BIN_DIR}/o92util -fsp source.tmp -sta sta.tmp -haf haf.tmp -disp disp.tmp -prog || \
        { echo "insar.sh: error running o92util with FSP source" 1>&2; exit 1; }
elif [ $SRC_TYPE == "MT" ]
then
    ${BIN_DIR}/o92util -mag source.tmp -sta sta.tmp -haf haf.tmp -disp disp.tmp -prog -empirical ${EMPREL} || \
        { echo "insar.sh: error running o92util with MT source" 1>&2; exit 1; }
elif [ $SRC_TYPE == "FLT" ]
then
    ${BIN_DIR}/o92util -flt source.tmp -sta sta.tmp -haf haf.tmp -disp disp.tmp -prog || \
        { echo "insar.sh: error running o92util with FLT source" 1>&2; exit 1; }
else
    echo "insar.sh: no source type named $SRC_TYPE" 1>&2
    usage
fi


# Calculate LOS displacements and wrap them onto radar phase
$BIN_DIR/vec2los -f disp.tmp -o los.tmp -a ${AZ} -i ${INC}
$BIN_DIR/wraplos -f los.tmp -o phase.tmp -w ${WVL}

# Extract maximum LOS displacements and determine scale parameters for gridding
MINMAX=`awk '{print $4}' los.tmp |\
        awk 'BEGIN{mn=1e10;mx=-1e10}{
            if($1<mn){mn=$1}
            if($1>mx){mx=$1}
        }END{print mn,mx}'`
echo "Minimum LOS displacement: $(echo $MINMAX | awk '{printf("%.3f m\n"),$1}')"
echo "Maximum LOS displacement: $(echo $MINMAX | awk '{printf("%.3f m\n"),$2}')"
ULOS1=`echo $MINMAX | awk '{if($1<0){print -$1}else{print $1}}'`
ULOS2=`echo $MINMAX | awk '{if($2<0){print -$2}else{print $2}}'`
T=`echo $ULOS1 $ULOS2 | awk '{if($1>$2){print $1}else{print $2}}' | awk -f los_scale_max.awk`
DT=`echo $T | awk -f los_scale_lbl.awk`


#####
#	PLOT RESULTS
#####
gmt set PS_MEDIA 8.5ix11i

PORTRAIT=`echo $X $Y | awk '{if($1<$2){print "-P"}}'`
PROJ="-JM5i -P"
LIMS="-R$W/$E/$S/$N"


# Plot phase grid
awk '{print $1,$2,$4}' phase.tmp | gmt xyz2grd -Gphase.grd $LIMS -I${NN}+/${NN}+ || \
    { echo "insar.sh: xyz2grd error" 1>&2; exit 1; }
gmt grdimage phase.grd $PROJ $LIMS -Cphase.cpt -Y1.5i -K -nn+a > $PSFILE || \
    { echo "insar.sh: grdimage error" 1>&2; exit 1; }

# Scale bar for cyclic phase color palette
XCENTER=0.7
YCENTER=-0.8
DANGLE=5
PHASE_SCALE_DIAMETER=1.0
INNER_RING_FACTOR=0.60
$BIN_DIR/grid -x 0 355 -dx $DANGLE |\
    awk '{
        x = '$XCENTER'
        y = '$YCENTER'
        phase = ($1+'$DANGLE'/2)*0.0174533; if (phase>3.14) {phase = phase - 6.28}
        print x,y,phase,$1-1-'$DANGLE'/10,$1+'$DANGLE'+'$DANGLE'/10
    }' |\
    gmt psxy -JX1i -R0/1/0/1 -Sw${PHASE_SCALE_DIAMETER}i -Cphase.cpt -N -K -O >> $PSFILE || \
    { echo "insar.sh: error drawing phase scale bar" 1>&2; exit 1; }

# Scale bar outlines
echo $XCENTER $YCENTER $PHASE_SCALE_DIAMETER |\
    gmt psxy -JX1i -R0/1/0/1 -Sci -W1p -N -K -O >> $PSFILE || \
    { echo "insar.sh: error drawing outer outline of phase scale bar" 1>&2; exit 1; }
echo $XCENTER $YCENTER $PHASE_SCALE_DIAMETER |\
    awk '{print $1,$2,$3*'$INNER_RING_FACTOR'}' |\
    gmt psxy -JX1i -R0/1/0/1 -Sci -W1p -Gwhite -N -K -O >> $PSFILE || \
    { echo "insar.sh: error drawing inner outline of phase scale bar" 1>&2; exit 1; }

# Labels for the scale bar
echo $XCENTER $YCENTER $PHASE_SCALE_DIAMETER |\
    awk '{print $1+$3/2+0.2, $2+0.1, "16,0 LM LOS Displacement (Wrapped)"}' |\
    gmt pstext -JX1i -R0/1/0/1 -F+f+j -N -K -O >> $PSFILE || \
    { echo "insar.sh: error labeling phase scale bar" 1>&2; exit 1; }
echo $XCENTER $YCENTER $PHASE_SCALE_DIAMETER |\
    awk '{print $1+$3/2+0.2, $2-0.15, "10,2 LM (Positive: away from satellite)"}' |\
    gmt pstext -JX1i -R0/1/0/1 -F+f+j -N -K -O >> $PSFILE || \
    { echo "insar.sh: error labeling phase scale positive direction" 1>&2; exit 1; }

# Each cycle of colors represents one full wavelength, indicating the ground has moved 1/2 that
# amount (to account for the two-way travel of the signal)
PHASE_SCALE_YSHFT=2
$BIN_DIR/grid -x 0 355 -dx 1 |\
    awk 'BEGIN{
        x0 = '$XCENTER'
        y0 = '$YCENTER'+'$PHASE_SCALE_YSHFT'
        rad = 0.5*'$PHASE_SCALE_DIAMETER'*'$INNER_RING_FACTOR'-0.05
        d2r = 0.0174533
    }{
        x = x0 + rad*cos($1*d2r)
        y = y0 + rad*sin($1*d2r)
        print x,y
    }' |\
    gmt psxy -JX2i -R0/2/0/2 -W1p+ve5p+gblack+a40 -Ya-${PHASE_SCALE_YSHFT}i -K -O >> $PSFILE || \
    { echo "insar.sh: error drawing positive cycle arrow" 1>&2; exit 1; }
WVL2=$(echo $WVL | awk '{printf("%g"),$1/2*100}')
echo $XCENTER $YCENTER 12,2 CB $WVL2 |\
    gmt pstext -JX1i -R0/1/0/1 -F+f+j -D0/0.01i -N -K -O >> $PSFILE || \
    { echo "insar.sh: error labeling scale wavelength" 1>&2; exit 1; }
echo $XCENTER $YCENTER 12,2 CT cm |\
    gmt pstext -JX1i -R0/1/0/1 -F+f+j -D0/-0.02i -N -K -O >> $PSFILE || \
    { echo "insar.sh: error labeling scale wavelength" 1>&2; exit 1; }

# Positive direction
$BIN_DIR/grid -x -10 45 -dx 1 |\
    awk 'BEGIN{
        x0 = '$XCENTER'
        y0 = '$YCENTER'+'$PHASE_SCALE_YSHFT'
        rad = 0.5*'$PHASE_SCALE_DIAMETER'+0.05
        d2r = 0.0174533
    }{
        x = x0 + rad*cos($1*d2r)
        y = y0 + rad*sin($1*d2r)
        print x,y
    }' |\
    gmt psxy -JX2i -R0/2/0/2 -W1p+ve5p+gblack+a40 -Ya-${PHASE_SCALE_YSHFT}i -K -O >> $PSFILE || \
    { echo "insar.sh: error drawing positive direction arrow outside cyclic scale" 1>&2; exit 1; }
echo 45 |\
    awk 'BEGIN{
        x0 = '$XCENTER'
        y0 = '$YCENTER'
        rad = 0.5*'$PHASE_SCALE_DIAMETER'+0.05
        d2r = 0.0174533
    }{
        x = x0 + rad*cos($1*d2r)
        y = y0 + rad*sin($1*d2r)
        print x,y,"10,0 LM +"
    }' |\
    gmt pstext -JX1i -R0/1/0/1 -F+f+j -D0.090i/-0.045i -N -K -O >> $PSFILE || \
    { echo "insar.sh: error labeling scale wavelength" 1>&2; exit 1; }

# Map frame
ANNOT=`echo $W $E | awk '{if($2-$1<=10){print 1}else{print 2}}'`
gmt psbasemap $PROJ $LIMS -Bxa${ANNOT} -Bya1 -BWeSn -K -O --MAP_FRAME_TYPE=plain >> $PSFILE || \
    { echo "insar.sh: psbasemap error" 1>&2; exit 1; }


# Coastline
GS_VERSION=$(gs --version)
if [ "$GS_VERSION" == "9.24" -o "$GS_VERSION" == "9.25" -o "$GS_VERSION" == "9.51" -o "$GS_VERSION" == "9.52" ]
then
    echo "Just a heads up - some versions of Ghostscript do not support GMT transparency" 1>&2
    echo "Ghostscript 9.50 works fine for me" 1>&2
    echo "See: https://github.com/GenericMappingTools/gmt/issues/2903" 1>&2
fi
gmt pscoast $PROJ $LIMS -W1p,105@85 -G205@85 -S255@5 -N1/0.5p,black@85 -Dh -K -O >> $PSFILE || \
    { echo "insar.sh: pscoast error" 1>&2; exit 1; }


# Plot FFM slip contours
if [ $SRC_TYPE == "FFM" -o $SRC_TYPE == "FSP" ]
then
    case $SRC_TYPE in
        FFM) OPT="-ffm source.tmp";;
        FSP) OPT="-fsp source.tmp";;
    esac
    if [ $SEG -eq 0 ]
    then
        ${BIN_DIR}/ff2gmt $OPT -slip slip.tmp -clip clip.tmp -epi epi.tmp || \
            { echo "insar.sh: ff2gmt error" 1>&2; exit 1; }
    else
        ${BIN_DIR}/ff2gmt $OPT -slip slip.tmp -clipseg clip.tmp -epi epi.tmp || \
            { echo "insar.sh: ff2gmt error" 1>&2; exit 1; }
    fi
    MAXSLIP=`awk '{print $3}' slip.tmp |\
             awk 'BEGIN{mx=0}{if($1>mx){mx=$1}}END{print mx}' |\
             awk '{print $1}'`
    CONT=`echo $MAXSLIP |\
          awk '{
            if ($1>=50) {print 10}
            else if ($1>=20) {print 5}
            else if ($1>=10) {print 2}
            else if ($1>=2) {print 1}
            else {print 0.5}
          }'`
    echo $CONT $MAXSLIP | awk '{for (i=$1;i<=$2;i=i+$1){print i,"C"}}' > junk || \
        { echo "insar.sh: error making contour definition file" 1>&2; exit 1; }
    awk '{print $1,$2,$3}' slip.tmp |\
        gmt surface -Gslip.grd -I0.10/0.10 -Tb1 -Ti0.25 $LIMS || \
        { echo "insar.sh: GMT surface error" 1>&2; exit 1; }
    gmt psclip clip.tmp $PROJ $LIMS -K -O >> $PSFILE || \
        { echo "insar.sh: psclip error" 1>&2; exit 1; }
    gmt grdcontour slip.grd $PROJ $LIMS -W1p,205/205/205 -Cjunk -K -O -t40 >> $PSFILE || \
        { echo "insar.sh: grdcontour error" 1>&2; exit 1; }
    gmt psclip -C -K -O >> $PSFILE || \
        { echo "insar.sh: psclip error" 1>&2; exit 1; }
    gmt psxy clip.tmp $PROJ $LIMS -W1p,205/205/205 -K -O -t40 >> $PSFILE || \
        { echo "insar.sh: psxy error" 1>&2; exit 1; }
    rm junk
# else
    # awk '{print $1,$2,$4,$5,$6}' rect.out |\
    #     gmt psxy $PROJ $LIMS -SJ -W1p,205/205/205 -K -O -t40 >> $PSFILE
fi


# Plot epicenter
if [ $SRC_TYPE == "FFM" -o $SRC_TYPE == "FSP" ]
then
    LONX=`awk '{print $1}' epi.tmp`
    LATX=`awk '{print $2}' epi.tmp`
    echo $LONX $LATX |\
        gmt psxy $PROJ $LIMS -Sa0.15i -W1p,55/55/55 -K -O -t50 >> $PSFILE || \
        { echo "insar.sh: psxy error plotting epicenter" 1>&2; exit 1; }
fi


# Plot focal mechanisms
if [ $SRC_TYPE == "MT" ]
then
    MAX_MAG=$(awk '{print $7}' $SRC_FILE | gmt gmtinfo -C | awk '{print $2}')
    MAX_MAG_DIAM=0.2
    MEC_SCALE=$(echo $MAX_MAG | awk '{print 5/($1*$1*$1)}')
    awk '{print $1,$2,$3,$4,$5,$6,$7*$7*$7*'"$MEC_SCALE"'}' $SRC_FILE |\
        gmt psmeca $PROJ $LIMS -Sa${MAX_MAG_DIAM}i -L0.5p,55 -G105 -t25 -K -O >> $PSFILE
fi


# Legend (all coordinates are in cm from the bottom left)
# Legend box
echo 0.2 0.2 > legend.tmp
echo 0.2 2.0 >> legend.tmp
echo 2.0 2.0 >> legend.tmp
echo 2.0 0.2 >> legend.tmp
echo 0.2 0.2 >> legend.tmp
gmt psxy legend.tmp -JX10c -R0/10/0/10 -W1p -Gwhite -K -O >> $PSFILE || \
    { echo "insar.sh: psxy error plotting legend outline" 1>&2; exit 1; }

# Legend vectors
AZMINUS90=$(echo $AZ | awk '{print $1-90}')
echo 1.1 1.1 ${AZMINUS90} 1.5 | gmt psxy -JX10c -R0/10/0/10 -SV6p+e+jc+a40 -W1p -Gblack -K -O >> $PSFILE
echo 1.1 1.1 ${AZ} 0.7 | gmt psxy -JX10c -R0/10/0/10 -SV6p+e+jb+a40 -W1p -Gblack -K -O >> $PSFILE
echo 1.1 1.1 10,0 CM ${INC} | gmt pstext -JX10c -R0/10/0/10 -F+f+j -Gwhite -W1p -N -K -O >> $PSFILE

if [ $SRC_TYPE == "FFM" -o $SRC_TYPE == "FSP" ]
then
    echo 0 0 $CONT |\
        awk '{
          if($3==1) {print '$XLEG'+0.1,0.2,"10,2 LB FFM Slip Contours: "$3" meter"}
          else      {print '$XLEG'+0.1,0.2,"10,2 LB FFM Slip Contours: "$3" meters"}
        }' |\
        gmt pstext -JX10c -R0/10/0/10 -F+f+j -N -K -O >> $PSFILE || \
        { echo "insar.sh: pstext error plotting legend contour scale" 1>&2; exit 1; }
fi










# Plot displacement grid
echo 0 0 | gmt psxy $PROJ $LIMS -Y${Y}i -K -O >> $PSFILE


# Plot displacement grid
$BIN_DIR/colortool -lightness 50,100 -hue 130,90 -chroma 50,0 -gmt -T-1/0/0.05 > los0.cpt
$BIN_DIR/colortool -lightness 100,40 -hue -30,-60 -chroma 0,50 -gmt -T0/1/0.05 >> los0.cpt
gmt makecpt -T-${T}/${T} -C./los0.cpt -D > los.cpt || \
	{ echo "insar.sh: makecpt error" 1>&2; exit 1; }
awk '{print $1,$2,$4}' los.tmp | gmt xyz2grd -Glos.grd $LIMS -I${NN}+/${NN}+ || \
    { echo "insar.sh: xyz2grd error" 1>&2; exit 1; }
gmt grdimage los.grd $PROJ $LIMS -Clos.cpt -Y2.4i -nn+a -K -O >> $PSFILE || \
    { echo "insar.sh: grdimage error" 1>&2; exit 1; }

# Scale bar
gmt psscale -Dx0i/-1.1i+w5.0i/0.2i+h+mal -Clos.cpt -Ba${DT}g${DT}+l"LOS Displacement" -K -O >> $PSFILE || \
    { echo "insar.sh: psscale error on unwrapped LOS displacement" 1>&2; exit 1; }
echo "2.5 -1.2 10,2 CT (Positive: away from satellite)" |\
    gmt pstext -JX1i -R0/1/0/1 -F+f+j -N -K -O >> $PSFILE

# Map frame
ANNOT=`echo $W $E | awk '{if($2-$1<=10){print 1}else{print 2}}'`
gmt psbasemap $PROJ $LIMS -Bxa${ANNOT} -Bya1 -BWeSn -K -O --MAP_FRAME_TYPE=plain >> $PSFILE || \
    { echo "insar.sh: psbasemap error" 1>&2; exit 1; }


# Coastline
gmt pscoast $PROJ $LIMS -W1p,105@85 -G205@85 -S255@5 -N1/0.5p,black@85 -Dh -K -O >> $PSFILE || \
    { echo "insar.sh: pscoast error" 1>&2; exit 1; }


# Plot FFM slip contours
if [ $SRC_TYPE == "FFM" -o $SRC_TYPE == "FSP" ]
then
    case $SRC_TYPE in
        FFM) OPT="-ffm source.tmp";;
        FSP) OPT="-fsp source.tmp";;
    esac
    if [ $SEG -eq 0 ]
    then
        ${BIN_DIR}/ff2gmt $OPT -slip slip.tmp -clip clip.tmp -epi epi.tmp || \
            { echo "insar.sh: ff2gmt error" 1>&2; exit 1; }
    else
        ${BIN_DIR}/ff2gmt $OPT -slip slip.tmp -clipseg clip.tmp -epi epi.tmp || \
            { echo "insar.sh: ff2gmt error" 1>&2; exit 1; }
    fi
    MAXSLIP=`awk '{print $3}' slip.tmp |\
             awk 'BEGIN{mx=0}{if($1>mx){mx=$1}}END{print mx}' |\
             awk '{print $1}'`
    CONT=`echo $MAXSLIP |\
          awk '{
            if ($1>=50) {print 10}
            else if ($1>=20) {print 5}
            else if ($1>=10) {print 2}
            else if ($1>=2) {print 1}
            else {print 0.5}
          }'`
    echo $CONT $MAXSLIP | awk '{for (i=$1;i<=$2;i=i+$1){print i,"C"}}' > junk || \
        { echo "insar.sh: error making contour definition file" 1>&2; exit 1; }
    awk '{print $1,$2,$3}' slip.tmp |\
        gmt surface -Gslip.grd -I0.10/0.10 -Tb1 -Ti0.25 $LIMS || \
        { echo "insar.sh: GMT surface error" 1>&2; exit 1; }
    gmt psclip clip.tmp $PROJ $LIMS -K -O >> $PSFILE || \
        { echo "insar.sh: psclip error" 1>&2; exit 1; }
    gmt grdcontour slip.grd $PROJ $LIMS -W1p,205/205/205 -Cjunk -K -O -t40 >> $PSFILE || \
        { echo "insar.sh: grdcontour error" 1>&2; exit 1; }
    gmt psclip -C -K -O >> $PSFILE || \
        { echo "insar.sh: psclip error" 1>&2; exit 1; }
    gmt psxy clip.tmp $PROJ $LIMS -W1p,205/205/205 -K -O -t40 >> $PSFILE || \
        { echo "insar.sh: psxy error" 1>&2; exit 1; }
    rm junk
# else
    # awk '{print $1,$2,$4,$5,$6}' rect.out |\
    #     gmt psxy $PROJ $LIMS -SJ -W1p,205/205/205 -K -O -t40 >> $PSFILE
fi


# Plot epicenter
if [ $SRC_TYPE == "FFM" -o $SRC_TYPE == "FSP" ]
then
    LONX=`awk '{print $1}' epi.tmp`
    LATX=`awk '{print $2}' epi.tmp`
    echo $LONX $LATX |\
        gmt psxy $PROJ $LIMS -Sa0.15i -W1p,55/55/55 -K -O -t50 >> $PSFILE || \
        { echo "insar.sh: psxy error plotting epicenter" 1>&2; exit 1; }
fi


# Plot focal mechanisms
if [ $SRC_TYPE == "MT" ]
then
    MAX_MAG=$(awk '{print $7}' $SRC_FILE | gmt gmtinfo -C | awk '{print $2}')
    MAX_MAG_DIAM=0.2
    MEC_SCALE=$(echo $MAX_MAG | awk '{print 5/($1*$1*$1)}')
    awk '{print $1,$2,$3,$4,$5,$6,$7*$7*$7*'"$MEC_SCALE"'}' $SRC_FILE |\
        gmt psmeca $PROJ $LIMS -Sa${MAX_MAG_DIAM}i -L0.5p,55 -G105 -t25 -K -O >> $PSFILE
fi


# Legend
gmt psxy legend.tmp -JX10c -R0/10/0/10 -W1p -Gwhite -K -O >> $PSFILE || \
    { echo "insar.sh: psxy error plotting legend outline" 1>&2; exit 1; }
echo 1.1 1.1 ${AZMINUS90} 1.5 | gmt psxy -JX10c -R0/10/0/10 -SV6p+e+jc+a40 -W1p -Gblack -K -O >> $PSFILE
echo 1.1 1.1 ${AZ} 0.7 | gmt psxy -JX10c -R0/10/0/10 -SV6p+e+jb+a40 -W1p -Gblack -K -O >> $PSFILE
echo 1.1 1.1 10,0 CM ${INC} | gmt pstext -JX10c -R0/10/0/10 -F+f+j -Gwhite -W1p -N -K -O >> $PSFILE



echo 0 0 | gmt psxy $PROJ $LIMS -O >> $PSFILE || \
    { echo "insar.sh: psxy error finalizing PostScript file" 1>&2; exit 1; }

#####
#	CLEAN UP
#####
gmt psconvert $PSFILE -Tf
gmt psconvert $PSFILE -Tg -A