#!/bin/bash

###############################################################################
# Script for automatically computing and plotting surface displacements
# generated by an earthquake in an elastic half-space.
###############################################################################

###############################################################################
#	PARSE COMMAND LINE
###############################################################################
function usage() {
    echo "Usage: surf_disp.sh SRC_TYPE SRC_FILE [...options...]" 1>&2
    echo 1>&2
    echo "Required arguments" 1>&2
    echo "SRC_TYPE            MT, FLT, FFM, or FSP" 1>&2
    echo "SRC_FILE            Name of input fault file" 1>&2
    echo "                      MT:  evlo evla evdp str dip rak mag" 1>&2
    echo "                      FLT: evlo evla evdp str dip rak slip wid len" 1>&2
    echo "                      FFM: finite fault model in USGS .param format" 1>&2
    echo "                      FSP: finite fault model in SRCMOD FSP format" 1>&2
    echo 1>&2
    echo "Optional arguments (many of these are defined automatically)" 1>&2
    echo "-Rw/e/s/n           Map limits" 1>&2
    echo "-Tvmin/vmax/dv      Vertical displacement color palette" 1>&2
    echo "-vec_scale SCALE    Horizontal vector scale" 1>&2
    echo "-vec_legend LENGTH  Legend vector length (m)" 1>&2
    echo "-autothr THR        Auto displacement threshold (default: 1 mm)" 1>&2
    echo "-fade DISP_THR      Fade displacements below DISP_THR meters (default: 50 mm)" 1>&2
    echo "-nvert NN           Number of vertical contour grid points (default: 100/dimension)" 1>&2
    echo "-nvec NN_SAMP       Number of vectors (default: 20/dimension)" 1>&2
    echo "-seg                Plot segmented finite faults" 1>&2
    echo "-novector           Do not plot horizontal vectors" 1>&2
    echo "-gps GPS_FILE       Add observed displacements" 1>&2
    echo "-sta STA_FILE       Add only station locations" 1>&2
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
elif [ $# -lt 2 ]
then
    echo "surf_disp.sh: SRC_TYPE and SRC_FILE arguments required" 1>&2
    usage
fi
SRC_TYPE="$1"
SRC_FILE="$2"
shift
shift

# Check that the source type is an available option
if [ $SRC_TYPE != "FFM" -a $SRC_TYPE != "MT" -a $SRC_TYPE != "FSP" -a $SRC_TYPE != "FLT" ]
then
    echo "surf_disp.sh: source type must be FFM, FSP, MT, or FLT" 1>&2
    usage
fi

# Check that input file exists
if [ ! -f $SRC_FILE ]
then
    echo "surf_disp.sh: no source file $SRC_FILE found" 1>&2
    usage
fi


# Parse optional arguments
LIMS=""
VERT_CPT_RANGE=""
VEC_SCALE=""
DISP_LBL=""
NN="100" # Background vertical displacement grid is (NN x NN) points
NN_SAMP="20" # Horizontal vectors grid is (NN_SAMP x NN_SAMP) points
SEG="0"
PLOT_VECT="Y"
GPS_FILE=""
STA_FILE=""
EMPREL="WC"
DISP_THR="" # Horizontal displacements below DISP_THR meters will be faded
AUTO_THR="0.001" # Map limits defined based on this displacement threshold
HDEF_BIN_DIR=""
OFILE="surf_disp"
CLEAN="Y"
while [ "$1" != "" ]
do
    case $1 in
        -R*) LIMS="$1";;
        -T*) VERT_CPT_RANGE="$1";;
        -vec_scale) shift;VEC_SCALE="$1" ;;
        -vec_legend) shift;DISP_LBL="$1" ;;
        -autothr) shift;AUTO_THR="$1";;
        -fade) shift;DISP_THR="$1";;
        -nvert) shift;NN="$1";NN=$(echo $NN | awk '{printf("%d"),$1}');;
        -nvec) shift;NN_SAMP="$1";NN_SAMP=$(echo $NN_SAMP | awk '{printf("%d"),$1}');;
        -seg) SEG="1" ;;
        -novec*) PLOT_VECT="N" ;;
        -gps) shift;GPS_FILE="$1" ;;
        -sta) shift;STA_FILE="$1" ;;
        -emprel) shift;EMPREL="$1";;
        -hdefbin) shift;HDEF_BIN_DIR="$1";;
        -o) shift;OFILE="$1" ;;
        -noclean) CLEAN="N";;
        *) echo "surf_disp.sh: no option \"$1\"" 1>&2; usage;;
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
        echo "surf_disp.sh: Hdef executables not found in user-specified HDEF_BIN_DIR=$HDEF_BIN_DIR" 1>&2
        usage
    fi
else
    BIN_DIR=$(which o92util | xargs dirname)
    echo "surf_disp.sh: using Hdef executables found in BIN_DIR=$BIN_DIR" 1>&2
    if [ "$BIN_DIR" == "" ]
    then
        echo "surf_disp.sh: unable to find Hdef executables in PATH; exiting" 1>&2
        usage
    fi
fi


# GMT executables are required for this script!
GMT_DIR=$(which gmt | xargs dirname)
if [ "$GMT_DIR" == "" ]
then
    echo "surf_disp.sh: unable to find GMT executables; exiting" 1>&2
    exit 1
fi


###############################################################################
#	CLEAN UP FUNCTION
###############################################################################
function cleanup () {
    rm -f polar_mwh.cpt
    rm -f vert.cpt
    rm -f vert_scale_max.awk
    rm -f vert_scale_lbl.awk
    rm -f vect_label.awk
    rm -f vect_scale.awk
    rm -f vert.grd
    rm -f slip.grd
    rm -f o92util_auto_lims.dat
    rm -f gmt.*
    rm -f *.tmp
}
if [ "$CLEAN" == "Y" ]
then
    trap "cleanup" 0 1 2 3 8 9
fi


###############################################################################
#	DEFINE APPEARANCE OF DISPLACEMENT VECTORS AND CONTOURS
###############################################################################
# Vertical displacement color palette
cat > polar_mwh.cpt << EOF
# Simulates the POLAR colormap in Matlab
# Modified to make small values white
-1    blue    -0.1    white
-0.1    white    0.1    white
0.1    white    1    red
EOF


# Define the value at which the color bar for vertical displacements
# will saturate, based on maximum vertical displacements.
# IF (MAXIMUM VERTICAL DISPLACEMENT >= THRESHOLD) {USE THIS SATURATION VALUE}
cat > vert_scale_max.awk << EOF
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
  else if (\$1>=0.002) {print 0.001}
  else {print 0.0005}
}
EOF

# Define the annotation increment for the vertical displacement scale bar,
# based on the saturation value above.
# IF (MAXIMUM VERTICAL DISPLACEMENT >= THRESHOLD) {USE THIS ANNOTATION INCREMENT}
cat > vert_scale_lbl.awk << EOF
{
  if (\$1>=5) {print 2}
  else if (\$1>=2) {print 1}
  else if (\$1>=1) {print 0.5}
  else if (\$1>=0.5) {print 0.2}
  else if (\$1>=0.2) {print 0.1}
  else if (\$1>=0.1) {print 0.05}
  else if (\$1>=0.05) {print 0.02}
  else if (\$1>=0.02) {print 0.01}
  else if (\$1>=0.01) {print 0.005}
  else if (\$1>=0.005) {print 0.002}
  else if (\$1>=0.002) {print 0.001}
  else {print 0.0005}
}
EOF

# Use the maximum horizontal displacement to define the length of the
# vector in the legend.
# IF (MAXIMUM HORIZONTAL DISPLACEMENT >= THRESHOLD) {USE THIS LENGTH IN METERS AS LEGEND VECTOR}
cat > vect_label.awk << EOF
{
  if (\$1>10) {print 5}
  else if (\$1>5) {print 2}
  else if (\$1>1) {print 1}
  else if (\$1>0.5) {print 0.5}
  else if (\$1>0.1) {print 0.1}
  else if (\$1>0.05) {print 0.05}
  else if (\$1>0.02) {print 0.02}
  else if (\$1>0.01) {print 0.01}
  else if (\$1>0.005) {print 0.005}
  else {print 0.002}
}
EOF

# Use the maximum horizontal displacement to define the vector scaling.
# Larger earthquakes should have a smaller scale factor for all of the
# vectors to fit on the map.
# IF (MAXIMUM HORIZONTAL DISPLACEMENT >= THRESHOLD) {USE THIS VECTOR SCALING}
cat > vect_scale.awk << EOF
{
  if (\$1>10) {print 0.3}
  else if (\$1>5) {print 0.8}
  else if (\$1>1) {print 1.6}
  else if (\$1>0.5) {print 3.2}
  else if (\$1>0.2) {print 5}
  else if (\$1>0.1) {print 10}
  else if (\$1>0.05) {print 15}
  else if (\$1>0.02) {print 30}
  else if (\$1>0.01) {print 60}
  else if (\$1>0.005) {print 120}
  else if (\$1>0.002) {print 250}
  else if (\$1>0.001) {print 500}
  else {print 1000}
}
EOF


###############################################################################
###############################################################################
# Everything below this point *should* be automated. This script requires the
# tools O92UTIL, GRID, and FF2GMT from Matt's codes, and creates the figure
# using GMT 5/6 commands. All of the work is performed in the same directory
# that the script is run from.
###############################################################################
###############################################################################

#####
#	INPUT FILES FOR DISPLACEMENT CALCULATION
#####
# Copy source file to temporary file
cp $SRC_FILE ./source.tmp || { echo "surf_disp.sh: error copying EQ source file" 1>&2; exit 1; }

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
    D="10"  # Large initial increment, to get map limits without taking much time
    if [ $SRC_TYPE == "FFM" ]
    then
        ${BIN_DIR}/o92util -ffm source.tmp -auto h $Z $D -auto:thr $AUTO_THR -haf haf.tmp -disp disp.tmp || \
            { echo "surf_disp.sh: error running o92util with FFM source" 1>&2; exit 1; }
    elif [ $SRC_TYPE == "FSP" ]
    then
        ${BIN_DIR}/o92util -fsp source.tmp -auto h $Z $D -auto:thr $AUTO_THR -haf haf.tmp -disp disp.tmp  || \
            { echo "surf_disp.sh: error running o92util with FSP source" 1>&2; exit 1; }
    elif [ $SRC_TYPE == "MT" ]
    then
        ${BIN_DIR}/o92util -mag source.tmp -auto h $Z $D -auto:thr $AUTO_THR -haf haf.tmp -disp disp.tmp  || \
            { echo "surf_disp.sh: error running o92util with MT source" 1>&2; exit 1; }
    elif [ $SRC_TYPE == "FLT" ]
    then
        ${BIN_DIR}/o92util -flt source.tmp -auto h $Z $D -auto:thr $AUTO_THR -haf haf.tmp -disp disp.tmp  || \
            { echo "surf_disp.sh: error running o92util with FLT source" 1>&2; exit 1; }
    else
        echo "surf_disp.sh: no source type named \"$SRC_TYPE\"" 1>&2
        usage
    fi

    gmt gmtinfo -C disp.tmp > lims.tmp || \
        { echo "surf_disp.sh: error determining disp.tmp limits" 1>&2; exit 1; }
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
if [ $NN_SAMP -le 1 ]; then echo "NN_SAMP must be greater than 1" 1>&2; usage; fi
${BIN_DIR}/grid -x $W $E -nx $NN -y $S $N -ny $NN -z $Z -o sta.tmp || \
    { echo "surf_disp.sh: error in program grid" 1>&2; exit 1; }
if [ -z $GPS_FILE ]
then
    # Create (NN x NN) point horizontal grid for vectors
    ${BIN_DIR}/grid -x $W $E -nx $NN_SAMP -y $S $N -ny $NN_SAMP -z $Z -o sta_samp.tmp || \
        { echo "surf_disp.sh: error in program grid" 1>&2; exit 1; }
else
    # Take points from GPS file for vectors
    awk '{print $1,$2,0}' $GPS_FILE > sta_samp.tmp || \
        { echo "surf_disp.sh: error extracting points from GPS file" 1>&2; exit 1; }
fi


#####
#	COMPUTE SURFACE DISPLACEMENTS
#####
if [ $SRC_TYPE == "FFM" ]
then
    ${BIN_DIR}/o92util -ffm source.tmp -sta sta.tmp -haf haf.tmp -disp disp.tmp -prog || \
        { echo "surf_disp.sh: error running o92util with FFM source" 1>&2; exit 1; }
    ${BIN_DIR}/o92util -ffm source.tmp -sta sta_samp.tmp -haf haf.tmp -disp disp_samp.tmp -prog || \
        { echo "surf_disp.sh: error running o92util with FFM source" 1>&2; exit 1; }
elif [ $SRC_TYPE == "FSP" ]
then
    ${BIN_DIR}/o92util -fsp source.tmp -sta sta.tmp -haf haf.tmp -disp disp.tmp -prog || \
        { echo "surf_disp.sh: error running o92util with FSP source" 1>&2; exit 1; }
    ${BIN_DIR}/o92util -fsp source.tmp -sta sta_samp.tmp -haf haf.tmp -disp disp_samp.tmp -prog || \
        { echo "surf_disp.sh: error running o92util with FSP source" 1>&2; exit 1; }
elif [ $SRC_TYPE == "MT" ]
then
    ${BIN_DIR}/o92util -mag source.tmp -sta sta.tmp -haf haf.tmp -disp disp.tmp -prog -empirical ${EMPREL} || \
        { echo "surf_disp.sh: error running o92util with MT source" 1>&2; exit 1; }
    ${BIN_DIR}/o92util -mag source.tmp -sta sta_samp.tmp -haf haf.tmp -disp disp_samp.tmp -prog -empirical $EMPREL || \
        { echo "surf_disp.sh: error running o92util with MT source" 1>&2; exit 1; }
elif [ $SRC_TYPE == "FLT" ]
then
    ${BIN_DIR}/o92util -flt source.tmp -sta sta.tmp -haf haf.tmp -disp disp.tmp -prog || \
        { echo "surf_disp.sh: error running o92util with FLT source" 1>&2; exit 1; }
    ${BIN_DIR}/o92util -flt source.tmp -sta sta_samp.tmp -haf haf.tmp -disp disp_samp.tmp -prog || \
        { echo "surf_disp.sh: error running o92util with FLT source" 1>&2; exit 1; }
else
    echo "surf_disp.sh: no source type named $SRC_TYPE" 1>&2
    usage
fi

# Extract maximum vertical displacements and determine scale parameters for gridding
MINMAX=`awk '{print $6}' disp.tmp |\
        awk 'BEGIN{mn=1e10;mx=-1e10}{
            if($1<mn){mn=$1}
            if($1>mx){mx=$1}
        }END{print mn,mx}'`
echo "Minimum vertical displacement: $(echo $MINMAX | awk '{printf("%.3f m\n"),$1}')"
echo "Maximum vertical displacement: $(echo $MINMAX | awk '{printf("%.3f m\n"),$2}')"
V1=`echo $MINMAX | awk '{if($1<0){print -$1}else{print $1}}'`
V2=`echo $MINMAX | awk '{if($2<0){print -$2}else{print $2}}'`
T=`echo $V1 $V2 | awk '{if($1>$2){print $1}else{print $2}}' | awk -f vert_scale_max.awk`
COLOR_STEP=$(echo $T | awk '{print $1/100}')
DT=`echo $T | awk -f vert_scale_lbl.awk`


#####
#	PLOT RESULTS
#####
gmt set PS_MEDIA 8.5ix11i

PORTRAIT=`echo $X $Y | awk '{if($1<$2){print "-P"}}'`
PROJ="-JM5i -P"
LIMS="-R$W/$E/$S/$N"

# Colored grid of vertical displacements plotted under horizontal displacement vectors
if [ -z $VERT_CPT_RANGE ]
then
    gmt makecpt -T-${T}/${T}/${COLOR_STEP} -C./polar_mwh.cpt -D > vert.cpt || \
        { echo "surf_disp.sh: makecpt error" 1>&2; exit 1; }
else
    gmt makecpt $VERT_CPT_RANGE -C./polar_mwh.cpt -D > vert.cpt || \
        { echo "surf_disp.sh: makecpt error" 1>&2; exit 1; }
    DT=`echo $VERT_CPT_RANGE | awk -F/ '{print $2/2}'`
fi
awk '{print $1,$2,$6}' disp.tmp | gmt xyz2grd -Gvert.grd $LIMS -I${NN}+/${NN}+ || \
    { echo "surf_disp.sh: xyz2grd error" 1>&2; exit 1; }
gmt grdimage vert.grd $PROJ $LIMS -Cvert.cpt -Y1.5i -K > $PSFILE || \
    { echo "surf_disp.sh: grdimage error" 1>&2; exit 1; }
gmt psscale -Dx0i/-0.9i+w5.0i/0.2i+h+ml -Cvert.cpt -Ba$DT -Bg$DT -B+l"Vertical Displacement (m)" -K -O >> $PSFILE || \
    { echo "surf_disp.sh: psscale error" 1>&2; exit 1; }

# Map stuff
ANNOT=`echo $W $E | awk '{if($2-$1<=10){print 1}else{print 2}}'`
gmt psbasemap $PROJ $LIMS -Bxa${ANNOT} -Bya1 -BWeSn -K -O --MAP_FRAME_TYPE=plain >> $PSFILE || \
    { echo "surf_disp.sh: psbasemap error" 1>&2; exit 1; }

GS_VERSION=$(gs --version)
if [ "$GS_VERSION" == "9.24" -o "$GS_VERSION" == "9.25" -o "$GS_VERSION" == "9.51" -o "$GS_VERSION" == "9.52" ]
then
    echo "Just a heads up - some versions of Ghostscript do not support GMT transparency" 1>&2
    echo "Ghostscript 9.50 works fine for me" 1>&2
    echo "See: https://github.com/GenericMappingTools/gmt/issues/2903" 1>&2
fi
gmt pscoast $PROJ $LIMS -W1p,105/105/105 -G205/205/205 -N1/0.5p -Df -K -O -t85 >> $PSFILE || \
    { echo "surf_disp.sh: pscoast error" 1>&2; exit 1; }

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
            { echo "surf_disp.sh: ff2gmt error" 1>&2; exit 1; }
    else
        ${BIN_DIR}/ff2gmt $OPT -slip slip.tmp -clipseg clip.tmp -epi epi.tmp || \
            { echo "surf_disp.sh: ff2gmt error" 1>&2; exit 1; }
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
        { echo "surf_disp.sh: error making contour definition file" 1>&2; exit 1; }
    awk '{print $1,$2,$3}' slip.tmp |\
        gmt surface -Gslip.grd -I0.10/0.10 -Tb1 -Ti0.25 $LIMS || \
        { echo "surf_disp.sh: GMT surface error" 1>&2; exit 1; }
    gmt psclip clip.tmp $PROJ $LIMS -K -O >> $PSFILE || \
        { echo "surf_disp.sh: psclip error" 1>&2; exit 1; }
    gmt grdcontour slip.grd $PROJ $LIMS -W1p,205/205/205 -Cjunk -K -O -t40 >> $PSFILE || \
        { echo "surf_disp.sh: grdcontour error" 1>&2; exit 1; }
    gmt psclip -C -K -O >> $PSFILE || \
        { echo "surf_disp.sh: psclip error" 1>&2; exit 1; }
    gmt psxy clip.tmp $PROJ $LIMS -W1p,205/205/205 -K -O -t40 >> $PSFILE || \
        { echo "surf_disp.sh: psxy error" 1>&2; exit 1; }
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
        { echo "surf_disp.sh: psxy error plotting epicenter" 1>&2; exit 1; }
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


# Plot vectors
if [ $PLOT_VECT == "Y" ]
then
    if [ -z $GPS_FILE ]
    then
        # The maximum displacement can sometimes be much larger than other displacements if it sits
        # near one of the fault edges. Scaling this vector to a reasonable size will make all of the
        # other vectors too small to see.

        # If max displacement is much larger than the second largest displacement, don't use it
        MAXLN=`awk '{print sqrt($4*$4+$5*$5)}' disp_samp.tmp |\
               awk 'BEGIN{m1=0;m2=0}
                    {if($1>m1){m2=m1;m1=$1;ln=NR}}
                    END{if(m1>2*m2){print ln}else{print 0}}'`
    else
        MAXLN=0
    fi

    # Scale vectors based on maximum horizontal displacement
    MAX=`awk '{if(NR!='"$MAXLN"'){print sqrt($4*$4+$5*$5)}}' disp_samp.tmp |\
         awk 'BEGIN{mx=0}{if($1>mx){mx=$1}}END{print mx}' | awk '{printf("%.3f"),$1}'`
    echo "Maximum horizontal displacement: $MAX m"
    if [ -z $DISP_LBL ]
    then
        DISP_LBL=`echo $MAX | awk -f vect_label.awk`
    fi
    if [ -z $VEC_SCALE ]
    then
        VEC_SCALE=`echo $MAX | awk -f vect_scale.awk`
    fi

    # Set the faded displacement threshold based on the maximum horizontal displacement
    if [ "$DISP_THR" == "" ]
    then
        DISP_THR=$(echo $MAX |\
                   awk '{
                       thr = $1/10
                       if (thr>=0.01) {thr=0.01}
                       else if (thr>=0.005) {thr=0.005}
                       else if (thr>=0.002) {thr=0.002}
                       else {thr=0.001}
                       print thr
                   }')
    fi


    # Reduce scale for vectors shorter than VNORM
    VNORM=0.5

    # Vector head angle
    VANGLE=40

    # Vector head length
    VHEAD=9p


    # Ploting differs slightly depending on whether data are gridded or at defined GPS stations
    if [ -z $GPS_FILE ]
    then

        # No GPS file: plot vectors on grid

        # Plot displacements smaller than DISP_THR faded
        awk '{
            if (sqrt($4*$4+$5*$5)<'"$DISP_THR"') {
              print $1,$2,atan2($4,$5)/0.01745,'"$VEC_SCALE"'*sqrt($4*$4+$5*$5)
            }
        }' disp_samp.tmp |\
            gmt psxy $PROJ $LIMS -SV${VHEAD}+e+a${VANGLE}+n${VNORM} -W2p,175/175/175 -K -O >> $PSFILE || \
            { echo "surf_disp.sh: psxy error plotting modeled vectors" 1>&2; exit 1; }

        # Plot larger displacements in black
        awk '{
            if (sqrt($4*$4+$5*$5)>='"$DISP_THR"'&&NR!='"$MAXLN"') {
              print $1,$2,atan2($4,$5)/0.01745,'"$VEC_SCALE"'*sqrt($4*$4+$5*$5)
            }
        }' disp_samp.tmp |\
            gmt psxy $PROJ $LIMS -SV${VHEAD}+e+a${VANGLE}+n${VNORM} -W2p,black -K -O >> $PSFILE || \
            { echo "surf_disp.sh: psxy error plotting modeled vectors" 1>&2; exit 1; }

    else

        # There is a GPS file: plot vectors at GPS stations

        # Color vertical motions at stations same as background synthetic verticals
        awk '{print $1,$2,$5}' $GPS_FILE |\
            gmt psxy $PROJ $LIMS -Sc0.06i -W0.5p -Cvert.cpt -K -O >> $PSFILE || \
            { echo "surf_disp.sh: psxy error plotting observed vertical GPS" 1>&2; exit 1; }
        awk '{if(NF==6)print $1,$2,"4,0 LM",$6}' $GPS_FILE |\
            gmt pstext $PROJ $LIMS -F+f+j -D0.04i/0 -K -O >> $PSFILE || \
            { echo "surf_disp.sh: pstext error plotting GPS station name" 1>&2; exit 1; }

        # Plot horizontal GPS displacements in black
        awk '{print $1,$2,atan2($3,$4)/0.01745,'"$VEC_SCALE"'*sqrt($3*$3+$4*$4)}' $GPS_FILE |\
            gmt psxy $PROJ $LIMS -SV${VHEAD}+e+a${VANGLE}+n${VNORM} -W2p,black -K -O >> $PSFILE || \
            { echo "surf_disp.sh: psxy error plotting GPS vectors" 1>&2; exit 1; }

        # Plot synthetic displacements in another color
        awk '{print $1,$2,atan2($4,$5)/0.01745,'"$VEC_SCALE"'*sqrt($4*$4+$5*$5)}' disp_samp.tmp |\
            gmt psxy $PROJ $LIMS -SV${VHEAD}+e+a${VANGLE}+n${VNORM} -W2p,orange -K -O >> $PSFILE || \
            { echo "surf_disp.sh: psxy error plotting modeled vectors" 1>&2; exit 1; }
    fi
fi

# Plot only station locations
if [ "$STA_FILE" != "" ]
then
    gmt psxy $STA_FILE $PROJ $LIMS -Sc0.06i -Ggreen -W0.5p -K -O >> $PSFILE || \
	    { echo "surf_disp.sh: psxy error plotting GPS stations" 1>&2; exit 1; }
fi


# Legend (all coordinates are in cm from the bottom left)
if [ $PLOT_VECT == "Y" ]
then
    XLEG_ARROW=$(echo $VEC_SCALE $DISP_LBL | awk '{print $1*$2+0.6}')
    XLEG_TEXT=3.0
    XLEG_TEXT=$(echo $DISP_LBL |\
                awk '{
                    if ($1<0.01) {
                        print 3.1
                    } else if ($1<0.1) {
                        print 2.9
                    } else if ($1<1) {
                        print 2.7
                    } else if ($1<10) {
                        print 2.4
                    } else {
                        print 2.6
                    }
                }')
    XLEG=$(echo $XLEG_ARROW $XLEG_TEXT | awk '{if($1>$2){print $1}else{print $2}}')
    XMID=$(echo $XLEG 0.2 | awk '{print ($1+$2)/2}')

    # Legend box
    echo 0.2 0.2 > legend.tmp
    echo 0.2 1.5 >> legend.tmp
    echo $XLEG 1.5 >> legend.tmp
    echo $XLEG 0.2 >> legend.tmp
    echo 0.2 0.2 >> legend.tmp
    gmt psxy legend.tmp -JX10c -R0/10/0/10 -W1p -Gwhite -K -O >> $PSFILE || \
        { echo "surf_disp.sh: psxy error plotting legend outline" 1>&2; exit 1; }

    # Legend vector
    echo $VEC_SCALE $DISP_LBL |\
        awk '{print '$XMID',0.5,0,$1*$2}' |\
        gmt psxy -JX -R -Sv${VHEAD}+e+a${VANGLE}+jc -W2p,black -N -K -O >> $PSFILE || \
        { echo "surf_disp.sh: psxy error plotting legend vector" 1>&2; exit 1; }

    # Legend vector label
    echo $VEC_SCALE $DISP_LBL |\
        awk '{
            if ($2!=1) {
                print '$XMID',1.0,12","0,"CM",$2,"meters"
            } else {
                print '$XMID',1.0,12","0,"CM",$2,"meter"
            }
        }' |\
        gmt pstext -JX -R -F+f+j -N -K -O >> $PSFILE || \
        { echo "surf_disp.sh: pstext error plotting legend vector scale" 1>&2; exit 1; }

    # Color explanation text
    if [ -z $GPS_FILE ]
    then
        echo $VEC_SCALE $DISP_LBL |\
            awk '{print '$XLEG'+0.1,"0.2 10,2 LB Displacements less than '"$DISP_THR"' m are in light grey"}' |\
            gmt pstext -JX -R -F+f+j -Gwhite -N -K -O >> $PSFILE || \
            { echo "surf_disp.sh: pstext error plotting legend text" 1>&2; exit 1; }
    else
        echo $VEC_SCALE $DISP_LBL |\
            awk '{print '$XLEG'+0.1,"0.2 10,2 LB Observed=black; Synthetic=color"}' |\
            gmt pstext -JX -R -F+f+j -Gwhite -N -K -O >> $PSFILE || \
            { echo "surf_disp.sh: pstext error plotting legend text" 1>&2; exit 1; }
    fi
else
    VEC_SCALE=0
    DISP_LBL=0
fi


if [ $SRC_TYPE == "FFM" -o $SRC_TYPE == "FSP" ]
then
    echo $VEC_SCALE $DISP_LBL $CONT |\
        awk '{
          if($3==1) {print '$XLEG'+0.1,0.6,"10,2 LB FFM Slip Contours: "$3" meter"}
          else      {print '$XLEG'+0.1,0.6,"10,2 LB FFM Slip Contours: "$3" meters"}
        }' |\
        gmt pstext -JX10c -R0/10/0/10 -F+f+j -N -K -O >> $PSFILE || \
        { echo "surf_disp.sh: pstext error plotting legend contour scale" 1>&2; exit 1; }
fi

echo 0 0 | gmt psxy $PROJ $LIMS -O >> $PSFILE || \
    { echo "surf_disp.sh: psxy error finalizing PostScript file" 1>&2; exit 1; }

#####
#	CLEAN UP
#####
gmt psconvert $PSFILE -Tf
gmt psconvert $PSFILE -Tg -A
