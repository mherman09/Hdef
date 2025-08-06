#!/bin/bash


# Initialize script name and log file
SCRIPT=`basename $0`
TIMESTAMP=`date "+%s"`

LOG_FILE=$SCRIPT.log

echo "$SCRIPT [`date "+%H:%M:%S"`]: starting" | tee $LOG_FILE





####################################################################################################
# USAGE STATEMENT
####################################################################################################
# Usage statement
function usage() {
    echo "Usage: $0 -f SEIS_FILE -start DATE_START -end DATE_END -dday DDAY" 1>&2
    echo "       [...<options>...]" 1>&2
    echo "" 1>&2
    echo "Required Arguments" 1>&2
    echo "-f SEIS_FILE                  Seismicity file (origin_time lon lat dep mag [str dip rak])" 1>&2
    echo "-start DATE_START             Starting time" 1>&2
    echo "-end DATE_END                 Ending time" 1>&2
    echo "-dday DDAY                    Time increment between frames (days)" 1>&2
    echo "" 1>&2
    echo "Optional arguments" 1>&2
    echo "-J<opt>                       Map projection in GMT format (default: Mercator)" 1>&2
    echo "-Rw/e/s/n                     Map limits in GMT format (default: plot all seismicity in file)" 1>&2
    echo "-mag:color MAG_COLOR_CHANGE   Change color every time earthquake above MAG_COLOR_CHANGE occurs" 1>&2
    echo "-mag:min MAG_MIN              Minimum magnitude to include on map/graph (default: all)" 1>&2
    echo "-mag:bold MAG_BOLD            Make earthquakes larger than MAG_BOLD bolder" 1>&2
    echo "-mag:label MAG_LABEL          Label earthquakes larger than MAG_LABEL" 1>&2
    echo "-scale SCALE                  Adjust earthquake scaling (default: 0.001)" 1>&2
    echo "-dday:change DDAY,NDAYS       Change DDAY at NDAYS from start" 1>&2
    echo "-timecpt TIME_CPT             Earthquake timing (by day) color palette" 1>&2
    echo "-datecpt DATE_CPT             Earthquake timing (by date) color palette" 1>&2
    echo "-color-by-col4                Color earthquakes by value in column 4 of SEIS_FILE instead of date" 1>&2
    echo "-col4cpt CPT                  Color palette file for custom coloring" 1>&2
    echo "-fade FADE_TIME               Time to fade out symbols (days)" 1>&2
    echo "-trans:max MAX_TRANS          Maximum transparency of faded symbols (default:90)" 1>&2
    echo "-mech                         Plot focal mechanisms (default: no mechanisms)" 1>&2
    echo "-plates PB_FILE               Plot plate boundaries (default: no boundaries)" 1>&2
    echo "-topo TOPO_FILE               Plot shaded topography in background" 1>&2
    echo "-topo:color-mode MODE         BW, COLOR-LAND-ONLY, COLOR-OCEAN-ONLY, COLOR" 1>&2
    echo "-topo:clean ON|OFF            Remove topo and gradient file after script finishes (default:ON)" 1>&2
    echo "-other:psxy FILE:OPT          Other psxy file to plot (can repeat)" 1>&2
    echo "-other:pstext WORDS:OPT       Other pstext file to plot (can repeat)" 1>&2
    echo "-other:psimage FILE:OPT       Other psimage file to plot (can repeat)" 1>&2
    echo "-pscoast:options OPT          Additional flags for pscoast (separate by \":\")" 1>&2
    echo "-add-north-arrow JUST         Add north arrow in GMT justified corner (e.g., TR)" 1>&2
    echo "-mag-axis:min MAG_MIN         Minimum magnitude on y-axis of mag vs time frame" 1>&2
    echo "-mag-axis:max MAG_MAX         Maximum magnitude on y-axis of mag vs time frame" 1>&2
    echo "-no-mag-vs-time               Do not plot magnitude versus time" 1>&2
    echo "-o OUTPUT_NAME                Change name of output video file (default: seis_movie.mp4)" 1>&2
    echo "-parallel N                   Make frames on N processes (default: N=1)" 1>&2
    echo "-clean                        Remove frame files after running (default is to keep them)" 1>&2
    exit 1
}

if [ $# -eq 0 ]
then
    usage
fi




####################################################################################################
# SCRIPT VARIABLES
####################################################################################################


echo | tee -a $LOG_FILE
echo "$SCRIPT [`date "+%H:%M:%S"`]: parsing command line arguments" | tee -a $LOG_FILE



# Seismicity file: origin_time(YYYY-MM-DDTHH:MM:SS)  lon  lat  dep  mag  [str  dip  rak]
SEIS_FILE=
SEIS_FILE=seis.dat

# Animation parameters
DATE_START=                             # Starting date
DATE_END=                               # Ending date
DDAY=                                   # Number of days between frames
FADE_TIME=10000000                      # Time to fade earthquakes out (days)
MAX_TRANS=90                            # Maximum transparency level (out of 100)
DDAY_NDAY_CHANGE_LIST=

# Map parameters
MAP_PROJ="-JM5i"
MAP_LIMS=
MAG_MIN=-10
MAG_AXIS_MIN=2
MAG_AXIS_MAX=8
MAG_COLOR_CHANGE=
MAG_BOLD=15
MAG_LABEL=15
TIME_CPT=
DATE_CPT=
COL4_CPT=
COLOR_BY=DATE
PLOT_MECH=N
PLATE_BOUNDARY_FILE=
TOPO_FILE=
TOPO_CLEAN=ON
COLOR_MODE=BW
SEIS_SCALE=0.001
PLOT_MAG_VS_TIME=Y
ADD_NORTH_ARROW=
GRIDLINES=OFF

# Other mapping commands
PSXY_LIST=""
PSTEXT_LIST=""
PSIMAGE_LIST=""
PSCOAST_OPTIONS=""
PSCOAST_PEN="0.75p"

# Miscellaneous
NPROC=1
CLEAN="N"
OUTPUT_MP4_FILE=seis_movie.mp4


# Parse command line
while [ "$1" != "" ]
do
    case $1 in
        -f) shift; SEIS_FILE=$1;;
        -start) shift; DATE_START=$1;;
        -end) shift; DATE_END=$1;;
        -dday) shift; DDAY=$1;;
        -dday:change) shift; DDAY_NDAY_CHANGE_LIST="$DDAY_NDAY_CHANGE_LIST $1";;
        -R*) MAP_LIMS=$1;;
        -J*) MAP_PROJ=$1;;
        -mag:min) shift; MAG_MIN=$1;;
        -mag-axis:min) shift; MAG_AXIS_MIN=$1;;
        -mag-axis:max) shift; MAG_AXIS_MAX=$1;;
        -mag:color) shift; MAG_COLOR_CHANGE=$1;;
        -mag:bold) shift; MAG_BOLD=$1;;
        -mag:label) shift; MAG_LABEL=$1;;
        -scale) shift; SEIS_SCALE=$1;;
        -timecpt) shift; TIME_CPT=$1;;
        -datecpt) shift; DATE_CPT=$1;;
        -col4cpt) shift; COL4_CPT=$1;;
        -color-by-col4) COLOR_BY=COL4;;
        -no-mag-vs-time) PLOT_MAG_VS_TIME=N;;
        -fade) shift; FADE_TIME=$1;;
        -trans:max) shift; MAX_TRANS=$1;;
        -mech) PLOT_MECH=Y;;
        -plates) shift; PLATE_BOUNDARY_FILE=$1;;
        -topo) shift; TOPO_FILE=$1;;
        -topo:color-mode) shift; COLOR_MODE=$1;;
        -topo:clean) shift; TOPO_CLEAN=$1;;
        -other:psxy) shift;PSXY_LIST="$PSXY_LIST;$1";;
        -other:pstext) shift;PSTEXT_LIST="$PSTEXT_LIST;$1";;
        -other:psimage) shift;PSIMAGE_LIST="$PSIMAGE_LIST;$1";;
        -pscoast:options) shift; PSCOAST_OPTIONS=`echo $1 | sed -e "s/:/ /g"` ;;
        -pscoast:pen) shift; PSCOAST_PEN="$1" ;;
        -add-gridlines) shift; GRIDLINES=ON ;;
        -add-north-arrow) shift; ADD_NORTH_ARROW="$1" ;;
        -o) shift; OUTPUT_MP4_FILE="$1" ;;
        -parallel) shift; NPROC="$1" ;;
        -clean) CLEAN="Y";;
        *) ;;
    esac
    shift
done


# Check that file variables are defined and files exist
echo "$SCRIPT [`date "+%H:%M:%S"`]: checking that seismicity file variable is set" | tee -a $LOG_FILE
for VAR in SEIS_FILE
do
    if [ "${!VAR}" == "" ]
    then
        echo "$0: File $VAR is not defined" 1>&2
        usage
    fi
    if [ ! -f ${!VAR} ]
    then
        echo "$SCRIPT: Could not find $VAR named \"${!VAR}\"" 1>&2
        usage
    fi
done
echo "$SCRIPT [`date "+%H:%M:%S"`]: using seismicity file $SEIS_FILE" >> $LOG_FILE


# If MAP_LIMS is not defined, set it here from SEIS_FILE
if [ "$MAP_LIMS" == "" ]
then
    echo "$SCRIPT [`date "+%H:%M:%S"`]: getting map limits from the input file" | tee -a $LOG_FILE
    MAP_LIMS=$(gmt gmtinfo $SEIS_FILE -I0.01 -i1:2)
    echo "$SCRIPT [`date "+%H:%M:%S"`]: using map limits: MAP_LIMS=$MAP_LIMS" | tee -a $LOG_FILE
else
    echo "$SCRIPT [`date "+%H:%M:%S"`]: setting map limits from command line: MAP_LIMS=$MAP_LIMS" | tee -a $LOG_FILE
fi
echo $MAP_LIMS >> $LOG_FILE


# If plate boundary file is defined but does not exist, set to none
if [ "$PLATE_BOUNDARY_FILE" != "" ]
then
    if [ ! -f $PLATE_BOUNDARY_FILE ]
    then
        echo "$SCRIPT [`date "+%H:%M:%S"`]: could not find plate boundary file named \"$PLATE_BOUNDARY_FILE\"" | tee -a $LOG_FILE
        echo "$SCRIPT [`date "+%H:%M:%S"`]: not plotting plate boundaries" | tee -a $LOG_FILE
        PLATE_BOUNDARY_FILE=
    fi
fi


# Check that variables are defined
for VAR in DATE_START DATE_END DDAY MAP_LIMS
do
    if [ "${!VAR}" == "" ]
    then
        echo "$SCRIPT: Variable $VAR is not defined" 1>&2
        usage
    fi
done


# Update date format
if [[ $DATE_START == ????-??-?? ]]
then
    DATE_START=${DATE_START}T00:00:00
fi
if [[ $DATE_END == ????-??-?? ]]
then
    DATE_END=${DATE_END}T00:00:00
fi
echo "$SCRIPT [`date "+%H:%M:%S"`]: starting animation at $DATE_START" | tee -a $LOG_FILE
echo "$SCRIPT [`date "+%H:%M:%S"`]: ending animation at $DATE_END" | tee -a $LOG_FILE


# Determine magnitude for color changes
if [ "$TIME_CPT" == "" ]
then
    if [ "$MAG_COLOR_CHANGE" == "" ]
    then
        echo "$SCRIPT [`date "+%H:%M:%S"`]: determining magnitude to change colors" | tee -a $LOG_FILE
        MAG_COLOR_CHANGE=$(gmt gmtinfo $SEIS_FILE -i4 -C | awk '{print $2-1.0}')
        echo "$SCRIPT [`date "+%H:%M:%S"`]: changing colors every magnitude ${MAG_COLOR_CHANGE}+" | tee -a $LOG_FILE
    else
        echo "$SCRIPT [`date "+%H:%M:%S"`]: changing colors every magnitude ${MAG_COLOR_CHANGE}+ (from cmdln)" | tee -a $LOG_FILE
    fi
fi

# Look for user-defined cpt files
if [ "$DATE_CPT" != "" ]
then
    test -f $DATE_CPT || { echo "$SCRIPT: could not find date color palette named \"$DATE_CPT\""; usage; }
fi
if [ "$TIME_CPT" != "" ]
then
    test -f $TIME_CPT || { echo "$SCRIPT: could not find timing color palette named \"$TIME_CPT\""; usage; }
fi







####################################################################################################
#	SET ANIMATION PARAMETERS
####################################################################################################


echo | tee -a $LOG_FILE
echo "$SCRIPT [`date "+%H:%M:%S"`]: setting animation frame parameters" | tee -a $LOG_FILE


# Determine number of days in animation
NDAYS=$(echo $DATE_START $DATE_END | dateutil -nday -format "YYYY-MM-DDTHH:MM:SS")
echo "$SCRIPT [`date "+%H:%M:%S"`]: animation is a total of $NDAYS days" | tee -a $LOG_FILE


# Calculate number of frames
if [ "$DDAY_NDAY_CHANGE_LIST" == "" ]
then
    NFRAMES=$(echo $NDAYS $DDAY | awk '{printf("%d"),$1/$2+1}')
    IFRAME_CHANGE=1
    DDAY_CHANGE=$DDAY
    echo "$SCRIPT [`date "+%H:%M:%S"`]: plotting frames every $DDAY days" | tee -a $LOG_FILE
else
    NDAY0=0
    DDAY0=$DDAY
    NFRAMES=0
    IFRAME_CHANGE=1
    DDAY_CHANGE=$DDAY
    for DDAY_NDAY_CHANGE in $DDAY_NDAY_CHANGE_LIST
    do
        DDAY_NEW=`echo $DDAY_NDAY_CHANGE | awk -F, '{print $1}'`
        NDAY_CHANGE=`echo $DDAY_NDAY_CHANGE | awk -F, '{print $2}'`
        DFRAMES=`echo $NDAY0 $NDAY_CHANGE $DDAY0 | awk '{printf("%d"),($2-$1)/$3}'`
        NFRAMES=`echo $NFRAMES $DFRAMES | awk '{print $1+$2+1}'`
        IFRAME_CHANGE="$IFRAME_CHANGE $NFRAMES"
        DDAY_CHANGE="$DDAY_CHANGE $DDAY_NEW"
        NDAY0=$NDAY_CHANGE
        DDAY0=$DDAY_NEW
    done
    NCHANGE=`echo $IFRAME_CHANGE | awk '{print NF}'`
    DFRAMES=`echo $NDAY0 $NDAYS $DDAY0 | awk '{printf("%d"),($2-$1)/$3}'`
    NFRAMES=`echo $NFRAMES $DFRAMES | awk '{print $1+$2}'`
fi
echo "$SCRIPT [`date "+%H:%M:%S"`]: generating $NFRAMES frames" | tee -a $LOG_FILE



# GMT parameters
gmt set PS_MEDIA 100ix100i
gmt set FORMAT_GEO_OUT D
gmt set FORMAT_GEO_MAP D


# Map parameters
MAP_WID=`echo $MAP_PROJ | awk -F/ '{print $NF}'`
MAP_PROJ_FLAG=`echo $MAP_PROJ | grep -o  "\-J[A-Za-z]*"`
MAP_PROJ_PARAM=`echo $MAP_PROJ | sed -e "s/-J[A-Za-z]*//" | awk -F"/" '{for(i=1;i<=NF-1;i++){printf("%s/"),$i}}END{printf("'$MAP_WID'\n")}'`
echo "$SCRIPT [`date "+%H:%M:%S"`]: using MAP_PROJ=$MAP_PROJ" | tee -a $LOG_FILE
MAP_TIKS=$(echo $MAP_LIMS | sed -e "s/-R//" |\
           awk -F/ '{
                if (/\+r/) {
                   dlon = $3-$1
                } else {
                   dlon = $2-$1
                }
                if (dlon<0.1) {
                    print 0.02
                } else if (dlon<0.2) {
                    print 0.05
                } else if (dlon<0.5) {
                    print 0.1
                } else if (dlon<1.0) {
                    print 0.2
                } else if (dlon<2.0) {
                    print 0.5
                } else if (dlon<5.0) {
                    print 1.0
                } else if (dlon<10.0) {
                    print 2.0
                } else {
                    print 5.0
                }
            }')
echo "$SCRIPT [`date "+%H:%M:%S"`]: map projection: $MAP_PROJ" | tee -a $LOG_FILE
echo "$SCRIPT [`date "+%H:%M:%S"`]: map tick interval: $MAP_TIKS" | tee -a $LOG_FILE


# Magnitude versus time parameters
GRAPH_WID=`echo $MAP_WID | sed -e "s/i//"`
TIME_PROJ="-JX${GRAPH_WID}iT/1.5i"
DAY_TIME_PROJ=$(echo $TIME_PROJ | sed -e "s/T//")
TIME_LIMS="-R${DATE_START}/${DATE_END}/${MAG_AXIS_MIN}/${MAG_AXIS_MAX}"
DAY_END=$(echo $DATE_START $DATE_END | dateutil -nday -format "YYYY-MM-DDTHH:MM:SS" | awk '{print $1}') # LOOKS LIKE I CAN REPLACE THIS WITH $NDAYS
echo "$SCRIPT [`date "+%H:%M:%S"`]: ending day of animation: $DAY_END (REPLACE WITH \$NDAYS?)" | tee -a $LOG_FILE
YEAR_TIKS=$(echo $NDAYS |\
           awk '{
               if ($1<365*5) {
                   print 1
               } else if ($1<365*10) {
                   print 2
               } else if ($1<365*20) {
                   print 5
               } else if ($1<365*50) {
                   print 10
               } else {
                   print 20
               }
           }')
MONTH_TIKS=$(echo $NDAYS |\
           awk '{
               if ($1<365) {
                   print 1
               } else if ($1<365*2) {
                   print 2
               } else if ($1<365*3) {
                   print 3
               } else if ($1<365*4) {
                   print 6
               } else if ($1<365*6) {
                   print 12
               } else {
                   print 0
               }
           }')
DAY_TIKS=$(echo $NDAYS |\
           awk '{
               if ($1<2) {
                   print 0.25
               } else if ($1<5) {
                   print 0.5
               } else if ($1<10) {
                   print 1
               } else if ($1<20) {
                   print 2
               } else if ($1<50) {
                   print 5
               } else if ($1<100) {
                   print 10
               } else if ($1<200) {
                   print 20
               } else if ($1<500) {
                   print 50
               } else if ($1<1000) {
                   print 100
               } else if ($1<2000) {
                   print 200
               } else if ($1<5000) {
                   print 500
               } else if ($1<10000) {
                   print 1000
               } else if ($1<20000) {
                   print 2000
               } else {
                   print 5000
               }
           }')
MAG_TIKS=$(echo $TIME_LIMS | sed -e "s/-R//" |\
           awk -F/ '{
               dmag = $4-$3
               if (dmag<0.1) {
                   print 0.02
               } else if (dmag<0.3) {
                   print 0.05
               } else if (dmag<0.7) {
                   print 0.1
               } else if (dmag<1.0) {
                   print 0.2
               } else if (dmag<3.0) {
                   print 0.5
               } else if (dmag<7.0) {
                   print 1.0
               } else if (dmag<10.0) {
                   print 2.0
               } else {
                   print 5.0
               }
           }')
echo "$SCRIPT [`date "+%H:%M:%S"`]: mag vs time projection: $TIME_PROJ" | tee -a $LOG_FILE
echo "$SCRIPT [`date "+%H:%M:%S"`]: mag vs time limits: $TIME_LIMS" | tee -a $LOG_FILE
echo "$SCRIPT [`date "+%H:%M:%S"`]: mag vs time year tick interval: $YEAR_TIKS" | tee -a $LOG_FILE
echo "$SCRIPT [`date "+%H:%M:%S"`]: mag vs time month tick interval: $MONTH_TIKS" | tee -a $LOG_FILE
echo "$SCRIPT [`date "+%H:%M:%S"`]: mag vs time day tick interval: $DAY_TIKS" | tee -a $LOG_FILE
echo "$SCRIPT [`date "+%H:%M:%S"`]: mag vs time magnitude tick interval: $MAG_TIKS" | tee -a $LOG_FILE


# Select only data inside map frame
echo "$SCRIPT [`date "+%H:%M:%S"`]: selecting data within map frame" | tee -a $LOG_FILE
NCOLS=`head -1 $SEIS_FILE | awk '{print NF-1}'`
CMD="gmt gmtselect $SEIS_FILE $MAP_LIMS -i1:2,0,3:$NCOLS -o2,0:1,3:$NCOLS"
echo $CMD >> $LOG_FILE
$CMD | awk '{if($5>='$MAG_MIN'){print $0}}' > make_seis_movie_seis.tmp


# Calculate number of days since starting time for each earthquake
echo "$SCRIPT [`date "+%H:%M:%S"`]: calculating number of days since start time for each event" | tee -a $LOG_FILE
awk '{print "'$DATE_START'",$1}' make_seis_movie_seis.tmp |\
    dateutil -nday -format "YYYY-MM-DDTHH:MM:SS" > nday.tmp


# Color palette
EQ_CPT="${SCRIPT}.${TIMESTAMP}.cpt"
if [ "$TIME_CPT" == "" ]
then
    echo "$SCRIPT [`date "+%H:%M:%S"`]: generating timing color palette" | tee -a $LOG_FILE
    paste nday.tmp make_seis_movie_seis.tmp |\
        awk '{if($1>0 && $6>='$MAG_COLOR_CHANGE'){printf("%.5f\n"),$1}}' > make_seis_movie_nday_color.tmp
    NCOLORS=$(wc make_seis_movie_nday_color.tmp | awk '{print $1+1}')
    echo green > color_list.tmp
    echo blue >> color_list.tmp
    echo orange >> color_list.tmp
    echo 105/0/155 >> color_list.tmp
    echo yellow >> color_list.tmp
    echo magenta >> color_list.tmp
    colortool -hue -40,90 -lightness 15,95 -chroma 100,100 -ncolors $NCOLORS >> color_list.tmp
    gmt makecpt -Ccategorical -T0/$NCOLORS/1 | awk '{if(NR<='$NCOLORS'){print $2}}' > color_list.tmp
    awk '{if(NR<='$NCOLORS'){print $0}}' color_list.tmp > j; mv j color_list.tmp
    awk '{
        if (NR==1) {
            print 0,$1
        } else {
            print nday,$1
        }
        nday = $1
    }' make_seis_movie_nday_color.tmp |\
        paste - color_list.tmp |\
        awk '{if(NF>1){print $1,$3,$2,$3}}END{print "B black";print "F",$1}' > ${EQ_CPT}
    rm make_seis_movie_nday_color.tmp
    echo "$SCRIPT [`date "+%H:%M:%S"`]: built timing color palette ${EQ_CPT}" | tee -a $LOG_FILE
    cat ${EQ_CPT} >> $LOG_FILE
else
    echo "$SCRIPT [`date "+%H:%M:%S"`]: copying timing color palette $TIME_CPT to ${EQ_CPT}" | tee -a $LOG_FILE
    cp $TIME_CPT ${EQ_CPT}
fi

if [ "$DATE_CPT" != "" ]
then
    echo "$SCRIPT [`date "+%H:%M:%S"`]: using date-specified timing color palette $DATE_CPT" | tee -a $LOG_FILE
    awk '{if(NF==4){print "'"$DATE_START"'",$1;print "'"$DATE_START"'",$3}}' $DATE_CPT |\
        dateutil -nday -format "YYYY-MM-DDTHH:MM:SS" |\
        awk '{d1=$1;getline;d2=$1;print d1,d2}' > j
    paste j $DATE_CPT |\
        awk '{if(NF==6){print $1,$4,$2,$6}else{print $0}}' > ${EQ_CPT}
    rm j
    echo "$SCRIPT [`date "+%H:%M:%S"`]: created ${EQ_CPT}" | tee -a $LOG_FILE
fi

if [ "$COLOR_BY" == "COL4" ]
then
    if [ "$COL4_CPT" != "" ]
    then
        cp $COL4_CPT $EQ_CPT
        echo "$SCRIPT [`date "+%H:%M:%S"`]: copied $COL4_CPT to $EQ_CPT" | tee -a $LOG_FILE
    else
        CMD="gmt makecpt -T0/100/10 -Cplasma -D -I"
        $CMD > ${EQ_CPT}
        echo $CMD >> $LOG_FILE
        echo "$SCRIPT [`date "+%H:%M:%S"`]: made $EQ_CPT with $CMD" | tee -a $LOG_FILE
    fi
fi





####################################################################################################
#	CREATE ANIMATION BACKGROUND
####################################################################################################


echo | tee -a $LOG_FILE
echo "$SCRIPT [`date "+%H:%M:%S"`]: generating background features" | tee -a $LOG_FILE



PSFILE=frame_base.ps



# Initialize plot
gmt psxy -T -K -Y3.5i > $PSFILE



# Topography
if [ "$TOPO_FILE" != "" ]
then

    # Topography gradient for shaded hillslopes
    if [ ! -f topo_cut_grad.grd ]
    then
        echo "$SCRIPT [`date "+%H:%M:%S"`]: could not find topo gradient file...generating" | tee -a $LOG_FILE
        gmt grdcut $TOPO_FILE $MAP_PROJ $MAP_LIMS -Gtopo_cut.grd
        gmt grdgradient topo_cut.grd -A0/270 -Ne0.4 -Gtopo_cut_grad.grd
    else
        echo "$SCRIPT [`date "+%H:%M:%S"`]: using existing topo gradient file" | tee -a $LOG_FILE

    fi

    echo "$SCRIPT [`date "+%H:%M:%S"`]: calculating min/max elevation for color palette" | tee -a $LOG_FILE
    TOPO_MINMAX=`gmt grdinfo topo_cut.grd -C | awk '{print $6,$7}'`
    TOPO_MIN=`echo $TOPO_MINMAX | awk '{printf("%d"),$1}'`
    TOPO_MAX=`echo $TOPO_MINMAX | awk '{printf("%d"),$2}'`
    echo "$SCRIPT [`date "+%H:%M:%S"`]: minimum elevation: $TOPO_MIN" | tee -a $LOG_FILE
    echo "$SCRIPT [`date "+%H:%M:%S"`]: maximum elevation: $TOPO_MAX" | tee -a $LOG_FILE
    echo "$SCRIPT [`date "+%H:%M:%S"`]: using color mode $COLOR_MODE for topography cpt" | tee -a $LOG_FILE
    if [ "$COLOR_MODE" == "BW" ]
    then
        colortool -hue 300,180 -chroma 0,0 -lightness 75,100 -gmt -T${TOPO_MIN}/0/1 > make_seis_movie_topo.cpt
        colortool -hue 150,60 -chroma 0,0 -lightness 95,100 -gmt -T0/${TOPO_MAX}/1 >> make_seis_movie_topo.cpt
    elif [ "$COLOR_MODE" == "COLOR-LAND-ONLY" ]
    then
        colortool -hue 300,180 -chroma 0,0 -lightness 75,100 -gmt -T${TOPO_MIN}/0/1 > make_seis_movie_topo.cpt
        colortool -hue 150,60 -chroma 40,40 -lightness 60,100 -gmt -T0/${TOPO_MAX}/1 >> make_seis_movie_topo.cpt
    elif [ "$COLOR_MODE" == "COLOR-OCEAN-ONLY" ]
    then
        colortool -hue 300,180 -chroma 40,2 -lightness 45,100 -gmt -T${TOPO_MIN}/0/1 > make_seis_movie_topo.cpt
        colortool -hue 150,60 -chroma 0,0 -lightness 95,100 -gmt -T0/${TOPO_MAX}/1 >> make_seis_movie_topo.cpt
    elif [ "$COLOR_MODE" == "COLOR" ]
    then
        colortool -hue 300,180 -chroma 40,40 -lightness 45,100 -gmt -T${TOPO_MIN}/0/1 > make_seis_movie_topo.cpt
        colortool -hue 150,60 -chroma 40,40 -lightness 60,100 -gmt -T0/${TOPO_MAX}/1 >> make_seis_movie_topo.cpt
    else
        colortool -hue 300,180 -chroma 0,0 -lightness 75,100 -gmt -T${TOPO_MIN}/0/1 > make_seis_movie_topo.cpt
        colortool -hue 150,60 -chroma 0,0 -lightness 95,100 -gmt -T0/${TOPO_MAX}/1 >> make_seis_movie_topo.cpt
    fi
    tail -1 make_seis_movie_topo.cpt | awk '{print "F",$2}' >> make_seis_movie_topo.cpt
    head -1 make_seis_movie_topo.cpt | awk '{print "B",$2}' >> make_seis_movie_topo.cpt
    echo "$SCRIPT [`date "+%H:%M:%S"`]: plotting topography" | tee -a $LOG_FILE
    CMD="gmt grdimage topo_cut.grd $MAP_PROJ $MAP_LIMS -Itopo_cut_grad.grd -Cmake_seis_movie_topo.cpt -K -O"
    echo $CMD >> $LOG_FILE
    $CMD >> $PSFILE
fi



# Gridlines
if [ "$GRIDLINES" == "OFF" ]
then
    echo "$SCRIPT [`date "+%H:%M:%S"`]: plotting gridlines" | tee -a $LOG_FILE
    gmt psbasemap $MAP_PROJ $MAP_LIMS -Bxg -Byg -K -O --MAP_GRID_PEN=0.25p,65@60 >> $PSFILE
fi



# Coastline
echo "$SCRIPT [`date "+%H:%M:%S"`]: plotting coastline" | tee -a $LOG_FILE
CMD="gmt pscoast $MAP_PROJ $MAP_LIMS -Dh -W${PSCOAST_PEN} -C245 $PSCOAST_OPTIONS -K -O"
echo $CMD >> $LOG_FILE
$CMD >> $PSFILE



# Plate boundaries
if [ "$PLATE_BOUNDARY_FILE" != "" ]
then
    echo "$SCRIPT [`date "+%H:%M:%S"`]: plotting plate boundaries" | tee -a $LOG_FILE
    gmt psxy $PLATE_BOUNDARY_FILE $MAP_PROJ $MAP_LIMS -W1p -K -O >> $PSFILE
fi



# Other datasets
if [ "$PSXY_LIST" != "" ]
then
    echo $PSXY_LIST | awk -F";" '{for(i=2;i<=NF;i++){print $i}}' > make_seis_movie_psxy_list.tmp
    while read PSXY
    do
        PSXY_FILE=`echo $PSXY | awk -F: '{print $1}'`
        PSXY_OPTIONS=`echo $PSXY | awk -F: '{for(i=2;i<=NF;i++){print $i}}'`
        echo "$SCRIPT [`date "+%H:%M:%S"`]: plotting file $PSXY_FILE with psxy" | tee -a $LOG_FILE
        CMD="gmt psxy $MAP_PROJ $MAP_LIMS $PSXY_FILE $PSXY_OPTIONS -K -O"
        echo $CMD >> $LOG_FILE
        $CMD >> $PSFILE
    done < make_seis_movie_psxy_list.tmp
    rm make_seis_movie_psxy_list.tmp
fi


# Other text
if [ "$PSTEXT_LIST" != "" ]
then
    echo $PSTEXT_LIST | awk -F";" '{for(i=2;i<=NF;i++){print $i}}' > make_seis_movie_pstext_list.tmp
    while read PSTEXT
    do
        PSTEXT_FILE=`echo $PSTEXT | awk -F: '{print $1}'`
        PSTEXT_OPTIONS=`echo $PSTEXT | awk -F: '{for(i=2;i<=NF;i++){print $i}}'`
        PSTEXT_FILE_EXIST=`test -f "$PSTEXT_FILE" && echo Y || echo N`
        if [ "$PSTEXT_FILE_EXIST" == "Y" ]
        then
            echo "$SCRIPT [`date "+%H:%M:%S"`]: plotting file $PSTEXT_FILE with pstext" | tee -a $LOG_FILE
            CMD="gmt pstext $PSTEXT_FILE $MAP_PROJ $MAP_LIMS $PSTEXT_OPTIONS -K -O"
        else
            echo "$SCRIPT [`date "+%H:%M:%S"`]: plotting text $PSTEXT_FILE with pstext" | tee -a $LOG_FILE
            CMD="echo $PSTEXT_FILE | gmt pstext $MAP_PROJ $MAP_LIMS $PSTEXT_OPTIONS -K -O"
        fi
        echo $CMD >> $LOG_FILE
        eval $CMD >> $PSFILE
    done < make_seis_movie_pstext_list.tmp
    rm make_seis_movie_pstext_list.tmp
fi






# North arrow
if [ "$ADD_NORTH_ARROW" != "" ]
then
    case $ADD_NORTH_ARROW in
        TR) XSHFT=-0.1i; YSHFT=-0.1i;;
        *) XSHFT=0; YSHFT=0;;
    esac
    CMD="gmt psbasemap $MAP_PROJ $MAP_LIMS -Tdj${ADD_NORTH_ARROW} -Xa$XSHFT -Ya$YSHFT -K -O"
    $CMD >> $PSFILE
    echo $CMD >> $LOG_FILE
fi



# Basemap
echo "$SCRIPT [`date "+%H:%M:%S"`]: plotting basemap" | tee -a $LOG_FILE
CMD="gmt psbasemap $MAP_PROJ $MAP_LIMS -Bxa${MAP_TIKS} -Bya${MAP_TIKS} -BWeSn -K -O"
echo $CMD >> $LOG_FILE
$CMD >> $PSFILE







if [ "$PLOT_MAG_VS_TIME" == "Y" ]
then

    echo "$SCRIPT [`date "+%H:%M:%S"`]: plotting magnitude versus time frame" | tee -a $LOG_FILE

    # Initialize magnitude versus time origin
    gmt psxy -T -K -O -Y-2.5i >> $PSFILE

    # Magnitude versus time gridlines
    CMD="gmt psbasemap $TIME_PROJ $TIME_LIMS -Bsxg${MONTH_TIKS}o -Byg1 -K -O --MAP_GRID_PEN=0.25p,225,4_2:0"
    echo $CMD >> $LOG_FILE
    $CMD >> $PSFILE

    # Magnitude versus time frame
    NDAYS_INT=`echo $NDAYS | awk '{printf("%d"),$1}'`
    if [ $NDAYS_INT -le 7 ]
    then
        gmt psbasemap $TIME_PROJ $TIME_LIMS \
            -Bsxa${MONTH_TIKS}O+l"Date" -Bpxa${DAY_TIKS}D -Bya${MAG_TIKS}+l"Magnitude" \
            -BWeN -K -O \
            --MAP_TICK_LENGTH_PRIMARY=2.0p --MAP_TICK_LENGTH_SECONDARY=6.0p \
            --FORMAT_DATE_MAP=o-dd --FORMAT_TIME_MAP=a >> $PSFILE
    elif [ $NDAYS_INT -le 31 ]
    then
        gmt psbasemap $TIME_PROJ $TIME_LIMS \
            -Bsxa${MONTH_TIKS}O+l"Date" -Bpxa${DAY_TIKS}d -Bya${MAG_TIKS}+l"Magnitude" \
            -BWeN -K -O \
            --MAP_TICK_LENGTH_PRIMARY=2.0p --MAP_TICK_LENGTH_SECONDARY=6.0p \
            --FORMAT_DATE_MAP=o --FORMAT_TIME_MAP=a >> $PSFILE
    else
        gmt psbasemap $TIME_PROJ $TIME_LIMS \
            -Bsxa${YEAR_TIKS}Y+l"Date" -Bpxa${MONTH_TIKS}o -Bya${MAG_TIKS}+l"Magnitude" \
            -BWeN -K -O \
            --MAP_TICK_LENGTH_PRIMARY=2.0p --MAP_TICK_LENGTH_SECONDARY=6.0p >> $PSFILE
    fi
    if [ $NDAYS_INT -le 5 ]
    then
        gmt psbasemap $DAY_TIME_PROJ -R0/${DAY_END}/${MAG_AXIS_MIN}/${MAG_AXIS_MAX} \
            -Bxa${DAY_TIKS}+l"Number of Days Since ${DATE_START}" -BS -K -O >> $PSFILE
    else
        gmt psbasemap $DAY_TIME_PROJ -R0/${DAY_END}/${MAG_AXIS_MIN}/${MAG_AXIS_MAX} \
            -Bxa${DAY_TIKS}+l"Number of Days Since $(echo ${DATE_START} | awk -FT '{print $1}')" -BS -K -O >> $PSFILE
    fi

fi

# Finalize plot
gmt psxy -T -O >> $PSFILE

echo "$SCRIPT [`date "+%H:%M:%S"`]: converting background PostScript image to PNG" | tee -a $LOG_FILE
gmt psconvert -Tg -A frame_base.ps
rm frame_base.ps






####################################################################################################
#	CREATE REST OF ANIMATION
####################################################################################################




echo | tee -a $LOG_FILE


function generate_frames() {

    IFRAME=$1
    ENDFRAME=$2

    while [ $IFRAME -le $ENDFRAME ]
    do
        echo "$SCRIPT [`date "+%H:%M:%S"`]: Working on frame $IFRAME of $ENDFRAME" | tee -a $LOG_FILE

        # Postscript file name
        IFRAME_ZEROS=$(echo $IFRAME | awk '{printf("%05d"),$1}')
        PSFILE=frame_${IFRAME_ZEROS}.ps

        # Initialize frame
        gmt psxy -T -K -Y3.5i > $PSFILE


        # Date of frame
        if [ "$DDAY_NDAY_CHANGE_LIST" == "" ]
        then
            DATE=$(echo $DATE_START $IFRAME $DDAY | awk '{print $1,($2-1)*$3}' | dateutil -date -format "YYYY-MM-DDTHH:MM:SS")
        else
            #echo $IFRAME $IFRAME_CHANGE
            IDDAY=`echo $IFRAME $IFRAME_CHANGE |\
                awk '{
                    iframe = $1
                    idday = $2
                    for(i=NF; i>=2; i--) {
                        if (iframe > $i) {
                            idday = i-1
                            break
                        }
                    }
                    print idday
                }'`
            DDAY_CURRENT=`echo $IFRAME_CHANGE $DDAY_CHANGE |\
                awk '{
                    for (i=1;i<='$NCHANGE';i++) {
                        iframe_change[i] = $i
                        dday_array[i] = $(i+'$NCHANGE')
                        #print iframe_change[i],dday_array[i]
                    }
                } END {
                    dday = 0
                    i = 1
                    last_frame_change = 1
                    while (i<'$IDDAY') {
                        dday = dday + (iframe_change[i+1] - last_frame_change)*dday_array[i]
                        #print dday"..."
                        last_frame_change = iframe_change[i+1]
                        i++
                    }
                    dday = dday + ('$IFRAME'-last_frame_change)*dday_array[i]
                    print dday
                }'`
            DATE=$(echo $DATE_START $DDAY_CURRENT | dateutil -date -format "YYYY-MM-DDTHH:MM:SS")
            echo DDAY_CURRENT=$DDAY_CURRENT
            echo DATE=$DATE
        fi

        # Seismicity in time range of interest
        SEIS_FILE_FRAME=seis_range_${IFRAME}.tmp
        DAY=$(echo $DATE_START $DATE | dateutil -nday -format "YYYY-MM-DDTHH:MM:SS" | awk '{print $1}')
        paste nday.tmp make_seis_movie_seis.tmp |\
            awk '{
                if (0<=$1 && $1<='$DAY') {
                    trans = ('$DAY'-$1)/'$FADE_TIME'
                    trans = '$MAX_TRANS'*trans
                    if (trans>'$MAX_TRANS') {
                        trans = '$MAX_TRANS'
                    } else if (trans<0) {
                        trans = 0
                    }
                    print $1,$2,$3,$4,$6,trans,$7,$8,$9,$5
                }
            }' > $SEIS_FILE_FRAME # day date lon lat mag fade [str dip rak] dep

        # Seismicity
        awk '{
            if ($5>='$MAG_BOLD') {
                print "> -W2p"
            } else {
                print "> -W0.5p"
            }
            if ("'$COLOR_BY'"=="DATE") {
                print $3,$4,$1,$5*$5*$5*'$SEIS_SCALE',$6
            } else if ("'$COLOR_BY'"=="COL4" && NF==7) {
                print $3,$4,$7,$5*$5*$5*'$SEIS_SCALE',$6
            } else if ("'$COLOR_BY'"=="COL4" && NF==10) {
                print $3,$4,$10,$5*$5*$5*'$SEIS_SCALE',$6
            }
        }' $SEIS_FILE_FRAME |\
            gmt psxy $MAP_PROJ $MAP_LIMS -Sci -C${EQ_CPT} -t -K -O >> $PSFILE


        # Label large events
        awk '{if ($5>='$MAG_LABEL') {print $3,$4,$6,$5}}' $SEIS_FILE_FRAME |\
            gmt mapproject $MAP_PROJ $MAP_LIMS |\
            awk '{printf("%.3f %.3f %.3f %.1f\n"),$1,$2+2.54*$4*$4*$4*'$SEIS_SCALE'/2+0.08,$3,$4}' |\
            gmt mapproject $MAP_PROJ $MAP_LIMS -I |\
            awk '{printf("%.4f %.4f %.3f %.1f\n"),$1,$2,$3,$4}' |\
            gmt pstext $MAP_PROJ $MAP_LIMS -F+f12,1,black=0.5p,white+jCB -t -K -O >> $PSFILE


        # Focal mechanisms
        if [ "$PLOT_MECH" == "Y" ]
        then
            awk '{if($7!=0||$8!=0||$9!=0){print $3,$4,$1,$7,$8,$9,$5*$5*$5*0.001,$6/1.2}}' $SEIS_FILE_FRAME |\
                gmt psmeca $MAP_PROJ $MAP_LIMS -Sa5i -L0.5p -C${EQ_CPT} -t -K -O >> $PSFILE
        fi

        # Basemap
        gmt psbasemap $MAP_PROJ $MAP_LIMS -Bxa${MAP_TIKS} -Bya${MAP_TIKS} -BWeSn -K -O >> $PSFILE


        # Date in top left corner
        # echo $MAP_LIMS | sed -e "s/-R//" | awk -F/ '{print $1,$4,"12,3 LT '$DATE'"}' |\
        echo $DATE | gmt pstext $MAP_PROJ $MAP_LIMS -F+f12,3+cTL -D0.03i/-0.04i -Gwhite@15 -N -K -O >> $PSFILE


        if [ "$PSIMAGE_LIST" != "" ]
        then
            echo $PSIMAGE_LIST | awk -F";" '{for(i=2;i<=NF;i++){print $i}}' > psimage_list_${FRAME}.tmp
            while read PSIMAGE
            do
                PSIMAGE_FILE=`echo $PSIMAGE | awk -F: '{print $1}'`
                PSIMAGE_OPTIONS=`echo $PSIMAGE | awk -F: '{for(i=2;i<=NF;i++){print $i}}'`
                echo plotting file $PSIMAGE_FILE
                gmt psimage $MAP_PROJ $MAP_LIMS $PSIMAGE_FILE $PSIMAGE_OPTIONS -K -O >> $PSFILE
            done < psimage_list_${FRAME}.tmp
        fi



        if [ "$PLOT_MAG_VS_TIME" == "Y" ]
        then

            # Initialize magnitude versus time origin
            gmt psxy -T -K -O -Y-2.5i >> $PSFILE

            # Seismicity
            awk '{
                if ($5>='$MAG_BOLD') {
                    print "> -W2p"
                } else {
                    print "> -W0.5p"
                }
                if ("'$COLOR_BY'"=="DATE") {
                    print $2,$5,$1,$5*$5*$5*'"$SEIS_SCALE"',$6
                } else if ("'$COLOR_BY'"=="COL4" && NF==7) {
                    print $2,$5,$7,$5*$5*$5*'"$SEIS_SCALE"',$6
                }
            }' $SEIS_FILE_FRAME |\
                gmt psxy $TIME_PROJ $TIME_LIMS -Sci -C${EQ_CPT} -t -K -O >> $PSFILE

            # Focal mechanisms
            if [ "$PLOT_MECH" == "Y" ]
            then
                awk '{if($7!=0||$8!=0||$9!=0){print $2,$5,$1,$7,$8,$9,$5*$5*$5*0.001,$6/1.2}}' $SEIS_FILE_FRAME |\
                    gmt psmeca $TIME_PROJ $TIME_LIMS -Sa5i -L0.5p -C${EQ_CPT} -t -K -O >> $PSFILE
            fi

            # Date line
            echo $TIME_LIMS | sed -e "s/-R//" | awk -F/ '{print "'$DATE'",$3;print "'$DATE'",$4}' |\
                gmt psxy $TIME_PROJ $TIME_LIMS -W1p -K -O >> $PSFILE

            # Magnitude versus time frame
            NDAYS_INT=`echo $NDAYS | awk '{printf("%d"),$1}'`
            if [ $NDAYS_INT -le 7 ]
            then
                gmt psbasemap $TIME_PROJ $TIME_LIMS \
                    -Bsxa${MONTH_TIKS}O+l"Date" -Bpxa${DAY_TIKS}D -Bya${MAG_TIKS}+l"Magnitude" \
                    -BWeN -K -O \
                    --MAP_TICK_LENGTH_PRIMARY=2.0p --MAP_TICK_LENGTH_SECONDARY=6.0p \
                    --FORMAT_DATE_MAP=o-dd --FORMAT_TIME_MAP=a >> $PSFILE
            elif [ $NDAYS_INT -le 31 ]
            then
                gmt psbasemap $TIME_PROJ $TIME_LIMS \
                    -Bsxa${MONTH_TIKS}O+l"Date" -Bpxa${DAY_TIKS}d -Bya${MAG_TIKS}+l"Magnitude" \
                    -BWeN -K -O \
                    --MAP_TICK_LENGTH_PRIMARY=2.0p --MAP_TICK_LENGTH_SECONDARY=6.0p \
                    --FORMAT_DATE_MAP=o --FORMAT_TIME_MAP=a >> $PSFILE
            else
                gmt psbasemap $TIME_PROJ $TIME_LIMS \
                    -Bsxa${YEAR_TIKS}Y+l"Date" -Bpxa${MONTH_TIKS}o -Bya${MAG_TIKS}+l"Magnitude" \
                    -BWeN -K -O \
                    --MAP_TICK_LENGTH_PRIMARY=2.0p --MAP_TICK_LENGTH_SECONDARY=6.0p >> $PSFILE
            fi
            if [ $NDAYS_INT -le 5 ]
            then
                gmt psbasemap $DAY_TIME_PROJ -R0/${DAY_END}/${MAG_AXIS_MIN}/${MAG_AXIS_MAX} \
                    -Bxa${DAY_TIKS}+l"Number of Days Since ${DATE_START}" -BS -K -O >> $PSFILE
            else
                gmt psbasemap $DAY_TIME_PROJ -R0/${DAY_END}/${MAG_AXIS_MIN}/${MAG_AXIS_MAX} \
                    -Bxa${DAY_TIKS}+l"Number of Days Since $(echo ${DATE_START} | awk -FT '{print $1}')" -BS -K -O >> $PSFILE
            fi

        fi

        # Finalize plot
        gmt psxy -T -O >> $PSFILE


        # Convert frame to PNG and superimpose on fixed background image
        PNG=`basename $PSFILE .ps`.png
        gmt psconvert $PSFILE -TG -A
        gm composite $PNG frame_base.png j${IFRAME} && mv j${IFRAME} $PNG


        # Remove seismicity file for this frame
        rm -f $SEIS_FILE_FRAME

        # Record frame as completed to track progress
        echo $IFRAME >> frame_progress.tmp

        # Update frame counter
        IFRAME=$(echo $IFRAME | awk '{print $1+1}')

    done

    echo done >> done_file.tmp

}


# Create frames in parallel
# Each process takes a chunk of the loop; set the start and end frames for each sub-loop
NFRAMES_LOCAL=`echo $NFRAMES $NPROC | awk '{print $1/$2}'`
seq $NFRAMES_LOCAL $NFRAMES_LOCAL $NFRAMES |\
    awk 'BEGIN{start=1}{
        if (NR=='$NPROC') {
            end = '$NFRAMES'
        } else {
            end = int($0)
        }
        print start,end
        start = end + 1
    }' > end_frame_list.tmp

# Remove processor count and total frame progress file if they exist
test -f done_file.tmp && rm -f done_file.tmp
test -f frame_progress.tmp && rm -f frame_progress.tmp
touch done_file.tmp
touch frame_progress.tmp

# Generate frames in parallel
while read START END
do
    generate_frames $START $END &
done < end_frame_list.tmp

# Have we finished yet?
DONE=`test -f done_file.tmp && wc done_file.tmp | awk '{print $1}' || echo 0`
while [ $DONE -lt $NPROC ]
do
    FRAMES_COMPLETED=$(wc frame_progress.tmp | awk '{print $1}')
    echo "$SCRIPT [`date "+%H:%M:%S"`]: finished $FRAMES_COMPLETED frames out of $NFRAMES total" | tee -a $LOG_FILE
    DONE=`wc done_file.tmp | awk '{print $1}'`
    sleep 10s
done




# Clean up a little bit
rm frame_*.ps
rm make_seis_movie_seis.tmp nday.tmp color_list.tmp ${EQ_CPT}
rm make_seis_movie.sh.*.cpt
if [ "$TOPO_CLEAN" == "ON" ]
then
    rm topo_cut.grd topo_cut_grad.grd
fi
rm topo_bw.cpt
rm end_frame_list.tmp done_file.tmp frame_progress.tmp




# Generate movie from frames
echo | tee -a $LOG_FILE
echo "$SCRIPT [`date "+%H:%M:%S"`]: running ffmpeg to create movie \"${OUTPUT_MP4_FILE}\"" | tee -a $LOG_FILE
echo ffmpeg -loglevel warning -framerate 24 -y -i "frame_%05d.png" -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -pix_fmt yuv420p seis_movie.mp4 >> $LOG_FILE
ffmpeg -loglevel warning -framerate 24 -y -i "frame_%05d.png" -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -pix_fmt yuv420p ${OUTPUT_MP4_FILE}


if [ "$CLEAN" == "Y" ]
then
    rm frame_*.png
fi



echo "$SCRIPT [`date "+%H:%M:%S"`]: see messages in log file $LOG_FILE"
echo "$SCRIPT [`date "+%H:%M:%S"`]: finished" | tee -a $LOG_FILE