#!/bin/bash

#####
#
#####
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
    echo "-Rw/e/s/n                     Map limits (default: plot all seismicity)" 1>&2
    echo "-mag:min MAG_MIN              Minimum magnitude for magnitude versus time frame" 1>&2
    echo "-mag:max MAG_MAX              Maximum magnitude for magnitude versus time frame" 1>&2
    echo "-mag:color MAG_COLOR_CHANGE   Change color every time earthquake above MAG_COLOR_CHANGE occurs" 1>&2
    echo "-mag:bold MAG_BOLD            Make earthquakes larger than MAG_BOLD bolder" 1>&2
    echo "-dday:change DDAY,NDAYS       Change DDAY at NDAYS from start" 1>&2
    echo "-timecpt TIME_CPT             Earthquake timing (by day) color palette" 1>&2
    echo "-fade FADE_TIME               Time to fade out symbols (days)" 1>&2
    echo "-trans:max MAX_TRANS          Maximum transparency of faded symbols (default:90)" 1>&2
    echo "-mech                         Plot focal mechanisms (default: no mechanisms)" 1>&2
    echo "-plates PB_FILE               Plot plate boundaries (default: no boundaries)" 1>&2
    echo "-topo TOPO_FILE               Plot shaded topography in background" 1>&2
    echo "-other:psxy FILE:OPT          Other psxy file to plot (can repeat)" 1>&2
    echo "-other:pstext WORDS:OPT       Other psxy file to plot (can repeat)" 1>&2
    echo "-other:psimage FILE:OPT       Other psimage file to plot (can repeat)" 1>&2
    echo "-clean                        Remove frame files after running (default is to keep them)" 1>&2
    exit 1
}

if [ $# -eq 0 ]
then
    usage
fi

# Seismicity file: origin_time(YYYY-MM-DDTHH:MM:SS)  lon  lat  dep  mag  str  dip  rak
SEIS_FILE=
SEIS_FILE=seis.dat

# Animation parameters
DATE_START=                             # Starting date
DATE_END=                               # Ending date
DDAY=                                   # Number of days between frames
FADE_TIME=10000000                      # Time to fade earthquakes out (days)
MAX_TRANS=90
DDAY_NDAY_CHANGE_LIST=

# Map parameters
MAP_LIMS=
MAG_MIN=2
MAG_MAX=8
MAG_COLOR_CHANGE=
MAG_BOLD=15
TIME_CPT=
PLOT_MECH=N
PLATE_BOUNDARY_FILE=
TOPO_FILE=

# Other commands
PSXY_LIST=""
PSTEXT_LIST=""
PSIMAGE_LIST=""
CLEAN="N"


while [ "$1" != "" ]
do
    case $1 in
        -f) shift; SEIS_FILE=$1;;
        -start) shift; DATE_START=$1;;
        -end) shift; DATE_END=$1;;
        -dday) shift; DDAY=$1;;
        -dday:change) shift; DDAY_NDAY_CHANGE_LIST="$DDAY_NDAY_CHANGE_LIST $1";;
        -R*) MAP_LIMS=$1;;
        -mag:min) shift; MAG_MIN=$1;;
        -mag:max) shift; MAG_MAX=$1;;
        -mag:color) shift; MAG_COLOR_CHANGE=$1;;
        -mag:bold) shift; MAG_BOLD=$1;;
        -timecpt) shift; TIME_CPT=$1;;
        -fade) shift; FADE_TIME=$1;;
        -trans:max) shift; MAX_TRANS=$1;;
        -mech) PLOT_MECH=Y;;
        -plates) shift; PLATE_BOUNDARY_FILE=$1;;
        -topo) shift; TOPO_FILE=$1;;
        -other:psxy) shift;PSXY_LIST="$PSXY_LIST;$1";;
        -other:pstext) shift;PSTEXT_LIST="$PSTEXT_LIST;$1";;
        -other:psimage) shift;PSIMAGE_LIST="$PSIMAGE_LIST;$1";;
        -clean) CLEAN="Y";;
        *) ;;
    esac
    shift
done

# Check that file variables are defined and files exist
for VAR in SEIS_FILE
do
    if [ "${!VAR}" == "" ]
    then
        echo "$0: File $VAR is not defined" 1>&2
        usage
    fi
    if [ ! -f ${!VAR} ]
    then
        echo "$0: Could not find $VAR named \"${!VAR}\"" 1>&2
        usage
    fi
done

# If MAP_LIMS is not defined, set it here from SEIS_FILE
if [ "$MAP_LIMS" == "" ]
then
    echo "$0: getting map limits from the input file"
    MAP_LIMS=$(gmt gmtinfo $SEIS_FILE -I0.01 -i1:2)
    echo "$0: using map limits $MAP_LIMS"
fi

# If plate boundary file is defined but does not exist, set to none
if [ "$PLATE_BOUNDARY_FILE" != "" ]
then
    if [ ! -f $PLATE_BOUNDARY_FILE ]
    then
        echo "$0: could not find plate boundary file named \"$PLATE_BOUNDARY_FILE\""
        echo "$0: not plotting plate boundaries"
        PLATE_BOUNDARY_FILE=
    fi
fi



# Check that variables are defined
for VAR in DATE_START DATE_END DDAY MAP_LIMS
do
    if [ "${!VAR}" == "" ]
    then
        echo "$0: Variable $VAR is not defined" 1>&2
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

# Determine magnitude for color changes
if [ "$MAG_COLOR_CHANGE" == "" -a "$TIME_CPT" == "" ]
then
    echo "$0: determining magnitude to change colors"
    MAG_COLOR_CHANGE=$(gmt gmtinfo $SEIS_FILE -i4 -C | awk '{print $2-1.0}')
    echo "$0: changing colors every magnitude ${MAG_COLOR_CHANGE}+"
elif [ "$TIME_CPT" != "" ]
then
    if [ ! -f $TIME_CPT ]
    then
        echo "$0: could not find timing color palette named \"$TIME_CPT\""
        usage
    fi
    echo "$0: using custom timing color palette \"$TIME_CPT\""
fi


#####
#	SET ANIMATION PARAMETERS
#####
NDAYS=$(echo $DATE_START $DATE_END | dateutil -nday -format "YYYY-MM-DDTHH:MM:SS")
if [ "$DDAY_NDAY_CHANGE_LIST" == "" ]
then
    NFRAMES=$(echo $NDAYS $DDAY | awk '{printf("%d"),$1/$2}')
    IFRAME_CHANGE=1
    DDAY_CHANGE=$DDAY
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

gmt set PS_MEDIA 11ix17i
gmt set FORMAT_GEO_OUT D
gmt set FORMAT_GEO_MAP D

MAP_PROJ="-JM5i"
MAP_TIKS=$(echo $MAP_LIMS | sed -e "s/-R//" |\
           awk -F/ '{
               dlon = $2-$1
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

TIME_PROJ="-JX5iT/1.5i"
DAY_TIME_PROJ=$(echo $TIME_PROJ | sed -e "s/T//")
TIME_LIMS="-R${DATE_START}/${DATE_END}/${MAG_MIN}/${MAG_MAX}"
DAY_END=$(echo $DATE_START $DATE_END | dateutil -nday -format "YYYY-MM-DDTHH:MM:SS" | awk '{print $1}')
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
DAY_TIKS=$(echo $NDAYS |\
           awk '{
               if ($1<2) {
                   print 0.5
               } else if ($1<5) {
                   print 1
               } else if ($1<10) {
                   print 2
               } else if ($1<20) {
                   print 5
               } else if ($1<50) {
                   print 10
               } else if ($1<100) {
                   print 20
               } else if ($1<200) {
                   print 50
               } else if ($1<500) {
                   print 100
               } else if ($1<1000) {
                   print 200
               } else if ($1<2000) {
                   print 500
               } else if ($1<5000) {
                   print 1000
               } else if ($1<10000) {
                   print 2000
               } else if ($1<20000) {
                   print 5000
               } else {
                   print 10000
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

# Select only data inside map frame
gmt gmtselect $SEIS_FILE $MAP_LIMS -i1:2,0,3:7 -o2,0:1,3:7 > make_seis_movie_seis.tmp

# Calculate number of days since starting time
awk '{print "'$DATE_START'",$1}' make_seis_movie_seis.tmp |\
    dateutil -nday -format "YYYY-MM-DDTHH:MM:SS" > nday.tmp

# Color palette
if [ "$TIME_CPT" == "" ]
then
    paste nday.tmp make_seis_movie_seis.tmp |\
        awk '{if($1>0 && $6>='$MAG_COLOR_CHANGE'){printf("%.5f\n"),$1}}' > make_seis_movie_nday_color.tmp
    NCOLORS=$(wc make_seis_movie_nday_color.tmp | awk '{print $1+1}')
    colortool -hue -20,90 -lightness 30,95 -chroma 100,100 -ncolors $NCOLORS > color_list.tmp
    awk '{
        if (NR==1) {
            print 0,$1
        } else {
            print nday,$1
        }
        nday = $1
    }' make_seis_movie_nday_color.tmp |\
        paste - color_list.tmp |\
        awk '{if(NF>1){print $1,$3,$2,$3}}END{print "B black";print "F",$1}' > make_seis_movie_time.cpt
    rm make_seis_movie_nday_color.tmp
    echo Built timing color palette make_seis_movie_time.cpt:
    cat make_seis_movie_time.cpt
else
    cp $TIME_CPT make_seis_movie_time.cpt
fi

#####
#	CREATE ANIMATION
#####
IFRAME=1
while [ $IFRAME -le $NFRAMES ]
do
    echo "Working on frame $IFRAME of $NFRAMES"

    # Postscript file name
    IFRAME_ZEROS=$(echo $IFRAME | awk '{printf("%05d"),$1}')
    PSFILE=frame_${IFRAME_ZEROS}.ps

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
                print $1,$2,$3,$4,$6,trans,$7,$8,$9
            }
        }' > seis_range.tmp # day date lon lat mag fade

    # Initialize plot
    gmt psxy -T -K -Y3.5i > $PSFILE

    # Topography
    if [ "$TOPO_FILE" != "" ]
    then
        if [ ! -f topo_cut_grad.grd ]
        then
            gmt grdcut $SRTM15 $MAP_LIMS -Gtopo_cut.grd
            gmt grdgradient topo_cut.grd -A0/270 -Ne0.4 -Gtopo_cut_grad.grd
        fi
        if [ ! -f topo_bw.cpt ]
        then
            TOPO_MINMAX=`gmt grdinfo topo_cut.grd -C | awk '{print $6,$7}'`
            TOPO_MIN=`echo $TOPO_MINMAX | awk '{printf("%d"),$1}'`
            TOPO_MAX=`echo $TOPO_MINMAX | awk '{printf("%d"),$2}'`
            colortool -hue 300,180 -chroma 0,0 -lightness 75,100 -gmt -T${TOPO_MIN}/0/1 > topo_bw.cpt
            colortool -hue 150,60 -chroma 0,0 -lightness 95,100 -gmt -T0/${TOPO_MAX}/1 >> topo_bw.cpt
            tail -1 topo_bw.cpt | awk '{print "F",$2}' >> topo_bw.cpt
            head -1 topo_bw.cpt | awk '{print "B",$2}' >> topo_bw.cpt
        fi
        gmt grdimage topo_cut.grd $MAP_PROJ $MAP_LIMS -Itopo_cut_grad.grd -Ctopo_bw.cpt -K -O >> $PSFILE
    fi

    # Coastline
    gmt pscoast $MAP_PROJ $MAP_LIMS -Dh -W1p -C245 -K -O >> $PSFILE

    # Plate boundaries
    if [ "$PLATE_BOUNDARY_FILE" != "" ]
    then
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
            echo plotting file $PSXY_FILE
            #echo $PSXY_FILE
            #echo $PSXY_OPTIONS
            gmt psxy $MAP_PROJ $MAP_LIMS $PSXY_FILE $PSXY_OPTIONS -K -O >> $PSFILE
        done < make_seis_movie_psxy_list.tmp
        rm make_seis_movie_psxy_list.tmp
    fi

    # Seismicity
    awk '{
        if ($5>='$MAG_BOLD') {
            print "> -W2p"
        } else {
            print "> -W0.5p"
        }
        print $3,$4,$1,$5*$5*$5*0.001,$6
    }' seis_range.tmp |\
        gmt psxy $MAP_PROJ $MAP_LIMS -Sci -Cmake_seis_movie_time.cpt -t -K -O >> $PSFILE

    # Focal mechanisms
    if [ "$PLOT_MECH" == "Y" ]
    then
        awk '{if($7!=0||$8!=0||$9!=0){print $3,$4,$1,$7,$8,$9,$5*$5*$5*0.001,$6/1.2}}' seis_range.tmp |\
            gmt psmeca $MAP_PROJ $MAP_LIMS -Sa5i -L0.5p -Cmake_seis_movie_time.cpt -t -K -O >> $PSFILE
    fi

    # Basemap
    gmt psbasemap $MAP_PROJ $MAP_LIMS -Bxa${MAP_TIKS} -Bya${MAP_TIKS} -BWeSn -K -O >> $PSFILE

    # Date in bottom left corner
    echo $MAP_LIMS | sed -e "s/-R//" | awk -F/ '{print $1,$4,"10,2 LT '$DATE'"}' |\
        gmt pstext $MAP_PROJ $MAP_LIMS -F+f+j -D0.05i/-0.05i -N -K -O >> $PSFILE

    # Other text
    if [ "$PSTEXT_LIST" != "" ]
    then
        echo $PSTEXT_LIST | awk -F";" '{for(i=2;i<=NF;i++){print $i}}' > make_seis_movie_pstext_list.tmp
        while read PSTEXT
        do
            PSTEXT_WORDS=`echo $PSTEXT | awk -F: '{print $1}'`
            PSTEXT_OPTIONS=`echo $PSTEXT | awk -F: '{for(i=2;i<=NF;i++){print $i}}'`
            echo plotting words $PSTEXT_WORDS
            #echo $PSTEXT_FILE
            #echo $PSTEXT_OPTIONS
            echo $PSTEXT_WORDS | gmt pstext $MAP_PROJ $MAP_LIMS $PSTEXT_OPTIONS -K -O >> $PSFILE
        done < make_seis_movie_pstext_list.tmp
        rm make_seis_movie_pstext_list.tmp
    fi



    if [ "$PSIMAGE_LIST" != "" ]
    then
        echo $PSIMAGE_LIST | awk -F";" '{for(i=2;i<=NF;i++){print $i}}' > psimage_list.tmp
        while read PSIMAGE
        do
            PSIMAGE_FILE=`echo $PSIMAGE | awk -F: '{print $1}'`
            PSIMAGE_OPTIONS=`echo $PSIMAGE | awk -F: '{for(i=2;i<=NF;i++){print $i}}'`
            echo plotting file $PSIMAGE_FILE
            #echo $PSIMAGE_FILE
            #echo $PSIMAGE_OPTIONS
            gmt psimage $MAP_PROJ $MAP_LIMS $PSIMAGE_FILE $PSIMAGE_OPTIONS -K -O >> $PSFILE
        done < psimage_list.tmp
    fi



    # Initialize magnitude versus time origin
    gmt psxy -T -K -O -Y-2.5i >> $PSFILE

    # Magnitude versus time gridlines
    gmt psbasemap $TIME_PROJ $TIME_LIMS -Bsxg${MONTH_TIKS}o -Byg1 -K -O --MAP_GRID_PEN=0.25p,225,4_2:0 >> $PSFILE

    # Seismicity
    awk '{
        if ($5>='$MAG_BOLD') {
            print "> -W2p"
        } else {
            print "> -W0.5p"
        }
        print $2,$5,$1,$5*$5*$5*0.001,$6
    }' seis_range.tmp |\
        gmt psxy $TIME_PROJ $TIME_LIMS -Sci -Cmake_seis_movie_time.cpt -t -K -O >> $PSFILE

    # Focal mechanisms
    if [ "$PLOT_MECH" == "Y" ]
    then
        awk '{if($7!=0||$8!=0||$9!=0){print $2,$5,$1,$7,$8,$9,$5*$5*$5*0.001,$6/1.2}}' seis_range.tmp |\
            gmt psmeca $TIME_PROJ $TIME_LIMS -Sa5i -L0.5p -Cmake_seis_movie_time.cpt -t -K -O >> $PSFILE
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
            --MAP_TICK_LENGTH_PRIMARY=2.0p --MAP_TICK_LENGTH_SECONDARY=6.0p >> $PSFILE
    else
        gmt psbasemap $TIME_PROJ $TIME_LIMS \
            -Bsxa${YEAR_TIKS}Y+l"Date" -Bpxa${MONTH_TIKS}o -Bya${MAG_TIKS}+l"Magnitude" \
            -BWeN -K -O \
            --MAP_TICK_LENGTH_PRIMARY=2.0p --MAP_TICK_LENGTH_SECONDARY=6.0p >> $PSFILE
    fi
    if [ $NDAYS_INT -le 5 ]
    then
        gmt psbasemap $DAY_TIME_PROJ -R0/${DAY_END}/${MAG_MIN}/${MAG_MAX} \
            -Bxa${DAY_TIKS}+l"Number of Days Since ${DATE_START}" -BS -K -O >> $PSFILE
    else
        gmt psbasemap $DAY_TIME_PROJ -R0/${DAY_END}/${MAG_MIN}/${MAG_MAX} \
            -Bxa${DAY_TIKS}+l"Number of Days Since $(echo ${DATE_START} | awk -FT '{print $1}')" -BS -K -O >> $PSFILE
    fi

    # Finalize plot
    gmt psxy -T -O >> $PSFILE
    IFRAME=$(echo $IFRAME | awk '{print $1+1}')
done

echo "$0: converting PostScript frames to PNG"
gmt psconvert -Tg -A frame_*.ps -Vt
rm frame_*.ps
rm make_seis_movie_seis.tmp nday.tmp color_list.tmp make_seis_movie_time.cpt seis_range.tmp
rm topo_cut.grd topo_cut_grad.grd topo_bw.cpt
if [ "$CLEAN" == "Y" ]
then
    rm frame_*.png
fi

echo "$0: creating movie \"seis_movie.mp4\""
ffmpeg -loglevel warning -framerate 24 -y -i "frame_%05d.png" -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -pix_fmt yuv420p seis_movie.mp4
