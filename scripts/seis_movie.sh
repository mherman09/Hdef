#!/bin/bash

function usage {
    echo "Usage: $0 ...options..." 1>&2
    echo "-f FILE        Seismicity data (YYYY-MM-DDTHH:MM:SS LON LAT DEP MAG [STR DIP RAK])" 1>&2
    echo "-tstart        Starting time for movie (YYYY-MM-DDTHH:MM:SS)" 1>&2
    echo "-tend          Ending time for movie (YYYY-MM-DDTHH:MM:SS)" 1>&2
    echo "-dt            Time between frames (days)" 1>&2
    echo "-fade FADE_T   Time to fade events to maximum transparency" 1>&2
    echo "-nfade NFADE   Number of fading increments" 1>&2
    echo "-max_trans TR  Maximum transparency" 1>&2
    echo "-mech          Plot focal mechanisms" 1>&2
    echo "-color DEP|MEC Color by depth or focal mechanism (1.5=th, 2.5=ss, 3.5=no, 4.5=other)" 1>&2
    echo "-cpt CPTFILE   Color palette file" 1>&2
    exit 1
}

SEIS_FILE=""
TSTART="2010-01-01T00:00:00"
TEND="2017-01-01T00:00:00"
DT="200"
FADE_DURATION="-1"
NFADE="3"
MAX_TRANS="80"
PLOT_MECH="N"
COLOR_BY="DEP"
CPT_FILE=""
while [ "$1" != "" ]
do
    case $1 in
        -f) shift; SEIS_FILE="$1";;
        -tstart) shift; TSTART="$1";;
        -tend) shift; TEND="$1";;
        -dt) shift; DT="$1";;
        -fade) shift; FADE_DURATION="$1";;
        -nfade) shift; NFADE="$1";;
        -max_trans) shift; MAX_TRANS="$1";;
        -mech) PLOT_MECH="Y";;
        -color) shift;COLOR_BY="$1";;
        -cpt) shift;CPT_FILE="$1";;
        *) echo "No option $1" 1>&2; usage;;
    esac
    shift
done

if [ "$SEIS_FILE" == "" ]; then echo "No seismicity file defined" 1>&2; usage; fi
if [ ! -f "$SEIS_FILE" ]; then echo "No seismicity file found named \"$SEIS_FILE\"" 1>&2; usage; fi


awk '{print "'"$TSTART"'",$1}' $SEIS_FILE | dateutil -ndays -long > nday.tmp
if [ $COLOR_BY == "MEC" ]
then
    awk '{print $6,$7,$8}' $SEIS_FILE | mtutil -sdr -ternary ternary.tmp
fi

TOTAL_DAYS=`dateutil -c $TSTART $TEND -ndays -long`
NFRAMES=`grid -x $DT $TOTAL_DAYS -dx $DT | wc | awk '{print $1}'`

DO_FADE=`echo $FADE_DURATION | awk '{if($1<=0){print "N"}else{print "Y"}}'`
FADE_DT=`echo $FADE_DURATION $NFADE | awk '{print $1/($2-1)}'`

#####
#	GENERATE FRAMES
#####
gmt set PS_MEDIA 8.5ix11i
gmt set MAP_FRAME_TYPE plain

PROJ="-JM6i -P"
LIMS="-R137/146/34/42"

if [ $COLOR_BY == "DEP" ]
then
    if [ -f "$CPT_FILE" ]
    then
        cp $CPT_FILE dep.cpt
    else
        gmt makecpt -Cmagma.cpt -I -T0/100/1 -D > dep.cpt
    fi
elif [ $COLOR_BY == "MEC" ]
then
    if [ -f "$CPT_FILE" ]
    then
        cp $CPT_FILE mech.cpt
    else
        echo 1 255/155/0 2 255/155/0 > mech.cpt
        echo 2 0/155/0 3 0/155/0 >> mech.cpt
        echo 3 0/0/55 4 0/0/55 >> mech.cpt
        echo 4 255/205/205 5 255/205/205 >> mech.cpt
    fi
else
    usage
fi

IFRAME="00001"
while [ $IFRAME -le $NFRAMES ]
do
    PSFILE="seis_movie_${IFRAME}.ps"
    DAY=`echo $IFRAME $DT | awk '{print $1*$2}'`

    echo IFRAME=$IFRAME DAY=$DAY

    echo 0 0 | gmt psxy $PROJ $LIMS -K > $PSFILE
    gmt pscoast $PROJ $LIMS -Dh -W1p -K -O >> $PSFILE
    gmt psxy $USGS15 $PROJ $LIMS -W1p -K -O >> $PSFILE

    if [ "$DO_FADE" == "Y" ]
    then

        # OLDER EVENTS FADE OVER TIME

        # PLOT OLDEST EVENTS WITH MAXIMUM TRANSPARENCY
        IFADE=$NFADE
        TRANS=`echo $IFADE $NFADE $MAX_TRANS | awk '{print ($1-1)/($2-1)*$3}'`
        echo IFADE=$IFADE NFADE=$NFADE MAX_TRANS=$MAX_TRANS TRANS=$TRANS
        if [ $PLOT_MECH == "N" ]
        then
            paste nday.tmp $SEIS_FILE |\
                awk 'BEGIN{day='"$DAY"';fade_duration='"$FADE_DURATION"'}{
                    if($1>=0 && $1<=day-fade_duration){
                        print $3,$4,$5,$6*$6*$6*0.001
                    }
                }' |\
                gmt psxy $PROJ $LIMS -Sci -Cdep.cpt -W0.5p -K -O -t$TRANS >> $PSFILE
        else
            if [ $COLOR_BY == "DEP" ]
            then
                paste nday.tmp $SEIS_FILE |\
                    awk 'BEGIN{day='"$DAY"';fade_duration='"$FADE_DURATION"'}{
                        if($1>=0 && $1<=day-fade_duration){
                            print $3,$4,$5,$7,$8,$9,$6*$6*$6*0.001
                        }
                    }' |\
                    gmt psmeca $PROJ $LIMS -Sa5i -Zdep.cpt -W0.5p -K -O -t$TRANS >> $PSFILE
            elif [ $COLOR_BY == "MEC" ]
            then
                paste nday.tmp $SEIS_FILE ternary.tmp |\
                    awk 'BEGIN{day='"$DAY"';fade_duration='"$FADE_DURATION"'}{
                        if($1>=0 && $1<=day-fade_duration){
                            if ($13>=0.55) {
                                mech = 1.5
                            } else if ($14>=0.75) {
                                mech = 2.5
                            } else if ($15>=0.60) {
                                mech = 3.5
                            } else {
                                mech = 4.5
                            }
                            print $3,$4,mech,$7,$8,$9,$6*$6*$6*0.001
                        }
                    }' |\
                    gmt psmeca $PROJ $LIMS -Sa5i -Zmech.cpt -W0.5p -K -O -t$TRANS >> $PSFILE
            fi
        fi

        # PLOT OTHER EVENTS WITH TRANSPARENCY INCREASING WITH AGE
        IFADE=`echo $IFADE | awk '{print $1-1}'`
        while [ $IFADE -ge 1 ]
        do
            TRANS=`echo $IFADE $NFADE $MAX_TRANS | awk '{print ($1-1)/($2-1)*$3}'`
            echo IFADE=$IFADE NFADE=$NFADE MAX_TRANS=$MAX_TRANS TRANS=$TRANS
            if [ $PLOT_MECH == "N" ]
            then
                paste nday.tmp $SEIS_FILE |\
                    awk 'BEGIN{
                        day='"$DAY"'
                        fade_duration='"$FADE_DURATION"'
                        nfade = '"$NFADE"'
                        ifade = '"$IFADE"'
                        tmin = day - fade_duration*(ifade)/(nfade-1)
                        tmax = day - fade_duration*(ifade-1)/(nfade-1)
                    }{
                        if ($1>=0 && tmin<$1 && $1<=tmax){ 
                            print $3,$4,$5,$6*$6*$6*0.001
                        }
                    }' |\
                    gmt psxy $PROJ $LIMS -Sci -Cdep.cpt -W0.5p -K -O -t$TRANS >> $PSFILE
            else
                if [ $COLOR_BY == "DEP" ]
                then
                    paste nday.tmp $SEIS_FILE |\
                        awk 'BEGIN{
                            day='"$DAY"'
                            fade_duration='"$FADE_DURATION"'
                            nfade = '"$NFADE"'
                            ifade = '"$IFADE"'
                            tmin = day - fade_duration*(ifade)/(nfade-1)
                            tmax = day - fade_duration*(ifade-1)/(nfade-1)
                        }{
                            if ($1>=0 && tmin<$1 && $1<=tmax){ 
                                print $3,$4,$5,$7,$8,$9,$6*$6*$6*0.001
                            }
                        }' |\
                        gmt psmeca $PROJ $LIMS -Sa5i -Zdep.cpt -W0.5p -K -O -t$TRANS >> $PSFILE
                elif [ $COLOR_BY == "MEC" ]
                then
                    paste nday.tmp $SEIS_FILE ternary.tmp |\
                        awk 'BEGIN{
                            day='"$DAY"'
                            fade_duration='"$FADE_DURATION"'
                            nfade = '"$NFADE"'
                            ifade = '"$IFADE"'
                            tmin = day - fade_duration*(ifade)/(nfade-1)
                            tmax = day - fade_duration*(ifade-1)/(nfade-1)
                        }{
                            if ($13>=0.55) {
                                mech = 1.5
                            } else if ($14>=0.75) {
                                mech = 2.5
                            } else if ($15>=0.60) {
                                mech = 3.5
                            } else {
                                mech = 4.5
                            }
                            if ($1>=0 && tmin<$1 && $1<=tmax){ 
                                print $3,$4,mech,$7,$8,$9,$6*$6*$6*0.001
                            }
                        }' |\
                        gmt psmeca $PROJ $LIMS -Sa5i -Zmech.cpt -W0.5p -K -O -t$TRANS >> $PSFILE
                fi
            fi
            IFADE=`echo $IFADE | awk '{print $1-1}'`
        done

    else

        # EVENTS REMAIN OPAQUE OVER TIME

        if [ $PLOT_MECH == "N" ]
        then
            paste nday.tmp $SEIS_FILE |\
                awk 'BEGIN{day='"$DAY"'}{
                    if ($1>=0 && $1<=day){ 
                        print $3,$4,$5,$6*$6*$6*0.0006
                    }
                }' |\
                gmt psxy $PROJ $LIMS -Sci -Cdep.cpt -W0.5p -K -O >> $PSFILE
        else
            if [ $COLOR_BY == "DEP" ]
            then
                paste nday.tmp $SEIS_FILE |\
                    awk 'BEGIN{day='"$DAY"'}{
                        if ($1>=0 && $1<=day){ 
                            print $3,$4,$5,$7,$8,$9,$6*$6*$6*0.0006
                        }
                    }' |\
                    gmt psmeca $PROJ $LIMS -Sa5i -Zdep.cpt -W0.5p -K -O >> $PSFILE
            elif [ $COLOR_BY == "MEC" ]
            then
                paste nday.tmp $SEIS_FILE ternary.tmp |\
                    awk 'BEGIN{day='"$DAY"'}{
                        if ($1>=0 && $1<=day){ 
                            if ($13>=0.55) {
                                mech = 1.5
                            } else if ($14>=0.75) {
                                mech = 2.5
                            } else if ($15>=0.60) {
                                mech = 3.5
                            } else {
                                mech = 4.5
                            }
                            print $3,$4,mech,$7,$8,$9,$6*$6*$6*0.0006
                        }
                    }' |\
                    gmt psmeca $PROJ $LIMS -Sa5i -Zmech.cpt -K -O >> $PSFILE
            fi
        fi
    fi

    gmt psbasemap $PROJ $LIMS -Bxa1 -Bya1 -K -O >> $PSFILE
    dateutil -c $TSTART $DAY -date -long |\
        awk '{printf("0.05 0.05 20,0 LB %04d-%02d-%02d %02d:00:00\n"),$1,$2,$3,$4}' |\
        gmt pstext -JX1i -R0/1/0/1 -F+f+j -Gwhite@20 -N -K -O >> $PSFILE
    echo 0 0 | gmt psxy $PROJ $LIMS -O >> $PSFILE

    gmt psconvert $PSFILE -Tj
    rm $PSFILE

    IFRAME=`echo $IFRAME | awk '{printf("%05d\n"), $1+1}'`
done

ffmpeg -y -framerate 7 -i seis_movie_%05d.jpg -c:v libx264 -r 7 -pix_fmt yuv420p out.mp4


#####
#	CLEAN UP
#####
rm *.tmp
rm mech.cpt dep.cpt
rm gmt.*

