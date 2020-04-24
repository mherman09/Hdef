#!/bin/bash

#####
#	SET PATH TO HDEF EXECUTABLE
#####
# Check if rangen is set in PATH
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which rangen | xargs dirname)
fi

# Check for rangen in same directory as script
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/rangen | xargs dirname)
fi

# Check for rangen in relative directory ../bin (assumes script is in Hdef/dir)
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/../bin/rangen | xargs dirname)
fi

# Check for rangen in relative directory ../build (assumes script is in Hdef/dir)
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/../build/rangen | xargs dirname)
fi

# Hdef executables are required!
if [ "$BIN_DIR" == "" ]
then
    echo "$0: unable to find Hdef executable rangen; exiting" 1>&2
    exit 1
fi


#####
#	SET PATH TO TEST_VALUES SCRIPT
#####
TEST_BIN_DIR=`echo $0 | xargs dirname`


#####
#	RUN TEST
#####
trap "rm -f *.tmp" 0 1 2 3 8 9


N=10
$BIN_DIR/rangen -uniform 2 3 -npts $N -seed 82766 > uniform.tmp
cat > answer.tmp << EOF
   2.6477572780914285
   2.8565734572519910
   2.4300967940973788
   2.6368187761700863
   2.0131716555970156
   2.3760156307059308
   2.6947056081220313
   2.9171563232379003
   2.6463254730091719
   2.7922254384972267
EOF
$TEST_BIN_DIR/test_values.sh uniform.tmp answer.tmp 1 "rangen: uniform distribution" || exit 1


$BIN_DIR/rangen -normal 8.6 1.5 -npts $N -seed 2742 > normal.tmp
cat > answer.tmp << EOF
   6.7501018934908483
   8.0178941774922432
   7.9494583010276116
   8.4785506395132035
   9.8949054230567022
   8.7041088310506201
   9.1059686525888495
   7.1789407725551984
   10.429223977511743
   8.1593962443817283
EOF
$TEST_BIN_DIR/test_values.sh normal.tmp answer.tmp 1 "rangen: normal distribution" || exit 1




# gmt set PS_MEDIA 8.5ix11i
#
# PSFILE="test_rangen.ps"
#
# N=10000
# BIN_WID=0.1
# HEIGHT=2000
# YTIKS=200
# $BIN_DIR/rangen -uniform 2 3 -npts $N -seed 82766 > uniform.tmp
# gmt pshistogram uniform.tmp -JX4i -R2/3/0/$HEIGHT -W$BIN_WID -G205/205/225 -L0.5p -Bxa0.1 -Bya$YTIKS -P -K > $PSFILE
#
# BIN_WID=0.05
# HEIGHT=300
# YTIKS=100
# $BIN_DIR/rangen -normal 8.6 1.5 -npts $N -seed 2742 > normal.tmp
# gmt pshistogram normal.tmp -JX6i/4i -R4/14/0/$HEIGHT -W$BIN_WID -G205/205/225 -L0.5p -Bxa1 -Bya$YTIKS -P -Y5i -O >> $PSFILE
#
# ps2pdf $PSFILE
#
# #####
# #	CLEAN UP
# #####
# rm $PSFILE
# rm gmt.*
# rm *.tmp
