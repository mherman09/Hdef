#!/bin/bash

#####
#	SET PATH TO HDEF EXECUTABLE
#####
# Check if sphfinrot is set in PATH
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which sphfinrot | xargs dirname)
fi

# Check for sphfinrot in same directory as script
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/sphfinrot | xargs dirname)
fi

# Check for sphfinrot in relative directory ../bin (assumes script is in Hdef/dir)
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/../bin/sphfinrot | xargs dirname)
fi

# Check for sphfinrot in relative directory ../build (assumes script is in Hdef/dir)
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/../build/sphfinrot | xargs dirname)
fi

# Hdef executables are required!
if [ "$BIN_DIR" == "" ]
then
    echo "$0: unable to find Hdef executable sphfinrot; exiting" 1>&2
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


echo 0 0 | $BIN_DIR/sphfinrot -pole 90 0 -angle 45 > coord.tmp
echo 0 -45 > answer.tmp
$TEST_BIN_DIR/test_values.sh coord.tmp answer.tmp 2 "sphfinrot: stdin, stdout" || exit 1

cat > lonlat.tmp << EOF
0 0
0 90
EOF
$BIN_DIR/sphfinrot -f lonlat.tmp -o newlonlat.tmp -pole 135 0 -angle 45
cat > answer.tmp << EOF
 -9.735610E+00 -3.000000E+01
  4.500000E+01  4.500000E+01
EOF
$TEST_BIN_DIR/test_values.sh newlonlat.tmp answer.tmp 2 "sphfinrot: file in, file out" || exit 1
