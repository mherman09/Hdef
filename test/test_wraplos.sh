#!/bin/bash

#####
#	SET PATH TO HDEF EXECUTABLE
#####
# Check if wraplos is set in PATH
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which wraplos | xargs dirname)
fi

# Check for wraplos in same directory as script
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/wraplos | xargs dirname)
fi

# Check for wraplos in relative directory ../bin (assumes script is in Hdef/dir)
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/../bin/wraplos | xargs dirname)
fi

# Check for wraplos in relative directory ../build (assumes script is in Hdef/dir)
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/../build/wraplos | xargs dirname)
fi

# Hdef executables are required!
if [ "$BIN_DIR" == "" ]
then
    echo "$0: unable to find Hdef executable wraplos; exiting" 1>&2
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


echo 1 2 3 4 > los.tmp
cat los.tmp | $BIN_DIR/wraplos -w 0.17 > phase.tmp
cat > answer.tmp << EOF
   1.0000000000000000        2.0000000000000000        3.0000000000000000       0.36959913571640612
EOF
$TEST_BIN_DIR/test_values.sh phase.tmp answer.tmp 4 "wraplos: standard in/out"

echo 5 6 7 8 >> los.tmp
cat > answer.tmp << EOF
   1.0000000000000000        2.0000000000000000        3.0000000000000000       0.36959913571640612
   5.0000000000000000        6.0000000000000000        7.0000000000000000       0.73919827143281225
EOF
$BIN_DIR/wraplos -f los.tmp -w 0.17 > phase.tmp
$TEST_BIN_DIR/test_values.sh phase.tmp answer.tmp 4 "wraplos: read input"
