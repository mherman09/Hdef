#!/bin/bash

#####
#	SET PATH TO HDEF EXECUTABLE
#####
# Check if lola2distaz is set in PATH
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which lola2distaz | xargs dirname)
fi

# Check for lola2distaz in same directory as script
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/lola2distaz | xargs dirname)
fi

# Check for lola2distaz in relative directory ../bin (assumes script is in Hdef/dir)
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/../bin/lola2distaz | xargs dirname)
fi

# Check for lola2distaz in relative directory ../build (assumes script is in Hdef/dir)
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/../build/lola2distaz | xargs dirname)
fi

# Hdef executables are required!
if [ "$BIN_DIR" == "" ]
then
    echo "$0: unable to find Hdef executable lola2distaz; exiting" 1>&2
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


echo 324.14147153 134.36681298 > answer.tmp
echo 35 -44 38 -46 | $BIN_DIR/lola2distaz -s > lola.tmp
$TEST_BIN_DIR/test_values.sh lola.tmp answer.tmp 2 "lola2distaz: stdin"
echo 35 -44 38 -46 > distaz.tmp
$BIN_DIR/lola2distaz -f distaz.tmp -o lola.tmp
$TEST_BIN_DIR/test_values.sh lola.tmp answer.tmp 2 "lola2distaz: file input"
$BIN_DIR/lola2distaz -c 35 -44 38 -46 -o lola.tmp
$TEST_BIN_DIR/test_values.sh lola.tmp answer.tmp 2 "lola2distaz: command line input"
