#!/bin/bash

#####
#	SET PATH TO HDEF EXECUTABLE
#####
# Check if distaz2lola is set in PATH
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which distaz2lola | xargs dirname)
fi

# Check for distaz2lola in same directory as script
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/distaz2lola | xargs dirname)
fi

# Check for distaz2lola in relative directory ../bin (assumes script is in Hdef/dir)
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/../bin/distaz2lola | xargs dirname)
fi

# Check for distaz2lola in relative directory ../build (assumes script is in Hdef/dir)
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/../build/distaz2lola | xargs dirname)
fi

# Hdef executables are required!
if [ "$BIN_DIR" == "" ]
then
    echo "$0: unable to find Hdef executable distaz2lola; exiting" 1>&2
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


echo 0.70711135 0.70705750 > answer.tmp

# Standard input
echo 0 0 111.19 45 | $BIN_DIR/distaz2lola -s > lola.tmp
$TEST_BIN_DIR/test_values.sh lola.tmp answer.tmp 2 "distaz2lola: stdin (-s)" || exit 1

echo 0 0 111.19 45 | $BIN_DIR/distaz2lola -f stdin > lola.tmp
$TEST_BIN_DIR/test_values.sh lola.tmp answer.tmp 2 "distaz2lola: stdin (-f stdin)" || exit 1

# File
echo 0 0 111.19 45 > distaz.tmp
$BIN_DIR/distaz2lola -f distaz.tmp -o lola.tmp
$TEST_BIN_DIR/test_values.sh lola.tmp answer.tmp 2 "distaz2lola: file input" || exit 1

# Command line input
$BIN_DIR/distaz2lola -c 0 0 111.19 45 -o lola.tmp
$TEST_BIN_DIR/test_values.sh lola.tmp answer.tmp 2 "distaz2lola: command line input" || exit 1

# Non-zero starting point
echo 237.76606699703069 47.001414599041539 > answer.tmp
$BIN_DIR/distaz2lola -c 5 52 7917 325 > lola.tmp
$TEST_BIN_DIR/test_values.sh lola.tmp answer.tmp 2 "distaz2lola: command line input, non-zero start" || exit 1

# Two answers
echo 0.70711135 0.70705750 > answer.tmp
echo 237.76606699703069 47.001414599041539 >> answer.tmp
echo 0 0 111.19 45 > distaz.tmp
echo 5 52 7917 325 >> distaz.tmp
$BIN_DIR/distaz2lola -f distaz.tmp > lola.tmp
$TEST_BIN_DIR/test_values.sh lola.tmp answer.tmp 2 "distaz2lola: file input, two entries" || exit 1
