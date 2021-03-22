#!/bin/bash

#####
#	SET PATH TO TEST_VALUES SCRIPT
#####
TEST_BIN_DIR=$(echo $0 | xargs dirname)


#####
#	SET PATH TO HDEF EXECUTABLE
#####
# Check for o92util
$TEST_BIN_DIR/test_find_hdef_exec.sh distaz2lola || { echo "$0: could not find distaz2lola; exiting" 1>&2; exit 1; }
BIN_DIR=$(cat hdefexec.tmp | xargs dirname)


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
