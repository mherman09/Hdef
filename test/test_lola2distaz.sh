#!/bin/bash

#####
#	SET PATH TO TEST_VALUES SCRIPT
#####
TEST_BIN_DIR=$(echo $0 | xargs dirname)


#####
#	SET PATH TO HDEF EXECUTABLE
#####
# Check for o92util
$TEST_BIN_DIR/test_find_hdef_exec.sh lola2distaz || { echo "$0: could not find lola2distaz; exiting" 1>&2; exit 1; }
BIN_DIR=$(cat hdefexec.tmp | xargs dirname)


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
