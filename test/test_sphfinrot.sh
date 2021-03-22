#!/bin/bash

#####
#	SET PATH TO TEST_VALUES SCRIPT
#####
TEST_BIN_DIR=$(echo $0 | xargs dirname)


#####
#	SET PATH TO HDEF EXECUTABLE
#####
# Check for o92util
$TEST_BIN_DIR/test_find_hdef_exec.sh sphfinrot || { echo "$0: could not find sphfinrot; exiting" 1>&2; exit 1; }
BIN_DIR=$(cat hdefexec.tmp | xargs dirname)


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
