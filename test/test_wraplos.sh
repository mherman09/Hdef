#!/bin/bash

#####
#	SET PATH TO TEST_VALUES SCRIPT
#####
TEST_BIN_DIR=$(echo $0 | xargs dirname)


#####
#	SET PATH TO HDEF EXECUTABLE
#####
# Check for o92util
$TEST_BIN_DIR/test_find_hdef_exec.sh wraplos || { echo "$0: could not find wraplos; exiting" 1>&2; exit 1; }
BIN_DIR=$(cat hdefexec.tmp | xargs dirname)


#####
#	RUN TEST
#####
# Clean up tmp files
trap "rm -f *.tmp" 0 1 2 3 8 9

# Standard input to wraplos
echo 1 2 3 4 > los.tmp
cat los.tmp | $BIN_DIR/wraplos -w 0.17 > phase.tmp
cat > answer.tmp << EOF
   1.0000000000000000        2.0000000000000000        3.0000000000000000       0.36959913571640612
EOF
$TEST_BIN_DIR/test_values.sh phase.tmp answer.tmp 4 "wraplos: read input from stdin" || \
    { echo "$0: error in wraplos reading from standard input" 1>&2; exit 1; }

# Read wraplos input from file
echo 5 6 7 8 >> los.tmp
cat > answer.tmp << EOF
   1.0000000000000000        2.0000000000000000        3.0000000000000000       0.36959913571640612
   5.0000000000000000        6.0000000000000000        7.0000000000000000       0.73919827143281225
EOF
$BIN_DIR/wraplos -f los.tmp -w 0.17 > phase.tmp
$TEST_BIN_DIR/test_values.sh phase.tmp answer.tmp 4 "wraplos: read input from file" || \
    { echo "$0: error in wraplos reading input from file" 1>&2; exit 1; }
