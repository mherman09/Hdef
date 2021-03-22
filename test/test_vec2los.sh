#!/bin/bash

#####
#	SET PATH TO TEST_VALUES SCRIPT
#####
TEST_BIN_DIR=$(echo $0 | xargs dirname)


#####
#	SET PATH TO HDEF EXECUTABLE
#####
# Check for o92util
$TEST_BIN_DIR/test_find_hdef_exec.sh vec2los || { echo "$0: could not find vec2los; exiting" 1>&2; exit 1; }
BIN_DIR=$(cat hdefexec.tmp | xargs dirname)


#####
#       RUN TEST
#####
# Clean up tmp files
trap "rm -f *.tmp" 0 1 2 3 8 9

# Standard input to vec2los
echo 4 5 6 1 2 3 > vec.tmp
cat vec.tmp | $BIN_DIR/vec2los -a 45 -i 30 > los.tmp
cat > answer.tmp << EOF
4 5 6 0.33711730708738386
EOF
$TEST_BIN_DIR/test_values.sh los.tmp answer.tmp 4 "vec2los: read input from stdin"  || \
    { echo "$0: error in vec2los reading from standard input" 1>&2; exit 1; }

# Read vec2los input from file
echo 4 5 6 -1 -2 -3 >> vec.tmp
echo 45 30 > look.tmp
echo 25 20 >> look.tmp
cat > answer.tmp << EOF
4 5 6 0.33711730708738386
4 5 6 -1.0743723112683896
EOF
$BIN_DIR/vec2los -f vec.tmp -look look.tmp > los.tmp
$TEST_BIN_DIR/test_values.sh los.tmp answer.tmp 4 "vec2los: read input from file" || \
    { echo "$0: error in vec2los reading input from file" 1>&2; exit 1; }
