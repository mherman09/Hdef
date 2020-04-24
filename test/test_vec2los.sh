#!/bin/bash

#####
#	SET PATH TO HDEF EXECUTABLE
#####
# Check if vec2los is set in PATH
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which vec2los | xargs dirname)
fi

# Check for vec2los in same directory as script
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/vec2los | xargs dirname)
fi

# Check for vec2los in relative directory ../bin (assumes script is in Hdef/dir)
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/../bin/vec2los | xargs dirname)
fi

# Check for vec2los in relative directory ../build (assumes script is in Hdef/dir)
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/../build/vec2los | xargs dirname)
fi

# Hdef executables are required!
if [ "$BIN_DIR" == "" ]
then
    echo "$0: unable to find Hdef executable vec2los; exiting" 1>&2
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


echo 4 5 6 1 2 3 > vec.tmp
cat vec.tmp | $BIN_DIR/vec2los -a 45 -i 30 > los.tmp
cat > answer.tmp << EOF
4 5 6 0.33711730708738386
EOF
$TEST_BIN_DIR/test_values.sh los.tmp answer.tmp 4 "vec2los: standard in/out"

echo 4 5 6 -1 -2 -3 >> vec.tmp
echo 45 30 > look.tmp
echo 25 20 >> look.tmp
cat > answer.tmp << EOF
4 5 6 0.33711730708738386
4 5 6 -1.0743723112683896
EOF
$BIN_DIR/vec2los -f vec.tmp -look look.tmp > los.tmp
$TEST_BIN_DIR/test_values.sh los.tmp answer.tmp 4 "vec2los: read look orientations"
