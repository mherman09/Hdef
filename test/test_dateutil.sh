#!/bin/bash

#####
#	SET PATH TO TEST_VALUES SCRIPT
#####
TEST_BIN_DIR=$(echo $0 | xargs dirname)


#####
#	SET PATH TO HDEF EXECUTABLE
#####
# Check for o92util
$TEST_BIN_DIR/test_find_hdef_exec.sh dateutil || { echo "$0: could not find dateutil; exiting" 1>&2; exit 1; }
BIN_DIR=$(cat hdefexec.tmp | xargs dirname)


#####
#	RUN TEST
#####
trap "rm -f *.tmp" 0 1 2 3 8 9


# Compute number of days from default format date, read from file
echo 1.0000231480225921 > answer.tmp
echo 2011-01-01T00:00:00 2011-01-02T00:00:02 > date.tmp
$BIN_DIR/dateutil -f date.tmp -nday > nday.tmp
$TEST_BIN_DIR/test_values.sh nday.tmp answer.tmp 1 "dateutil: default format (YYYY-MM-DDTHH:MM:SS), from file" || exit 1

# Compute number of days from long, space-separated date format, read from stdin
echo 2011 01 01 00 00 00 000 2011 01 02 00 00 02 000 |\
    $BIN_DIR/dateutil -nday -format "YYYY MM DD HH MM SS MMM" > nday.tmp
$TEST_BIN_DIR/test_values.sh nday.tmp answer.tmp 1 "dateutil: long format with spaces, stdin" || exit 1

# Compute date from starting date, number of days, with short, no-space date format, read from file
echo 20110101 13 > nday.tmp
$BIN_DIR/dateutil -f nday.tmp -date -format "YYYYMMDD" -o date.tmp
TEST=`awk '{if($1!="2011-01-14"){print "Failed"}else{print "Passed"}}' date.tmp`
if [ "$TEST" == "Failed" ]
then
    echo Failed test 1>&2
    exit 1
fi

# Compute two number of days with short, condensed date format, read from stdin
echo 1578 > answer.tmp
echo 1579 >> answer.tmp
echo 20190515 20230909 > date.tmp
echo 20190515 20230910 >> date.tmp
cat date.tmp | $BIN_DIR/dateutil -format "YYYYMMDD" -nday -o nday.tmp
$TEST_BIN_DIR/test_values.sh nday.tmp answer.tmp 1 "dateutil: short, condensed format, stdin, two entries" || exit 1
