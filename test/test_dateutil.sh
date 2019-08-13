#!/bin/bash

# Compute number of days from default format date, read from file
echo 1.0000231480225921 > answer.tmp
echo 2011-01-01T00:00:00 2011-01-02T00:00:02 > date.tmp
../bin/dateutil -f date.tmp -nday > nday.tmp
./test_values.sh nday.tmp answer.tmp 1 "dateutil: default format (YYYY-MM-DDTHH:MM:SS), from file" || exit 1

# Compute number of days from long, space-separated date format, read from stdin
echo 2011 01 01 00 00 00 000 2011 01 02 00 00 02 000 |\
    ../bin/dateutil -nday -format "YYYY MM DD HH MM SS MMM" > nday.tmp
./test_values.sh nday.tmp answer.tmp 1 "dateutil: long format with spaces, stdin" || exit 1

# Compute date from starting date, number of days, with short, no-space date format, read from file
echo 20110101 13 > nday.tmp
../bin/dateutil -f nday.tmp -date -format "YYYYMMDD" -o date.tmp
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
cat date.tmp | ../bin/dateutil -format "YYYYMMDD" -nday -o nday.tmp
./test_values.sh nday.tmp answer.tmp 1 "dateutil: short, condensed format, stdin, two entries" || exit 1

rm *.tmp
