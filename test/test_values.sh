#!/bin/bash

function usage {
    echo "Usage: $0 COMPUTED EXPECTED NFIELDS MESSAGE [-zero THR]" 1>&2
    exit 1
}

# Check for files and number of fields
COMPUTED="$1"
EXPECTED="$2"
NFIELDS="$3"
MESSAGE="$4"
if [ ! -f "$COMPUTED" ]; then echo "Missing computed value file" 1>&2; usage; fi
if [ ! -f "$EXPECTED" ]; then echo "Missing expected value file" 1>&2; usage; fi
if [ "$NFIELDS" == "" ]; then echo "Missing number of fields" 1>&2; usage; fi
if [ "$MESSAGE" == "" ]; then echo "Missing descriptive error message" 1>&2; usage; fi

# Parse optional arguments
THRESHOLD="1.0e-6"
if [ $# -ge 5 ]
then
    shift;shift;shift;shift
    while [ "$1" != "" ]
    do
        case $1 in
            -zero) shift;THRESHOLD=$1;echo "$0: new zero level set to $THRESHOLD";;
            *) echo "$0: no option \"$1\"" 1>&2;;
        esac
        shift
    done
fi

# Trap result.tmp file
trap "rm -f result.tmp" 0 1 2 3 8 9

# Compare computed and expected values
paste $COMPUTED $EXPECTED |\
    awk '{
        for (i=1;i<='$NFIELDS';i++) {
            computed = $i
            expected = $('$NFIELDS'+i)
            diff = computed - expected

            # Absolute value of difference
            if (diff<0) {
                diff = -diff
            }

            # Absolute value of expected (defined by user)
            if (expected<0) {
                expected_abs = -expected
            } else {
                expected_abs = expected
            }

            # Compare calculated difference with allowed threshold difference
            if (expected_abs<'"$THRESHOLD"'/10) {
                if (diff>'"$THRESHOLD"') {
                    print "FAIL",NR
                    printf("computed %14.6e\n"),computed > "/dev/stderr"
                    printf("expected %14.6e\n"),expected > "/dev/stderr"
                    exit
                } else {
                    print "PASS"
                }
            } else {
                if (diff/expected_abs>'"$THRESHOLD"') {
                    print "FAIL",NR
                    printf("computed %16.8e\n"),computed > "/dev/stderr"
                    printf("expected %16.8e\n"),expected > "/dev/stderr"
                    exit
                } else {
                    print "PASS"
                }
            }
        }
    }' > result.tmp
FAIL=`grep FAIL result.tmp`
if [ "$FAIL" != "" ]
then
    echo "Failed test $MESSAGE" 1>&2
    echo
    OFFENDING_LINE=$(echo $FAIL | awk '{print $2}')
    echo "Offending line of \"$COMPUTED\" is line $OFFENDING_LINE:" 1>&2
    sed -ne "${OFFENDING_LINE}p" $COMPUTED 1>&2
    echo
    echo "Expected to see:" 1>&2
    sed -ne "${OFFENDING_LINE}p" $EXPECTED 1>&2
    exit 1
fi
