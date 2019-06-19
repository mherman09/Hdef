#!/bin/bash

function usage {
    echo "Usage: $0 COMPUTED EXPECTED NFIELDS MESSAGE" 1>&2
    exit 1
}

# Check for files and number of fields
if [ ! -f "$1" ]; then echo "Missing computed value file" 1>&2; usage; fi
if [ ! -f "$2" ]; then echo "Missing expected value file" 1>&2; usage; fi
if [ $# -lt 3 ]; then echo "Missing number of fields" 1>&2; usage; fi

# Compare computed and expected values
paste $1 $2 |\
    awk '{
        for (i=1;i<='$3';i++) {
            computed = $i
            expected = $('$3'+i)
            diff = computed - expected
            if (diff<0) {
                diff = -diff
            }
            if (expected<0) {
                expected = -expected
            }
            if (expected<1e-8) {
                if (diff>1.0e-7) {
                    print "FAIL"
                    printf("computed %14.6e\n"),computed > "/dev/stderr"
                    printf("expected %14.6e\n"),expected > "/dev/stderr"
                    exit
                } else {
                    print "PASS"
                }
            } else {
                if (diff/expected>1.0e-7) {
                    print "FAIL"
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
    echo Failed test $4 1>&2
    exit 1
fi
