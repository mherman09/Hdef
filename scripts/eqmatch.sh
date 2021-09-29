#!/bin/bash

function usage() {
    echo "Usage: $0 INPUT_FILE_1 INPUT_FILE_2 [...options...]" 1>&2
    echo "" 1>&2
    echo "INPUT_FILE_1           YYYY-MM-DD HH-MM-SS ..." 1>&2
    echo "INPUT_FILE_2           YYYY-MM-DD HH-MM-SS ..." 1>&2
    echo "--threshold:sec SEC    Matching threshold (default: 5 seconds)" 1>&2
    exit 1
}

if [ $# -lt 2 ]
then
    usage
fi

INPUT_FILE_1=$1
INPUT_FILE_2=$2

if [ ! -f $INPUT_FILE_1 ]
then
    echo "$0: no file found named \"$INPUT_FILE_1\"" >&2
    usage
fi
if [ ! -f $INPUT_FILE_2 ]
then
    echo "$0: no file found named \"$INPUT_FILE_2\"" >&2
    usage
fi


MATCH_THRESHOLD=5


#####
#
#####
DATE_REF=$(tail -1 $INPUT_FILE_1 | awk '{print $1 "T" $2}')

awk '{print "'$DATE_REF'",$1 "T" $2}' $INPUT_FILE_1 |\
    dateutil -nday -format "YYYY-MM-DDTHH:MM:SS" |\
    awk '{print NR,$1}' |\
    sort -gk2 > nday_1.tmp
awk '{print "'$DATE_REF'",$1 "T" $2}' $INPUT_FILE_2 |\
    dateutil -nday -format "YYYY-MM-DDTHH:MM:SS" |\
    awk '{print NR,$1}' |\
     sort -gk2 > nday_2.tmp

N1=$(wc nday_1.tmp | awk '{print $1}')
N2=$(wc nday_2.tmp | awk '{print $1}')

cat nday_1.tmp nday_2.tmp |\
    awk '{
        if (NR<='$N1') {
            nday_1[$1] = $2
        } else {
            nday_2[$1] = $2
        }
     } END {
        thr = '$MATCH_THRESHOLD'/24/60/60
        i = 1
        j = 1
        while (i<='$N1' || j<='$N2') {
            dday = nday_1[i] - nday_2[j]
            if (dday<-thr) {
                #print "no match: nday_1 is too much smaller than nday_2" >> "/dev/stderr"
                i++
            } else if (dday>thr) {
                #print "no match: nday_1 is too much larger than nday_2" >> "/dev/stderr"
                j++
            } else {
                if (i<='$N1' && j<='$N2') {
                    #print "MATCH",i,j,nday_1[i],nday_2[j],'$N2' >> "/dev/stderr"
                    print i,j
                }
                i++
                j++
            }
        }
    }' > index_list.tmp

if [ -f match_1.tmp ]; then rm match_1.tmp; fi
if [ -f match_2.tmp ]; then rm match_2.tmp; fi
while read I J
do
    sed -n -e "${I}p" $INPUT_FILE_1 >> match_1.tmp
    sed -n -e "${J}p" $INPUT_FILE_2 >> match_2.tmp
done < index_list.tmp

paste match_1.tmp match_2.tmp

