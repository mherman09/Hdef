#!/bin/bash

set -e

#####
#	Script: eqmatch.sh
#
#	Match two earthquake catalogs by event origin time, plus other user-specified parameters
#####
# Usage statement
function usage() {
    echo "Usage: $0 INPUT_FILE_1 INPUT_FILE_2 [...options...]" 1>&2
    echo "" 1>&2
    echo "INPUT_FILE_1           YYYY-MM-DD  HH-MM-SS  ..." 1>&2
    echo "INPUT_FILE_2           YYYY-MM-DD  HH-MM-SS  ..." 1>&2
    echo "--threshold:sec SEC    Matching threshold (default: 5 seconds)" 1>&2
    echo "--field FIELD THR      Match by values in FIELD with threshold THR" 1>&2
    exit 1
}


# Two file names are required
if [ $# -lt 2 ]
then
    usage
fi

INPUT_FILE_1=$1
INPUT_FILE_2=$2
shift; shift

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


# Parse optional command line arguments
MATCH_THRESHOLD=5             # Threshold difference for matching events (in seconds)
EXTRA_FIELDS=                 # Extra parameters to define matching events (column)
EXTRA_THRESHOLDS=             # Extra parameters to define matching events (threshold difference)
while [ "$1" != "" ]
do
    case $1 in
        --threshold:sec) shift; MATCH_THRESHOLD=$1;;
        --field) shift
                 EXTRA_FIELDS="$EXTRA_FIELDS $1"
                 shift
                 EXTRA_THRESHOLDS="$EXTRA_THRESHOLDS $1";;
        *) echo "$0: no option $1" 1>&2; usage;;
    esac
    shift
done

# Check that number of extra fields and threshold values are the same
N_EXTRA_FIELDS=$(echo $EXTRA_FIELDS | awk '{print NF}')
N_EXTRA_THRESHOLDS=$(echo $EXTRA_THRESHOLDS | awk '{print NF}')
# echo N_EXTRA_FIELDS=$N_EXTRA_FIELDS
# echo N_EXTRA_THRESHOLDS=$N_EXTRA_THRESHOLDS
if [ $N_EXTRA_FIELDS -ne $N_EXTRA_THRESHOLDS ]
then
    echo "$0: number of extra fields does not equal the number of extra threshold values" 1>&2
    echo "EXTRA_FIELDS= $EXTRA_FIELDS" 1>&2
    echo "EXTRA_THRESHOLDS= $EXTRA_THRESHOLDS" 1>&2
    exit 1
fi


#####
#	CHECK FILES FOR COMMON ERRORS
#####
# Are there empty lines?
N_EMPTY_LINES=$(sed -ne "/^$/p" $INPUT_FILE_1 | wc | awk '{print $1}')
if [ $N_EMPTY_LINES -ne 0 ]
then
    echo "$0: input file 1 has empty lines" 1>&2; exit 1
fi
N_EMPTY_LINES=$(sed -ne "/^$/p" $INPUT_FILE_2 | wc | awk '{print $1}')
if [ $N_EMPTY_LINES -ne 0 ]
then
    echo "$0: input file 2 has empty lines" 1>&2; exit 1
fi


#####
#	SET UP FILES TO MATCH RECORDS
#####
# Set temporary file names
PREFIX=$(date +%s)
TMP_1=${PREFIX}_1.tmp
TMP_2=${PREFIX}_2.tmp
TMP_NDAY_1=${PREFIX}_nday_1.tmp
TMP_NDAY_2=${PREFIX}_nday_2.tmp
TMP_EXTRA_THRESHOLDS=${PREFIX}_extra_thresholds.tmp
TMP_EXTRA_1=${PREFIX}_extra_fields_1.tmp
TMP_EXTRA_2=${PREFIX}_extra_fields_2.tmp
TMP_INDEX=${PREFIX}_index_list.tmp
TMP_MATCH_1=${PREFIX}_match_1.tmp
TMP_MATCH_2=${PREFIX}_match_2.tmp


# Clean up function
function cleanup () {
    rm -f $TMP_1 \
          $TMP_2 \
          $TMP_NDAY_1 \
          $TMP_NDAY_2 \
          $TMP_EXTRA_THRESHOLDS \
          $TMP_EXTRA_1 \
          $TMP_EXTRA_2 \
          $TMP_INDEX \
          $TMP_MATCH_1 \
          $TMP_MATCH_2
}
trap "cleanup" 0 1 2 3 8 9


# Calculate number of days from the last date in the first file
DATE_REF=$(tail -1 $INPUT_FILE_1 | awk '{print $1 "T" $2}')

awk '{print "'$DATE_REF'",$1 "T" $2}' $INPUT_FILE_1 |\
    dateutil -nday -format "YYYY-MM-DDTHH:MM:SS" |\
    awk '{print NR,$1}' |\
    sort -gk2 > $TMP_NDAY_1
awk '{print "'$DATE_REF'",$1 "T" $2}' $INPUT_FILE_2 |\
    dateutil -nday -format "YYYY-MM-DDTHH:MM:SS" |\
    awk '{print NR,$1}' |\
     sort -gk2 > $TMP_NDAY_2


# Add extra fields (if specified on command line)
# Define a comma-separated list of fields to print from the input files with awk
AWK_FIELD_LIST=$(echo $EXTRA_FIELDS |\
    awk '{
        if(NF<=0){
            exit
        }
        for(i=1;i<NF;i++){
            printf("$%d,"),$i
        }
        printf("$%d\n"),$NF
    }')
# echo AWK_FIELD_LIST=$AWK_FIELD_LIST

if [ "$AWK_FIELD_LIST" != "" ]
then
    # Print the requested fields to temporary files
    awk '{print '$AWK_FIELD_LIST'}' $INPUT_FILE_1 > $TMP_EXTRA_1
    awk '{print '$AWK_FIELD_LIST'}' $INPUT_FILE_2 > $TMP_EXTRA_2

    # Check whether there were enough fields in the input files (awk leaves these blank)
    ERR=$(awk 'BEGIN{e=0}{if(NF!='$N_EXTRA_FIELDS'){e=1;exit}}END{print e}' $TMP_EXTRA_1)
    if [ "$ERR" == "1" ]
    then
        echo "$0: file 1 does not have enough input fields" 1>&2; usage
    fi
    ERR=$(awk 'BEGIN{e=0}{if(NF!='$N_EXTRA_FIELDS'){e=1;exit}}END{print e}' $TMP_EXTRA_2)
    if [ "$ERR" == "1" ]
    then
        echo "$0: file 2 does not have enough input fields" 1>&2; usage
    fi

    # Paste extra fields to nday files
    paste $TMP_NDAY_1 $TMP_EXTRA_1 > $TMP_1
    paste $TMP_NDAY_2 $TMP_EXTRA_2 > $TMP_2

    # Copy extra thresholds to a file
    echo $EXTRA_THRESHOLDS > $TMP_EXTRA_THRESHOLDS
else
    mv $TMP_NDAY_1 $TMP_1
    mv $TMP_NDAY_2 $TMP_2
    touch $TMP_EXTRA_THRESHOLDS
fi


#####
#	MATCH RECORDS
#####
# Number of events in each file
N1=$(wc $TMP_1 | awk '{print $1}')
N2=$(wc $TMP_2 | awk '{print $1}')


# Match records under given thresholds and collect file line numbers
cat $TMP_1 $TMP_2 $TMP_EXTRA |\
    awk '{
        if (NR<='$N1') {
            nday_1[$1] = $2
            for (n=1;n<='$N_EXTRA_FIELDS';n++){
                extra_1[$1,n] = $(n+2)
            }
        } else if (NR<='$N1'+'$N2') {
            nday_2[$1] = $2
            for (n=1;n<='$N_EXTRA_FIELDS';n++){
                extra_2[$1,n] = $(n+2)
            }
        } else {
            for (n=1;n<='$N_EXTRA_FIELDS';n++) {
                ethr[n] = $n
            }
        }
     } END {
        thr = '$MATCH_THRESHOLD'/60/60/24    # Origin time matching threshold in days
        i = 1
        j = 1
        while (i<='$N1' || j<='$N2') {
            dday = nday_1[i] - nday_2[j]
            for (n=1;n<='$N_EXTRA_FIELDS';n++) {
                dextra[n] = extra_1[i,n] - extra_2[j,n]
                if (dextra[n]<0) { dextra[n] = -dextra[n] }
            }
            if (dday<-thr) {
                #print "no match: nday_1 is too much smaller than nday_2" >> "/dev/stderr"
                i++
            } else if (dday>thr) {
                #print "no match: nday_1 is too much larger than nday_2" >> "/dev/stderr"
                j++
            } else {
                if (i<='$N1' && j<='$N2') {
                    #print "MATCH",i,j,nday_1[i],nday_2[j],'$N2' >> "/dev/stderr"
                    isMatch = 1
                    for (n=1;n<='$N_EXTRA_FIELDS';n++) {
                        if (dextra[n]>ethr[n]) {
                            #print "    actually, no match: extra field does not match:",extra_1[i,n],extra_2[j,n] >> "/dev/stderr"
                            isMatch = 0
                            j++
                            break
                        }
                    }
                    if (isMatch == 0) {continue}
                    print i,j
                }
                i++
                j++
            }
        }
    }' > $TMP_INDEX


# Collect matching lines and merge
if [ -f $TMP_MATCH_1 ]; then rm $TMP_MATCH_1; fi
if [ -f $TMP_MATCH_2 ]; then rm $TMP_MATCH_2; fi
while read I J
do
    sed -n -e "${I}p" $INPUT_FILE_1 >> $TMP_MATCH_1
    sed -n -e "${J}p" $INPUT_FILE_2 >> $TMP_MATCH_2
done < $TMP_INDEX

paste $TMP_MATCH_1 $TMP_MATCH_2
