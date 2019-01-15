#!/bin/bash

if [ $# -ne 2 ]
then
    echo Need to specify a FFM and threshold fraction
    exit
fi
FFM="$1"
THR="$2"

MAX=`awk 'BEGIN{mx=0}{if(NF==11&&substr($1,1,1)!="#"&&$4>mx){mx=$4}}END{print mx}' $FFM`
awk '{
    if (NF==11 && substr($1,1,1)!="#" && $4<'"$THR"'*'"$MAX"') {
        print $1,$2,$3,0,$5,$6,$7,$8,$9,$10,$11
    } else {
        print $0
    }
}' $FFM
