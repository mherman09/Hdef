#!/bin/bash

PSFILE="test_trg_schem.ps"

echo 0 0 | gmt psxy -JX1i -R0/1/0/1 -K > $PSFILE
../bin/trg_schem.sh -i 20,45,-100 -a $PSFILE
../bin/trg_schem.sh -i 103,85,148 -a $PSFILE -x 4,3,2
echo 0 0 | gmt psxy -JX1i -R0/1/0/1 -O >> $PSFILE

ps2pdf $PSFILE
rm $PSFILE
rm gmt.*


