#!/bin/bash

BIN_DIR=`./define_bin_dir.sh`

PSFILE="test_trg_schem.ps"

echo 0 0 | gmt psxy -JX1i -R0/1/0/1 -K > $PSFILE
$BIN_DIR/trg_schem.sh -i 20,45,-100 -a $PSFILE
$BIN_DIR/trg_schem.sh -i 103,85,148 -a $PSFILE -x 4,3,2
echo 0 0 | gmt psxy -JX1i -R0/1/0/1 -O >> $PSFILE

ps2pdf $PSFILE

echo "test_trg_schem.sh: check output figure in test_trg_schem.pdf"

rm $PSFILE
rm gmt.*
