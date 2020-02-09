#!/bin/bash

BIN_DIR=`./define_bin_dir.sh`

cat > evt.tmp << EOF
340 55  74 5.7 31
 23 50 -79 6.7 61
  2 37  93 5.0 31
354 73 -87 5.6 51
  7 70  99 5.2 36
334 66  66 5.3 12
EOF
${BIN_DIR}/mtutil -sdr evt.tmp -ternary ternary.tmp
paste ternary.tmp evt.tmp | awk '{print $1,$2,$3,$7,$8}' > j
mv j evt.tmp

PSFILE="test_ternary.ps"
echo 0 0 | gmt psxy -JX1i -R0/1/0/1 -K -P --PS_MEDIA=8.5ix11i > $PSFILE


${BIN_DIR}/ternary.sh -f evt.tmp -a $PSFILE

gmt makecpt -T0/50/1 -Cplasma -D > dep.cpt
${BIN_DIR}/ternary.sh -f evt.tmp -a $PSFILE -c dep.cpt -x 5.0,4.0,2.0 --frac_simple 1
${BIN_DIR}/ternary.sh -f evt.tmp -a $PSFILE -c dep.cpt -x 0.0,3.0,1.5 -s 0.1 -j 0.1 --frac_simple 0


echo 0 0 | gmt psxy -JX1i -R0/1/0/1 -O >> $PSFILE
ps2pdf $PSFILE

echo "test_ternary.sh: check plots in test_ternary.pdf"

rm $PSFILE
rm gmt.*
rm *.cpt
rm *.tmp
