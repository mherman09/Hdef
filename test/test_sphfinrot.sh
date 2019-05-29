#!/bin/bash

echo 0 0 | ../bin/sphfinrot -pole 90 0 -angle 45 > coord.tmp
echo 0 -45 > answer.tmp
test_values.sh coord.tmp answer.tmp 2 "sphfinrot: stdin, stdout" || exit 1

cat > lonlat.tmp << EOF
0 0
0 90
EOF
../bin/sphfinrot -f lonlat.tmp -o newlonlat.tmp -pole 135 0 -angle 45
cat > answer.tmp << EOF
 -9.735610E+00 -3.000000E+01
  4.500000E+01  4.500000E+01
EOF
test_values.sh newlonlat.tmp answer.tmp 2 "sphfinrot: file in, file out" || exit 1

rm *.tmp



