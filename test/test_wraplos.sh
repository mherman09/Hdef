#!/bin/bash


echo 1 2 3 4 > los.tmp
cat los.tmp | ../bin/wraplos -w 0.17 > phase.tmp
cat > answer.tmp << EOF
   1.0000000000000000        2.0000000000000000        3.0000000000000000       0.36959913571640612
EOF
./test_values.sh phase.tmp answer.tmp 4 "wraplos: standard in/out"

echo 5 6 7 8 >> los.tmp
cat > answer.tmp << EOF
   1.0000000000000000        2.0000000000000000        3.0000000000000000       0.36959913571640612
   5.0000000000000000        6.0000000000000000        7.0000000000000000       0.73919827143281225     
EOF
../bin/wraplos -f los.tmp -w 0.17 > phase.tmp
./test_values.sh phase.tmp answer.tmp 4 "wraplos: read input"

rm *.tmp
