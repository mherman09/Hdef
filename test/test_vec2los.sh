#!/bin/bash


echo 4 5 6 1 2 3 > vec.tmp
cat vec.tmp | ../bin/vec2los -a 45 -i 30 > los.tmp
cat > answer.tmp << EOF
4 5 6 0.33711730708738386
EOF
test_values.sh los.tmp answer.tmp 4 "vec2los: standard in/out"

echo 4 5 6 -1 -2 -3 >> vec.tmp
echo 45 30 > look.tmp
echo 25 20 >> look.tmp
cat > answer.tmp << EOF
4 5 6 0.33711730708738386
4 5 6 -1.0743723112683896
EOF
../bin/vec2los -f vec.tmp -look look.tmp > los.tmp
test_values.sh los.tmp answer.tmp 4 "vec2los: read look orientations"

rm *.tmp
