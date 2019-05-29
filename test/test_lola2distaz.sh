#!/bin/bash

echo 324.14147153 134.36681298 > answer.tmp
echo 35 -44 38 -46 | ../bin/lola2distaz -s > lola.tmp
test_values.sh lola.tmp answer.tmp 2 "lola2distaz: stdin"
echo 35 -44 38 -46 > distaz.tmp
../bin/lola2distaz -f distaz.tmp -o lola.tmp
test_values.sh lola.tmp answer.tmp 2 "lola2distaz: file input"
../bin/lola2distaz -c 35 -44 38 -46 -o lola.tmp
test_values.sh lola.tmp answer.tmp 2 "lola2distaz: command line input"

rm *.tmp
