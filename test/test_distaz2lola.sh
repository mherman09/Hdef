#!/bin/bash

echo 0.70711135 0.70705750 > answer.tmp

# Standard input
echo 0 0 111.19 45 | ../bin/distaz2lola -s > lola.tmp
./test_values.sh lola.tmp answer.tmp 2 "distaz2lola: stdin (-s)" || exit 1

echo 0 0 111.19 45 | ../bin/distaz2lola -f stdin > lola.tmp
./test_values.sh lola.tmp answer.tmp 2 "distaz2lola: stdin (-f stdin)" || exit 1

# File
echo 0 0 111.19 45 > distaz.tmp
../bin/distaz2lola -f distaz.tmp -o lola.tmp
./test_values.sh lola.tmp answer.tmp 2 "distaz2lola: file input" || exit 1

# Command line input
../bin/distaz2lola -c 0 0 111.19 45 -o lola.tmp
./test_values.sh lola.tmp answer.tmp 2 "distaz2lola: command line input" || exit 1

# Non-zero starting point
echo 237.76606699703069 47.001414599041539 > answer.tmp
../bin/distaz2lola -c 5 52 7917 325 > lola.tmp
./test_values.sh lola.tmp answer.tmp 2 "distaz2lola: command line input, non-zero start" || exit 1

# Two answers
echo 0.70711135 0.70705750 > answer.tmp
echo 237.76606699703069 47.001414599041539 >> answer.tmp
echo 0 0 111.19 45 > distaz.tmp
echo 5 52 7917 325 >> distaz.tmp
../bin/distaz2lola -f distaz.tmp > lola.tmp
./test_values.sh lola.tmp answer.tmp 2 "distaz2lola: file input, two entries" || exit 1


rm *.tmp
