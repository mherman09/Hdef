#!/bin/bash

#####
#	SET PATH TO TEST_VALUES SCRIPT
#####
TEST_BIN_DIR=$(echo $0 | xargs dirname)


#####
#	SET PATH TO HDEF EXECUTABLE
#####
# Check for o92util
$TEST_BIN_DIR/test_find_hdef_exec.sh platemotion || { echo "$0: could not find platemotion; exiting" 1>&2; exit 1; }
BIN_DIR=$(cat hdefexec.tmp | xargs dirname)


#####
#	RUN TEST
#####
trap "rm -f *.tmp" 0 1 2 3 8 9


# MORVEL56
echo -30 40 | $BIN_DIR/platemotion -plates NA/EU -model MORVEL56 -o morvel56.tmp
echo   -30.000000000000000        40.000000000000000        22.771560776120207       -2.0222476466178430 > answer.tmp
$TEST_BIN_DIR/test_values.sh morvel56.tmp answer.tmp 4 "platemotion: MORVEL56 EU wrt NA, output file" || exit 1

# MORVEL
echo -75 -35 | $BIN_DIR/platemotion -plates SA/NZ -model MORVEL > morvel.tmp
echo -75.000000000000000       -35.000000000000000        72.008593269570923        16.689689796741547 > answer.tmp
$TEST_BIN_DIR/test_values.sh morvel.tmp answer.tmp 4 "platemotion: MORVEL SA wrt NZ, stdout" || exit 1

# NUVEL1A
echo 115 -50 > loc.tmp
echo 120 -50 >> loc.tmp
$BIN_DIR/platemotion -f loc.tmp -plates AN/AU -model NUVEL1A > nuvel1a.tmp
cat > answer.tmp << EOF
   115.00000000000000       -50.000000000000000        22.881201729386667        68.378823029969183
   120.00000000000000       -50.000000000000000        18.268958516527089        69.521449139059513
EOF
$TEST_BIN_DIR/test_values.sh nuvel1a.tmp answer.tmp 4 "platemotion: NUVEL1A AU wrt AN, stdout" || exit 1

# ITRF08
$BIN_DIR/platemotion -f loc.tmp -plates AN/AU -model ITRF08 > itrf08.tmp
cat > answer.tmp << EOF
   115.00000000000000       -50.000000000000000        25.766441315239788        68.326351344076343
   120.00000000000000       -50.000000000000000        21.150858988809414        69.673756121784024
EOF
$TEST_BIN_DIR/test_values.sh itrf08.tmp answer.tmp 4 "platemotion: ITRF08 AU wrt AN, stdout" || exit 1

# Custom pole
echo -30 40 | $BIN_DIR/platemotion -pole 139.461/61.796/0.211 -o custom.tmp
echo   -30.000000000000000        40.000000000000000        22.846411219358291       -2.0281377828526641 > answer.tmp
$TEST_BIN_DIR/test_values.sh custom.tmp answer.tmp 4 "platemotion: Custom pole, file" || exit 1
