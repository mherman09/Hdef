#!/bin/bash

# MORVEL56
echo -30 40 | ../bin/platemotion -plates NA/EU -model MORVEL56 -o morvel56.tmp
echo   -30.000000000000000        40.000000000000000        22.771560776120207       -2.0222476466178430 > answer.tmp
./test_values.sh morvel56.tmp answer.tmp 4 "platemotion: MORVEL56 EU wrt NA, output file" || exit 1

# MORVEL
echo -75 -35 | ../bin/platemotion -plates SA/NZ -model MORVEL > morvel.tmp
echo -75.000000000000000       -35.000000000000000        72.008593269570923        16.689689796741547 > answer.tmp
./test_values.sh morvel.tmp answer.tmp 4 "platemotion: MORVEL SA wrt NZ, stdout" || exit 1

# NUVEL1A
echo 115 -50 > loc.tmp
echo 120 -50 >> loc.tmp
../bin/platemotion -f loc.tmp -plates AN/AU -model NUVEL1A > nuvel1a.tmp
cat > answer.tmp << EOF
   115.00000000000000       -50.000000000000000        22.881201729386667        68.378823029969183     
   120.00000000000000       -50.000000000000000        18.268958516527089        69.521449139059513 
EOF
./test_values.sh nuvel1a.tmp answer.tmp 4 "platemotion: NUVEL1A AU wrt AN, stdout" || exit 1

# ITRF08
../bin/platemotion -f loc.tmp -plates AN/AU -model ITRF08 > itrf08.tmp
cat > answer.tmp << EOF
   115.00000000000000       -50.000000000000000        25.766441315239788        68.326351344076343     
   120.00000000000000       -50.000000000000000        21.150858988809414        69.673756121784024  
EOF
./test_values.sh itrf08.tmp answer.tmp 4 "platemotion: ITRF08 AU wrt AN, stdout" || exit 1

# Custom pole
echo -30 40 | ../bin/platemotion -pole 139.461/61.796/0.211 -o custom.tmp
echo   -30.000000000000000        40.000000000000000        22.846411219358291       -2.0281377828526641 > answer.tmp 
./test_values.sh custom.tmp answer.tmp 4 "platemotion: Custom pole, file" || exit 1

# Clean up
rm *.tmp
