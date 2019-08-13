#!/bin/bash

# Displacement
echo -123 40.2 0 > sta.tmp
echo -122.93930377 40.34037337 10.6704 -122.88118251 40.34818434 19.3296 -122.86074437 40.25961846 19.3296 -6.73205e-06 10 0 > tri.tmp
../bin/triutil -flt tri.tmp -sta sta.tmp -disp disp.tmp
echo   -123.00000000000000        40.200000000000003        0.0000000000000000        9.6569791024380125E-003   2.0009010635221869E-004  -6.2046465664675754E-003 > answer.tmp
./test_values.sh disp.tmp answer.tmp 6 "triutil: disp" || exit 1

cat > sta.tmp << EOF
-123 40.2 0.0
-121 40.6 1.0
EOF
cat > tri.tmp << EOF
-122.93930377 40.34037337 10.6704 -122.88118251 40.34818434 19.3296 -122.86074437 40.25961846 19.3296 -6.73205e-06 10 0
-122.93930377 40.34037337 10.6704 -122.86074437 40.25961846 19.3296 -122.91878948 40.25180749 10.6704 -6.73205e-06 10 0
EOF
../bin/triutil -flt tri.tmp -sta sta.tmp -disp disp.tmp
cat > answer.tmp << EOF
  -123.00000000000000        40.200000000000003        0.0000000000000000        2.9966961388728952E-002   6.2629065689802088E-003  -2.2918691003340319E-002
  -121.00000000000000        40.600000000000001        1.0000000000000000       -5.4671675253995140E-003  -1.1078757892770959E-003   5.6537610980703784E-004
EOF
./test_values.sh disp.tmp answer.tmp 6 "triutil: disp, two outputs" || exit 1


# Strain
cat > tri.tmp << EOF
-122.93930377 40.34037337 10.6704 -122.88118251 40.34818434 19.3296 -122.86074437 40.25961846 19.3296 -1 10 0
-122.93930377 40.34037337 10.6704 -122.86074437 40.25961846 19.3296 -122.91878948 40.25180749 10.6704 -2 10 0
EOF
#triutil -flt tri.tmp -sta sta.tmp -strain strain.tmp
../bin/triutil -flt tri.tmp -sta sta.tmp -strain strain.tmp
cat > answer.tmp << EOF
  -123.00000000000000        40.200000000000003        0.0000000000000000       -2.0278065772745314E-005  -2.4945319983837185E-006   7.5908659237095965E-006  -1.2602271373136892E-005  -3.0443371261343393E-020  -5.9655329947237400E-021
  -121.00000000000000        40.600000000000001        1.0000000000000000        6.0244145766807270E-008  -1.8766841507615213E-008  -1.3823728862643908E-008   1.5015079437210632E-008  -5.9673505651432135E-010  -9.8567652500763217E-011
EOF
./test_values.sh strain.tmp answer.tmp 9 "triutil: strain" || exit 1


# Stress
#triutil -flt tri.tmp -sta sta.tmp -stress stress.tmp
../bin/triutil -flt tri.tmp -sta sta.tmp -stress stress.tmp
cat > answer.tmp << EOF
  -123.00000000000000        40.200000000000003        0.0000000000000000       -2229514.5357164023       -806831.83376747486       -9.6624717116355896E-009  -1008181.7098509513       -2.4354697009074713E-009  -4.7724263957789922E-010
  -121.00000000000000        40.600000000000001        1.0000000000000000        5925.6746772065080       -395.20430474729096       0.24470685041342222        1201.2063549768504       -47.738804521145710       -7.8854122000610571 
EOF
./test_values.sh stress.tmp answer.tmp 9 "triutil: stress" || exit 1


# Max shear stress
#triutil -flt tri.tmp -sta sta.tmp -estress estress.tmp
../bin/triutil -flt tri.tmp -sta sta.tmp -estress estress.tmp
cat > answer.tmp << EOF
  -123.00000000000000        40.200000000000003        0.0000000000000000        1513512.5208073591     
  -121.00000000000000        40.600000000000001        1.0000000000000000        3739.2521231165501
EOF
./test_values.sh estress.tmp answer.tmp 4 "triutil: max shear stress" || exit 1


# Normal traction
echo 15 80 -6 0.5 > trg.tmp
triutil -flt tri.tmp -sta sta.tmp -trg trg.tmp -normal normal.tmp
../bin/triutil -flt tri.tmp -sta sta.tmp -trg trg.tmp -normal normal.tmp
cat > answer.tmp << EOF
  -123.00000000000000        40.200000000000003        0.0000000000000000       -1580967.8173957772     
  -121.00000000000000        40.600000000000001        1.0000000000000000        4738.7840988880153 
EOF
./test_values.sh normal.tmp answer.tmp 4 "triutil: normal" || exit 1


# Shear traction
triutil -flt tri.tmp -sta sta.tmp -trg trg.tmp -shear shear.tmp
../bin/triutil -flt tri.tmp -sta sta.tmp -trg trg.tmp -shear shear.tmp
cat > answer.tmp << EOF
  -123.00000000000000        40.200000000000003        0.0000000000000000       -1232623.6797110836        1241807.6939554971     
  -121.00000000000000        40.600000000000001        1.0000000000000000        2655.0409356889468        2723.1861172609219 
EOF
./test_values.sh shear.tmp answer.tmp 5 "triutil: shear stress" || exit 1


# Coulomb stress change
triutil -flt tri.tmp -sta sta.tmp -trg trg.tmp -coul coul.tmp
../bin/triutil -flt tri.tmp -sta sta.tmp -trg trg.tmp -coul coul.tmp
cat > answer.tmp << EOF
  -123.00000000000000        40.200000000000003        0.0000000000000000       -2023107.5884089721     
  -121.00000000000000        40.600000000000001        1.0000000000000000        5024.4329851329549
EOF
./test_values.sh coul.tmp answer.tmp 4 "triutil: coulomb" || exit 1



# Half-space parameters
echo shear 60e9 lame 40e9 > haf.tmp
../bin/triutil -flt tri.tmp -sta sta.tmp -trg trg.tmp -haf haf.tmp -coul coul.tmp
cat > answer.tmp << EOF
  -123.00000000000000        40.200000000000003        0.0000000000000000       -2910810.4521128507     
  -121.00000000000000        40.600000000000001        1.0000000000000000        7004.0166289097569
EOF
./test_values.sh coul.tmp answer.tmp 4 "triutil: coulomb, half-space parameters" || exit 1


# Target parameters
cat > trg.tmp << EOF
15 80 -6 0.5
-5 55 -93 0.4
EOF
../bin/triutil -flt tri.tmp -sta sta.tmp -trg trg.tmp -coul coul.tmp -prog
cat > answer.tmp << EOF
  -123.00000000000000        40.200000000000003        0.0000000000000000       -2023107.5884089721     
  -121.00000000000000        40.600000000000001        1.0000000000000000        4460.4409590834784
EOF
./test_values.sh coul.tmp answer.tmp 4 "triutil: coulomb" || exit 1


# Cartesian coordinates
echo 0.5 0.6 -3 > sta.tmp
echo 0 0 6 1 0 5 1 1 7 -1 2 0 > tri.tmp
../bin/triutil -flt tri.tmp -sta sta.tmp -disp disp.tmp -xy
echo   0.50000000000000000       0.59999999999999998       -3.0000000000000000       -3.0663748564409593E-002  -6.0112932472348501E-002 -0.15985669393942395 > answer.tmp 
./test_values.sh disp.tmp answer.tmp 6 "triutil: disp, cartesian coordinates" || exit 1

#####
#	CLEAN UP
#####
rm *.tmp
