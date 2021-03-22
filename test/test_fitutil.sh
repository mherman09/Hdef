#!/bin/bash

#####
#	SET PATH TO TEST_VALUES SCRIPT
#####
TEST_BIN_DIR=$(echo $0 | xargs dirname)


#####
#	SET PATH TO HDEF EXECUTABLE
#####
# Check for o92util
$TEST_BIN_DIR/test_find_hdef_exec.sh fitutil || { echo "$0: could not find fitutil; exiting" 1>&2; exit 1; }
BIN_DIR=$(cat hdefexec.tmp | xargs dirname)


#####
#	RUN TEST
#####
trap "rm -f *.tmp" 0 1 2 3 8 9


# Fit quadratic polynomial
C0="4.6"
C1="-1.2"
cat > xy.tmp << EOF
1.0000000000000000 3.4
2.0000000000000000 2.2
3.0000000000000000 1
4.0000000000000000 -0.2
5.0000000000000000 -1.4
EOF
cat > answer.tmp << EOF
   4.6000000000000014
  -1.2000000000000004
EOF
$BIN_DIR/fitutil -f xy.tmp -poly 1 -label_coeff 0 > poly_coeff.tmp
$TEST_BIN_DIR/test_values.sh poly_coeff.tmp answer.tmp 1 "fitutil: 1st order polynomial" || exit 1

# Fit sinusoid
PERIOD="6.3"
AMPLITUDE="7.7"
PHASE="2.1"
cat > xy.tmp << EOF
-16.000000000000000 -7.4004
-15.000000000000000 -2.22834
-14.000000000000000 4.98245
-13.000000000000000 7.63476
-12.000000000000000 3.30197
-11.000000000000000 -4.05181
-10.000000000000000 -7.69856
-9.0000000000000000 -4.30184
-8.0000000000000000 3.03067
-7.0000000000000000 7.59039
-6.0000000000000000 5.20561
-5.0000000000000000 -1.94183
-4.0000000000000000 -7.31267
-3.0000000000000000 -5.9931
-2.0000000000000000 0.809603
-1.0000000000000000 6.87159
0.0000000000000000 6.64671
1.0000000000000000 0.340703
2.0000000000000000 -6.27702
3.0000000000000000 -7.15185
4.0000000000000000 -1.4834
5.0000000000000000 5.54222
6.0000000000000000 7.49722
7.0000000000000000 2.59296
EOF
cat > answer.tmp << EOF
   7.6999993020947040
   2.0999999661540500
EOF
$BIN_DIR/fitutil -f xy.tmp -sin $PERIOD -label_coeff 0 -o sin_coeff.tmp
$TEST_BIN_DIR/test_values.sh sin_coeff.tmp answer.tmp 1 "fitutil: sinusoid" || exit 1

# Fit exponential with free exponential constant
A="4.6"
C="-1.2"
cat > xy.tmp << EOF
1.0000000000000000 1.38549
2.0000000000000000 0.417303
3.0000000000000000 0.125689
4.0000000000000000 0.0378568
5.0000000000000000 0.0114023
EOF
$BIN_DIR/fitutil -f xy.tmp -exp -label_coeff n > exp_coeff.tmp
cat > answer.tmp << EOF
 4.59998641E+00
-1.19999901E+00
EOF
$TEST_BIN_DIR/test_values.sh exp_coeff.tmp answer.tmp 1 "fitutil: exponential, find exponential constant" || exit 1

# Fit exponential with fixed exponential constant
$BIN_DIR/fitutil -f xy.tmp -exp $C -label_coeff n > exp_coeff.tmp
cat > answer.tmp << EOF
   4.5999901531556864
EOF
$TEST_BIN_DIR/test_values.sh exp_coeff.tmp answer.tmp 1 "fitutil: exponential, fixed exponential constant" || exit 1

# Fit linear + sinusoid + exponential, print predicted values
C0=-5.2
C1=0.6
A=1.5
T=4.2
PHASE=0.6
AEXP=4.2
CEXP=-0.7
cat > xy.tmp << EOF
-3.0000000000000000 28.3164
-2.0000000000000000 9.60982
-1.0000000000000000 1.48651
0.0000000000000000 -0.153036
1.0000000000000000 -1.21651
2.0000000000000000 -3.61728
3.0000000000000000 -4.28111
4.0000000000000000 -2.10017
5.0000000000000000 -0.611316
6.0000000000000000 -1.76296
7.0000000000000000 -2.46435
8.0000000000000000 -0.382067
9.0000000000000000 1.7037
10.000000000000000 1.02502
11.000000000000000 -0.0610226
12.000000000000000 1.56111
13.000000000000000 3.99766
14.000000000000000 3.84889
15.000000000000000 2.49988
16.000000000000000 3.55706
17.000000000000000 6.17427
EOF
$BIN_DIR/fitutil -f xy.tmp -poly 1 -sin $T -exp $CEXP -label_coeff 0 -pre pre.tmp > multi.tmp
cat > answer.tmp << EOF
  -5.1999990119741932
  0.59999992890469644
   1.4999997983110380
  0.60000024366333293
   4.1999987142469690
EOF
$TEST_BIN_DIR/test_values.sh multi.tmp answer.tmp 1 "fitutil: linear + sinusoid + exponential" || exit 1
