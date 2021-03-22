#!/bin/bash

#####
#	SET PATH TO TEST_VALUES SCRIPT
#####
TEST_BIN_DIR=$(echo $0 | xargs dirname)


#####
#	SET PATH TO HDEF EXECUTABLE
#####
# Check for o92util
$TEST_BIN_DIR/test_find_hdef_exec.sh fltinv || { echo "$0: could not find fltinv; exiting" 1>&2; exit 1; }
BIN_DIR=$(cat hdefexec.tmp | xargs dirname)


#####
#	RUN TEST
#####
trap "rm -f *.tmp anneal.log" 0 1 2 3 8 9




echo "---------------------------------------------------------"
echo "Test #1: 1 strike-slip fault, 4 3-component displacements"
echo "----------"
# Input fault slip (Cartesian coordinate system)
X=0     # km
Y=0     # km
Z=10    # km
STR=0
DIP=90
RAK=0
SLIP=1  # m
WID=4   # km
LEN=6   # km
echo $X $Y $Z $STR $DIP $RAK $SLIP $WID $LEN > o92_flt.tmp

# Station locations
cat > o92_sta.tmp << EOF
-1 -1  0
 1 -1  0
 1  1  0
-1  1  0
EOF

# Calculate "observed" displacements with no noise
# $BIN_DIR/o92util -flt o92_flt.tmp -sta o92_sta.tmp -disp o92_disp.tmp -xy
cat > o92_disp.tmp << EOF
  -1.0000000000000000       -1.0000000000000000        0.0000000000000000       -5.5615559982043607E-004  -8.8477258863561543E-004   8.4543975578312263E-004
   1.0000000000000000       -1.0000000000000000        0.0000000000000000       -5.5609930276229918E-004   8.8517303568264017E-004  -8.4645683166874525E-004
   1.0000000000000000        1.0000000000000000        0.0000000000000000        5.5609930275325227E-004   8.8517303568264017E-004   8.4645683166874590E-004
  -1.0000000000000000        1.0000000000000000        0.0000000000000000        5.5615559984757679E-004  -8.8477258863563311E-004  -8.4543975578312989E-004
EOF

# Prepare displacement observation and fault geometry files for fltinv
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$6}' o92_disp.tmp > fltinv_disp.tmp
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$8*1e3,$9*1e3}' o92_flt.tmp > fltinv_flt.tmp

# Actual solution
#echo "Actual solution (strike-slip dip-slip):"
#awk '{printf("%16.8e%16.8e\n"),$7*cos($6*0.0174533),$7*sin($6*0.0174533)}' o92_flt.tmp

# Linear least squares solution
#echo "Least squares solution:"
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -o inversion.tmp || exit 1
echo "1.00000000E+00 -1.29382228E-16" > answer.tmp
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 1 fault, three-component disp" || exit 1

#echo ----------
#echo Finished Test #1
#echo ----------
#echo



echo "----------------------------------------------------------"
echo "Test #2: 4 strike-slip faults, 9 3-component displacements, inversion constraints"
echo "----------"
# Input fault slip
#  X  Y  Z STR DIP RAK SLIP WID LEN
cat > o92_flt.tmp << EOF
   4  0  4  90  90 170  2.0   6   8
   4  0 10  90  90 180  1.5   6   8
  -4  0  4  90  90 180  1.5   6   8
  -4  0 10  90  90 190  1.0   6   8
EOF

# Station locations
# $BIN_DIR/grid -x -6 6 -dx 6 -y -6 6 -dy 6 -z 0.0 -o o92_sta.tmp

# Calculate "observed" displacements with no noise
# $BIN_DIR/o92util -flt o92_flt.tmp -sta o92_sta.tmp -disp o92_disp.tmp -xy
cat > o92_disp.tmp << EOF
  -6.0000000000000000       -6.0000000000000000        0.0000000000000000      -0.20035016250035553      -0.15604240608572359        4.5958255739687107E-002
  -6.0000000000000000        0.0000000000000000        0.0000000000000000       -6.5466781066928636E-005  -7.8233895879369450E-002   1.6995235867557408E-005
  -6.0000000000000000        6.0000000000000000        0.0000000000000000       0.20033986533744855      -0.15603461936436319       -4.5948812017554172E-002
   0.0000000000000000       -6.0000000000000000        0.0000000000000000      -0.24311151977124124       -5.8490968107925534E-002   2.4734734081083933E-002
   0.0000000000000000        0.0000000000000000        0.0000000000000000       -1.0249819390403436E-004  -2.0293070361494761E-002   3.1741048190721941E-005
   0.0000000000000000        6.0000000000000000        0.0000000000000000       0.24310394708811836       -5.8493320065453450E-002  -2.4733618060377799E-002
   6.0000000000000000       -6.0000000000000000        0.0000000000000000      -0.18973221194885226       0.12131519562731848       -3.2253317716251804E-002
   6.0000000000000000        0.0000000000000000        0.0000000000000000       -7.6154061877362774E-005   7.8201785507845886E-002   1.5575775065860158E-005
   6.0000000000000000        6.0000000000000000        0.0000000000000000       0.18972281853858350       0.12130148637419379        3.2242003539650503E-002
EOF

# Prepare displacement observation and fault geometry files for fltinv
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$6}' o92_disp.tmp > fltinv_disp.tmp
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$8*1e3,$9*1e3}' o92_flt.tmp > fltinv_flt.tmp

# Actual solution
#echo "Actual solution (strike-slip dip-slip):"
#awk '{printf("%16.8e%16.8e\n"),$7*cos($6*0.0174533),$7*sin($6*0.0174533)}' o92_flt.tmp

# Linear least squares solutions
#echo "Least squares solution:"
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
 -1.96961551E+00  3.47296355E-01
 -1.50000000E+00 -4.38361161E-10
 -1.50000000E+00 -1.04931971E-10
 -9.84807752E-01 -1.73647971e-01
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, three-component disp" || exit 1

#echo Least squares solution + only horizontal displacements:
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -disp:components 12 \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
 -1.96961560E+00  3.47295266E-01
 -1.49999504E+00  3.26290455E-06
 -1.49999990E+00 -1.08904575E-06
 -9.84812714E-01 -1.73644915E-01
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, horizontal disp" || exit 1

#echo Least squares solution + only vertical displacements:
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -disp:components 3 \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
 -5.92420991E+00  5.83424688E-01
  2.05865282E+00 -5.30728300E-01
 -5.45459440E+00 -2.36128333E-01
  2.57384507E+00  3.57080123E-01
EOF
#$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, vertical disp" || exit 1

#echo "Least squares solution + fixed rake (actual rake value):"
cat > rake.tmp << EOF
170
180
180
190
EOF
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -flt:rake rake.tmp \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
 -1.96961551E+00  3.47296355E-01
 -1.50000000E+00  1.83697020E-16
 -1.50000000E+00  1.83697020E-16
 -9.84807753E-01 -1.73648178E-01
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, fixed rake (input value)" || exit 1

#echo "Least squares solution + fixed rake (all 180) + gels:"
cat > rake.tmp << EOF
180
180
180
180
EOF
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -flt:rake rake.tmp \
    -lsqr:mode gels \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
 -2.18026313E+00  2.67005226E-16
 -1.25371674E+00  1.53536020E-16
 -1.72480941E+00  2.11228232E-16
  7.35844003E-01 -9.01149003E-17
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, fixed rake (180)" || exit 1

#echo "Least squares solution + fixed rake (all 180) + nnls:"
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -flt:rake rake.tmp \
    -lsqr:mode nnls \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
 -2.13207226e+00  2.61238824E-16
 -1.37503158e+00  1.67855500E-16
 -1.55499585e+00  1.90464178E-16
 -0.00000000E+00  0.00000000E+00
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, fixed rake (180), nnls" || exit 1

#echo "Least squares solution + rotated rakes (135,225) + nnls:"
cat > rake.tmp << EOF
135 225
135 225
135 225
135 225
EOF
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -flt:rake rake.tmp \
    -lsqr:mode nnls \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
 -1.96961551E+00  3.47296355E-01
 -1.50000000E+00 -4.38365233E-10
 -1.50000000E+00 -1.04934617E-10
 -9.84807752E-01 -1.73647971e-01
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, rotated rakes, nnls" || exit 1

#echo "Least squares solution + 2 fixed slip magnitudes (ss1=-2.0, ss4=-1.0):"
cat > slip.tmp << EOF
 -2.0 99999
99999 99999
99999 99999
 -1.0 99999
EOF
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -flt:slip slip.tmp \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
 -2.00000000E+00  3.59056487E-01
 -1.37164031E+00 -3.09765641E-02
 -1.50515200E+00 -1.60242221E-02
 -1.00000000E+00 -1.17142278E-01
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, 2 fixed slip components" || exit 1

#echo Least squares solution + damping = 0.1:
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -damp 0.1 \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
 -2.08579403E+00  4.00222223E-01
 -1.02542077E+00 -2.01476020E-01
 -1.59764763E+00 -5.74055238E-02
 -6.61180279E-01  5.77337231E-02
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, damping=0.1" || exit 1

#echo Least squares solution + damping = 1.0:
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -damp 1.0 \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
 -9.42573327E-02  3.08437760E-03
 -2.42517153E-02  8.46777249E-03
 -7.77603451E-02  2.83108205E-02
 -1.96980984E-02  7.18217238E-03
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, damping=1.0" || exit 1

#echo Least squares solution + smoothing = 0.1:
cat > smooth.tmp << EOF
1 2 2 3
2 2 1 4
3 2 1 4
4 2 2 3
EOF
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -smooth 0.1 smooth.tmp \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
 -1.91156100E+00  3.24048634E-01
 -1.64834091E+00  1.07707954E-01
 -1.43189167E+00  1.10122449E-02
 -1.30272323E+00 -2.30797720E-01
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, smoothing=0.1" || exit 1

#echo Least squares solution + smoothing = 1.0:
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -smooth 1.0 smooth.tmp \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
 -1.59988382E+00  1.33773112E-01
 -1.59933238E+00  1.33078461E-01
 -1.59830262E+00  1.31169871E-01
 -1.59827693E+00  1.30915893E-01
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, smoothing=1.0" || exit 1

#echo Least squares solution + damping + smoothing:
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -damp 0.05 \
    -smooth 0.05 smooth.tmp \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
 -1.95967192E+00  3.35790450E-01
 -1.51533042E+00  3.93285795E-02
 -1.48272842E+00  9.86849084E-03
 -1.05891374E+00 -2.01470950E-01
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, damping=0.05 smoothing=0.05" || exit 1

#echo ----------
#echo Finished Test #2
#echo ----------
#echo



echo "---------------------------------------------------------"
echo "Test #3: 4 strike-slip faults, 16 line-of-sight displacements"
echo "----------"
# Input fault slip
#  X  Y  Z STR DIP RAK SLIP WID LEN
cat > o92_flt.tmp << EOF
   4  0  4  90  90 170  2.0   6   8
   4  0 10  90  90 180  1.5   6   8
  -4  0  4  90  90 180  1.5   6   8
  -4  0 10  90  90 190  1.0   6   8
EOF

# Station locations
# $BIN_DIR/grid -x -6 6 -dx 4 -y -6 6 -dy 4 -z 0.0 -o o92_sta.tmp

# Calculate "observed" displacements with no noise
# $BIN_DIR/o92util -flt o92_flt.tmp -sta o92_sta.tmp -disp o92_disp.tmp -xy
cat > o92_disp.tmp << EOF
  -6.0000000000000000       -6.0000000000000000        0.0000000000000000      -0.20035016250035553      -0.15604240608572359        4.5958255739687107E-002
  -6.0000000000000000       -2.0000000000000000        0.0000000000000000      -0.26189759414519564      -0.13133088717039326        6.7471905920287845E-002
  -6.0000000000000000        2.0000000000000000        0.0000000000000000       0.26186611946055371      -0.13131779209913871       -6.7450449634507859E-002
  -6.0000000000000000        6.0000000000000000        0.0000000000000000       0.20033986533744855      -0.15603461936436319       -4.5948812017554172E-002
  -2.0000000000000000       -6.0000000000000000        0.0000000000000000      -0.23985351202946126       -9.2105374442862292E-002   3.0984503618363555E-002
  -2.0000000000000000       -2.0000000000000000        0.0000000000000000      -0.38644233249399279       -7.0442500735258512E-002   4.4543479138611647E-002
  -2.0000000000000000        2.0000000000000000        0.0000000000000000       0.38639587672384712       -7.0436143752714203E-002  -4.4531981196420461E-002
  -2.0000000000000000        6.0000000000000000        0.0000000000000000       0.23984525795904715       -9.2102913348129201E-002  -3.0979960972422095E-002
   2.0000000000000000       -6.0000000000000000        0.0000000000000000      -0.22948586559377421       -1.1293533139646109E-002   1.2715759692174429E-002
   2.0000000000000000       -2.0000000000000000        0.0000000000000000      -0.39986884243685383       -6.3526843018681597E-002   8.7345818004774436E-002
   2.0000000000000000        2.0000000000000000        0.0000000000000000       0.39982668366087426       -6.3527086435346464E-002  -8.7346014419123263E-002
   2.0000000000000000        6.0000000000000000        0.0000000000000000       0.22947884231575769       -1.1300990219775549E-002  -1.2718750521928175E-002
   6.0000000000000000       -6.0000000000000000        0.0000000000000000      -0.18973221194885226       0.12131519562731848       -3.2253317716251804E-002
   6.0000000000000000       -2.0000000000000000        0.0000000000000000      -0.28469284590954791        7.4762110557563599E-002   1.1824895305332858E-004
   6.0000000000000000        2.0000000000000000        0.0000000000000000       0.28466441892087074        7.4745077966454576E-002  -1.4316941454659178E-004
   6.0000000000000000        6.0000000000000000        0.0000000000000000       0.18972281853858350       0.12130148637419379        3.2242003539650503E-002
EOF

# Prepare displacement observation and fault geometry files for fltinv
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$6}' o92_disp.tmp > fltinv_disp.tmp
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$8*1e3,$9*1e3}' o92_flt.tmp > fltinv_flt.tmp

# Calculate line-of-sight displacements
AZ="45"
INC="35"
# $BIN_DIR/vec2los -f o92_disp.tmp -o o92_los.tmp -a $AZ -i $INC
cat > o92_los.tmp << EOF
  -6.0000000000000000       -6.0000000000000000        0.0000000000000000      -0.23279311491233720
  -6.0000000000000000       -2.0000000000000000        0.0000000000000000      -0.26646922849756915
  -6.0000000000000000        2.0000000000000000        0.0000000000000000       0.11430523056767314
  -6.0000000000000000        6.0000000000000000        0.0000000000000000        5.2017993330452192E-002
  -2.0000000000000000       -6.0000000000000000        0.0000000000000000      -0.21005185155528006
  -2.0000000000000000       -2.0000000000000000        0.0000000000000000      -0.29018956237288968
  -2.0000000000000000        2.0000000000000000        0.0000000000000000       0.20855520833639754
  -2.0000000000000000        6.0000000000000000        0.0000000000000000       0.10334587327633581
   2.0000000000000000       -6.0000000000000000        0.0000000000000000      -0.14675962135377757
   2.0000000000000000       -2.0000000000000000        0.0000000000000000      -0.31851124305489031
   2.0000000000000000        2.0000000000000000        0.0000000000000000       0.24489374712544690
   2.0000000000000000        6.0000000000000000        0.0000000000000000       0.13366988895091877
   6.0000000000000000       -6.0000000000000000        0.0000000000000000       -2.1129306122340326E-002
   6.0000000000000000       -2.0000000000000000        0.0000000000000000      -0.12166557751326400
   6.0000000000000000        2.0000000000000000        0.0000000000000000       0.20826215020835376
   6.0000000000000000        6.0000000000000000        0.0000000000000000       0.16166072181469121
EOF

# Prepare displacement and fault files for fltinv
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,'"$AZ"','"$INC"'}' o92_los.tmp > fltinv_los.tmp

$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -los fltinv_los.tmp \
    -gf:model okada_rect \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
 -1.96961559E+00  3.47296349E-01
 -1.49999974E+00 -4.50134408E-07
 -1.50000043E+00 -1.39591787E-07
 -9.84805251E-01 -1.73647574E-01
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, 16 los disp" || exit 1

#echo ----------
#echo Finished Test #3
#echo ----------
#echo



echo "---------------------------------------------------------"
echo "Test #4: 4 strike-slip faults, 9 three-component, 16 line-of-sight displacements"
echo "----------"
# Input fault slip
#  X  Y  Z STR DIP RAK SLIP WID LEN
cat > o92_flt.tmp << EOF
   4  0  4  90  90 170  2.0   6   8
   4  0 10  90  90 180  1.5   6   8
  -4  0  4  90  90 180  1.5   6   8
  -4  0 10  90  90 190  1.0   6   8
EOF
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$8*1e3,$9*1e3}' o92_flt.tmp > fltinv_flt.tmp

# Station locations
# $BIN_DIR/grid -x -6 6 -dx 6 -y -6 6 -dy 6 -z 0.0 -o o92_sta_disp.tmp
# $BIN_DIR/grid -x -6 6 -dx 4 -y -6 6 -dy 4 -z 0.0 -o o92_sta_los.tmp

# Calculate "observed" displacements with no noise
# $BIN_DIR/o92util -flt o92_flt.tmp -sta o92_sta_disp.tmp -disp o92_disp.tmp -xy
cat > o92_disp.tmp << EOF
  -6.0000000000000000       -6.0000000000000000        0.0000000000000000      -0.20035016250035553      -0.15604240608572359        4.5958255739687107E-002
  -6.0000000000000000        0.0000000000000000        0.0000000000000000       -6.5466781066928636E-005  -7.8233895879369450E-002   1.6995235867557408E-005
  -6.0000000000000000        6.0000000000000000        0.0000000000000000       0.20033986533744855      -0.15603461936436319       -4.5948812017554172E-002
   0.0000000000000000       -6.0000000000000000        0.0000000000000000      -0.24311151977124124       -5.8490968107925534E-002   2.4734734081083933E-002
   0.0000000000000000        0.0000000000000000        0.0000000000000000       -1.0249819390403436E-004  -2.0293070361494761E-002   3.1741048190721941E-005
   0.0000000000000000        6.0000000000000000        0.0000000000000000       0.24310394708811836       -5.8493320065453450E-002  -2.4733618060377799E-002
   6.0000000000000000       -6.0000000000000000        0.0000000000000000      -0.18973221194885226       0.12131519562731848       -3.2253317716251804E-002
   6.0000000000000000        0.0000000000000000        0.0000000000000000       -7.6154061877362774E-005   7.8201785507845886E-002   1.5575775065860158E-005
   6.0000000000000000        6.0000000000000000        0.0000000000000000       0.18972281853858350       0.12130148637419379        3.2242003539650503E-002
EOF
# $BIN_DIR/o92util -flt o92_flt.tmp -sta o92_sta_los.tmp -disp o92_disp_los.tmp -xy
AZ="45"
INC="35"
# $BIN_DIR/vec2los -f o92_disp_los.tmp -o o92_los.tmp -a $AZ -i $INC
cat > o92_los.tmp << EOF
  -6.0000000000000000       -6.0000000000000000        0.0000000000000000      -0.23279311491233720
  -6.0000000000000000       -2.0000000000000000        0.0000000000000000      -0.26646922849756915
  -6.0000000000000000        2.0000000000000000        0.0000000000000000       0.11430523056767314
  -6.0000000000000000        6.0000000000000000        0.0000000000000000        5.2017993330452192E-002
  -2.0000000000000000       -6.0000000000000000        0.0000000000000000      -0.21005185155528006
  -2.0000000000000000       -2.0000000000000000        0.0000000000000000      -0.29018956237288968
  -2.0000000000000000        2.0000000000000000        0.0000000000000000       0.20855520833639754
  -2.0000000000000000        6.0000000000000000        0.0000000000000000       0.10334587327633581
   2.0000000000000000       -6.0000000000000000        0.0000000000000000      -0.14675962135377757
   2.0000000000000000       -2.0000000000000000        0.0000000000000000      -0.31851124305489031
   2.0000000000000000        2.0000000000000000        0.0000000000000000       0.24489374712544690
   2.0000000000000000        6.0000000000000000        0.0000000000000000       0.13366988895091877
   6.0000000000000000       -6.0000000000000000        0.0000000000000000       -2.1129306122340326E-002
   6.0000000000000000       -2.0000000000000000        0.0000000000000000      -0.12166557751326400
   6.0000000000000000        2.0000000000000000        0.0000000000000000       0.20826215020835376
   6.0000000000000000        6.0000000000000000        0.0000000000000000       0.16166072181469121
EOF

# Prepare displacement observation and files for fltinv
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$6}' o92_disp.tmp > fltinv_disp.tmp
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,'"$AZ"','"$INC"'}' o92_los.tmp > fltinv_los.tmp

$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -los fltinv_los.tmp \
    -gf:model okada_rect \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
 -1.96961540E+00  3.47296424E-01
 -1.49999953E+00 -3.83966638E-07
 -1.50000054E+00 -1.46434557E-07
 -9.84806772E-01 -1.73647802E-01
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, 9 3-comp disp, 16 los disp" || exit 1

#echo ----------
#echo Finished Test #4
#echo ----------
#echo


echo "----------------------------------------------------------"
echo "Test #5: 4 dip-slip faults, 25 3-component displacements, covariance"
echo "----------"
# Input fault slip
#  X  Y  Z STR DIP RAK SLIP WID LEN
cat > o92_flt.tmp << EOF
1.73205  3  1   0  30  65  2.0   4   6
5.19615  3  3   0  30 100  1.5   4   6
1.73205 -3  1   0  30  80  1.5   4   6
5.19615 -3  3   0  30 115  1.0   4   6
EOF
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$8*1e3,$9*1e3}' o92_flt.tmp > fltinv_flt.tmp
#echo "Actual solution (strike-slip dip-slip):"
#awk '{printf("%16.8e%16.8e\n"),$7*cos($6*0.0174533),$7*sin($6*0.0174533)}' o92_flt.tmp
#echo

# Station locations
# $BIN_DIR/grid -x -2 6 -nx 5 -y -7 7 -ny 5 -z 0.0 -o o92_sta.tmp

# Calculate "observed" displacements with small noise
# $BIN_DIR/o92util -flt o92_flt.tmp -sta o92_sta.tmp -disp o92_disp.tmp -xy
cat > o92_disp.tmp << EOF
  -2.0000000000000000       -7.0000000000000000        0.0000000000000000        9.2363264330245567E-002   1.1997320862263442E-002  -2.6759027063015094E-002
  -2.0000000000000000       -3.5000000000000000        0.0000000000000000       0.26563363177654353       -1.4098549560357256E-002  -1.9478349656136415E-002
  -2.0000000000000000        0.0000000000000000        0.0000000000000000       0.31359590532656739       -3.8874234337665070E-002  -2.2185088079829423E-002
  -2.0000000000000000        3.5000000000000000        0.0000000000000000       0.34257976047121147       -6.3892513570741480E-002  -3.4986197245589932E-003
  -2.0000000000000000        7.0000000000000000        0.0000000000000000       0.15858608694944923       -5.9003792733126677E-002  -6.4783443229684424E-003
   0.0000000000000000       -7.0000000000000000        0.0000000000000000        4.5771912399021186E-003  -5.5434008422655182E-002  -4.3851328742015777E-002
   0.0000000000000000       -3.5000000000000000        0.0000000000000000       -1.0978868052022694       0.18698321990661329       0.67922745926762684
   0.0000000000000000        0.0000000000000000        0.0000000000000000       -1.4164390793141619       0.29266260947764339       0.41134534913089860
   0.0000000000000000        3.5000000000000000        0.0000000000000000       -1.3304241268406902       0.72222663304432344       0.85806259330791590
   0.0000000000000000        7.0000000000000000        0.0000000000000000        4.3408212130634949E-002   8.2590750564958446E-002  -1.5277874870628167E-002
   2.0000000000000000       -7.0000000000000000        0.0000000000000000      -0.14807635651179163       -6.6910377700777152E-002   3.9175764679269615E-002
   2.0000000000000000       -3.5000000000000000        0.0000000000000000      -0.64956529903464455       0.10283404062936950       0.52582267290168272
   2.0000000000000000        0.0000000000000000        0.0000000000000000      -0.82202619691162360       0.29712592204123217       0.48615630689334566
   2.0000000000000000        3.5000000000000000        0.0000000000000000      -0.77571088711527281       0.52102669093855447       0.69404692889219433
   2.0000000000000000        7.0000000000000000        0.0000000000000000       -7.5960556910848756E-002  0.35273440070796258       0.15415135035315283
   4.0000000000000000       -7.0000000000000000        0.0000000000000000      -0.18523450564968194       -3.8557437048598389E-002   7.5709618597133668E-002
   4.0000000000000000       -3.5000000000000000        0.0000000000000000      -0.48500864302383345        8.4887554085128492E-003  0.23625675889031938
   4.0000000000000000        0.0000000000000000        0.0000000000000000      -0.62918571152584934        7.6439095752301123E-002  0.26090308831074188
   4.0000000000000000        3.5000000000000000        0.0000000000000000      -0.52837881264198572       0.18546542643434208       0.37919877618154901
   4.0000000000000000        7.0000000000000000        0.0000000000000000       -7.5033653108544082E-002  0.18351370478673845       0.10986768440522109
   6.0000000000000000       -7.0000000000000000        0.0000000000000000      -0.14395015301240430       -1.5037950166512673E-002   5.9034073156255043E-002
   6.0000000000000000       -3.5000000000000000        0.0000000000000000      -0.35616870327823746       -1.0196082753230890E-002   8.6092984192907188E-002
   6.0000000000000000        0.0000000000000000        0.0000000000000000      -0.45997116328091892       -2.3076039178626405E-002   7.9303836036420855E-002
   6.0000000000000000        3.5000000000000000        0.0000000000000000      -0.37977949011823853       -4.9292723514823000E-003  0.10568601091035872
   6.0000000000000000        7.0000000000000000        0.0000000000000000      -0.11876250050031774        1.3870294839854368E-002   3.8811230560353068E-004
EOF
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$6}' o92_disp.tmp > fltinv_disp.tmp

# Covariance matrix: iobs jobs icmp jcmp cov
cat > cov.tmp << EOF
0.01 0.01 0.04
0.02 0.01 0.05
0.01 0.02 0.03
0.02 0.02 0.08
0.01 0.02 0.05
0.01 0.01 0.04
0.02 0.01 0.03
0.02 0.01 0.05
0.02 0.02 0.06
0.01 0.01 0.04
0.02 0.01 0.05
0.01 0.02 0.03
0.02 0.02 0.08
0.01 0.02 0.05
0.01 0.01 0.04
0.02 0.01 0.03
0.02 0.01 0.05
0.02 0.02 0.06
0.02 0.01 0.03
0.02 0.01 0.05
0.02 0.02 0.06
0.01 0.01 0.04
0.02 0.01 0.05
0.01 0.02 0.03
0.02 0.02 0.08
EOF
awk '{print NR,NR,1,1,$1;print NR,NR,2,2,$2;print NR,NR,3,3,$3}' cov.tmp > j
mv j cov.tmp

cat > noise.tmp << EOF
-1.4152496444935702E-002 -1.2167028474564456E-003 -3.4635782241821289E-002
-8.4763468832385787E-004 7.5453025650004950E-004 -1.2063829950532134E-003
9.6945130095189933E-003 5.7152114352401430E-003 2.0121647387134782E-002
-2.6528798803991203E-003 1.9361163888658797E-002 2.7268521639765526E-002
-6.6578228558812820E-004 -1.6482822749079494E-002 -3.5189259417203009E-002
-1.0586362712237301E-003 -5.6500325397569303E-003 2.9778957366943363E-002
1.5334482399784789E-003 -2.7342840116851186E-003 -1.5969733802639708E-002
-5.4498953478676936E-003 8.5568838581746940E-004 -2.7436671816572850E-002
-7.9106420886759857E-003 -7.7488142616894788E-003 1.6351620153504973E-002
-3.3051888553463685E-003 6.3877203026596382E-003 8.3482303485578422E-004
-1.4367604742244799E-002 -7.8806167050283788E-004 -4.6090146108549473E-002
-2.3938993714293657E-003 4.3311392014123957E-004 9.2591178052279424E-003
-9.3238487535593460E-003 9.9525214457998473E-003 9.8477113915949465E-004
-9.0432519815406027E-003 1.0101622464705487E-002 4.2777417265639016E-002
-6.9264781718351405E-003 1.2579245530829138E-003 -2.0203347108802016E-002
-1.5681831508266683E-003 -4.1007816183323761E-003 -1.1826992643122770E-002
-8.2280915610644292E-003 2.1531838847666370E-003 5.1082676770735767E-002
-5.6215311799730572E-003 -1.5717176150302499E-003 3.1956130144547443E-002
1.4659825636416067E-002 -3.4483120757706313E-003 -1.6482395176984826E-002
-1.0679529637706523E-002 -3.3939806174258800E-003 -2.6680434844931780E-002
1.5619089408796661E-002 -1.0788891996656147E-002 -2.8875541930295984E-002
2.6976806776864187E-003 -1.3598107865878514E-003 1.4433164985812442E-002
6.9846425737653461E-003 5.9607180253583563E-004 9.9866441926177666E-003
-3.8827776300663854E-003 -9.1281533241271973E-003 6.0975697575783236E-003
9.3867809188609222E-003 -6.6163114139011934E-003 1.8654646922130976E-002
EOF
paste fltinv_disp.tmp noise.tmp |\
    awk '{print $1,$2,$3,$4+$7,$5+$8,$6+$9}' > j
mv j fltinv_disp.tmp

$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
  8.54029764E-01  1.82691412E+00
 -3.25610007E-01  1.46402442E+00
  2.56686931E-01  1.46217560E+00
 -3.69135504E-01  1.00378397E+00
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: 4 dip-slip faults, 25 observations, covariance matrix" || exit 1

#echo ----------
#echo Finished Test #5
#echo ----------
#echo



echo "----------------------------------------------------------------"
echo "Test #6: 1 strike-slip fault, pre-stresses from coincident fault"
echo "----------"
# Fault generating pre-stresses
#    X  Y  Z  STR DIP RAK SLIP WID LEN
echo 0  1 10    0  90   0    1   2   2 > flt.tmp

# Target fault
#  X  Y  Z STR DIP RAK SLIP WID LEN
echo 0  1 10 > sta.tmp

# Elastic half-space properties
echo shearmod 40e9 lame 40e9 > haf.tmp

# Calculate pre-stresses
# $BIN_DIR/o92util -flt flt.tmp -sta sta.tmp -stress stress.tmp -xy -haf haf.tmp
cat > stress.tmp << EOF
   0.0000000000000000        1.0000000000000000        10.000000000000000       -2.5723050739674535E-009   2.5723050739674535E-009   0.0000000000000000       -21004464.926200699        6.5406474328389001E-011  -1099.4169465314455
EOF

# Prepare fltinv input files
awk '{print $4,$5,$6,$7,$8,$9}' stress.tmp > j; mv j stress.tmp
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$8*1e3,$9*1e3}' flt.tmp > fltinv_flt.tmp

# Actual solution
#echo Actual solution
#echo -1 0

# Linear least squares solution
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -lsqr:mode gesv \
    -flt fltinv_flt.tmp \
    -prests stress.tmp \
    -gf:model okada_rect \
    -o inversion.tmp || exit 1
echo "-1.00000000E+00  0.00000000E+00" > answer.tmp
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 1 ss fault, pre-stresses from coincident fault" || exit 1

#echo ----------
#echo Finished Test #6
#echo ----------
#echo



echo "--------------------------------------------------------------"
echo "Test #7: 2 strike-slip faults, pre-stresses from central fault"
echo "----------"
# Fault generating pre-stresses
#    X  Y  Z STR DIP RAK SLIP WID LEN
echo "0  0 10   0  90   0    1   2   2" > flt.tmp

# Target faults
#    X   Y  Z STR DIP RAK SLIP WID LEN
echo "0   2 10   0  90   0    0   2   2" > sta.tmp
echo "0  -2 10   0  90   0    0   2   2" >> sta.tmp

# Elastic half-space properties
echo lame 40e9 shearmod 40e9 > haf.tmp

# Calculate pre-stresses
# $BIN_DIR/o92util -flt flt.tmp -sta sta.tmp -stress stress.tmp -xy -haf haf.tmp
cat > stress.tmp << EOF
   0.0000000000000000        2.0000000000000000        10.000000000000000        1.5667566622490962E-002  -9.2411234067602183E-002  0.31539102691881776        3572197.7065435522        1067.1251740383661        187.36879216043641
   0.0000000000000000       -2.0000000000000000        10.000000000000000       -1.5667567144590479E-002   9.2411238053089068E-002 -0.31539102688973752        3572197.7065435648       -1067.1251740382509        187.36879216043729
EOF

# Prepare fltinv input files
awk '{print $4,$5,$6,$7,$8,$9}' stress.tmp > j; mv j stress.tmp
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$8*1e3,$9*1e3}' sta.tmp > flt.tmp

# Actual solution
#echo Actual solution
#echo 0.173 0.000
#echo 0.173 0.000

# Linear least squares solution
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -lsqr:mode gesv \
    -flt flt.tmp \
    -prests stress.tmp \
    -gf:model okada_rect \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
  1.73043144E-01  6.59943338E-05
  1.73043144E-01 -6.59943338E-05
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "lsqr, minimize shear traction on 2 strike-slip faults" || exit 1

#echo ----------
#echo Finished Test #7
#echo ----------
#echo

#
#
#
#
#
#
#
#
#
#
#
#
#
#
rm *.tmp
#
#
#
#
#
#
#
#
#
#
#
#
#
#


if [ -f ../src/annealing_module.f90 ]
then
    echo "--------------------------------------------------------------"
    echo "Side Test: Are integer and double precision annealing algorithms the same?"
    echo "----------"
    awk 'BEGIN{p=0}{
        if (/^subroutine anneal_int_array/) {
            p = 1
        }
        if (p==1) {
            print $0
        }
        if (/^end subroutine anneal_int_array/) {
            p = 0
        }
    }' ../src/annealing_module.f90 > anneal_int_array.tmp
    awk 'BEGIN{p=0}{
        if (/^subroutine anneal_dp_array/) {
            p = 1
        }
        if (p==1) {
            print $0
        }
        if (/^end subroutine anneal_dp_array/) {
            p = 0
        }
    }' ../src/annealing_module.f90 > anneal_dp_array.tmp
    sed -e "s/dp/int/g" anneal_dp_array.tmp |\
        sed -e "/model/s/double precision/integer/g" > anneal_dp_array_mod.tmp
    PASSFAIL=`diff anneal_int_array.tmp anneal_dp_array_mod.tmp |\
        awk '{if(/[A-Z]/||/[a-z]/||/[0-9]/){print "FAIL";exit}}END{print "PASS"}'`
    if [ "$PASSFAIL" != "PASS" ]
    then
        echo "$0: WARNING: anneal_int_array() and anneal_dp_array() differ in more than variable type!"
        diff anneal_int_array.tmp anneal_dp_array.tmp
    fi
fi












echo "--------------------------------------------------------------"
echo "Test #8: 1 strike-slip fault, 4 3-component displacements, simulated annealing"
echo "----------"
# Input fault slip (Cartesian coordinate system)
X=0; Y=0; Z=10
STR=0; DIP=90; RAK=0
SLIP=1; WID=4; LEN=6
echo $X $Y $Z $STR $DIP $RAK $SLIP $WID $LEN > o92_flt.tmp

# Station locations
cat > o92_sta.tmp << EOF
-1 -1  0
 1 -1  0
 1  1  0
-1  1  0
EOF

# Calculate "observed" displacements with no noise
# $BIN_DIR/o92util -flt o92_flt.tmp -sta o92_sta.tmp -disp o92_disp.tmp -xy
cat > o92_disp.tmp << EOF
  -1.0000000000000000       -1.0000000000000000        0.0000000000000000       -5.5615559982043607E-004  -8.8477258863561543E-004   8.4543975578312263E-004
   1.0000000000000000       -1.0000000000000000        0.0000000000000000       -5.5609930276229918E-004   8.8517303568264017E-004  -8.4645683166874525E-004
   1.0000000000000000        1.0000000000000000        0.0000000000000000        5.5609930275325227E-004   8.8517303568264017E-004   8.4645683166874590E-004
  -1.0000000000000000        1.0000000000000000        0.0000000000000000        5.5615559984757679E-004  -8.8477258863563311E-004  -8.4543975578312989E-004
EOF

# Prepare displacement observation and fault geometry files for fltinv
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$6}' o92_disp.tmp > fltinv_disp.tmp
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$8*1e3,$9*1e3}' o92_flt.tmp > fltinv_flt.tmp

# Actual solution
#echo "Actual solution (strike-slip dip-slip):"
#awk '{printf("%16.8e%16.8e\n"),$7*cos($6*0.0174533),$7*sin($6*0.0174533)}' o92_flt.tmp

# Slip and rake constraints
echo "0 10" > slip.tmp
echo "-30 60" > rake.tmp
echo "0.25 1" > step.tmp

# Linear least squares solution
#echo "Least squares solution:"
$BIN_DIR/fltinv \
    -mode anneal \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -flt:slip slip.tmp \
    -flt:rake rake.tmp \
    -anneal:init_mode mean \
    -anneal:step step.tmp \
    -anneal:max_it 1000 \
    -anneal:reset_it 500 \
    -anneal:temp_start 0.5 \
    -anneal:temp_min 0.0 \
    -anneal:cool 0.98 \
    -anneal:log_file anneal.log \
    -anneal:seed 1 \
    -o inversion.tmp || exit 1
echo "1.18729953E+00 -1.69601249E-01" > answer.tmp
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: anneal, 1 fault, three-component disp" || exit 1

#echo ----------
#echo Finished Test #8
#echo ----------
#echo




echo "----------------------------------------------------------"
echo "Test #9: 4 dip-slip faults, 25 3-component displacements, covariance, annealing"
echo "----------"
# Input fault slip
#  X  Y  Z STR DIP RAK SLIP WID LEN
cat > o92_flt.tmp << EOF
1.73205  3  1   0  30  65  2.0   4   6
5.19615  3  3   0  30 100  1.5   4   6
1.73205 -3  1   0  30  80  1.5   4   6
5.19615 -3  3   0  30 115  1.0   4   6
EOF
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$8*1e3,$9*1e3}' o92_flt.tmp > fltinv_flt.tmp
#echo "Actual solution (strike-slip dip-slip):"
#awk '{printf("%16.8e%16.8e\n"),$7*cos($6*0.0174533),$7*sin($6*0.0174533)}' o92_flt.tmp
#echo

# Station locations
# $BIN_DIR/grid -x -2 6 -nx 5 -y -7 7 -ny 5 -z 0.0 -o o92_sta.tmp

# Calculate "observed" displacements with small noise
# $BIN_DIR/o92util -flt o92_flt.tmp -sta o92_sta.tmp -disp o92_disp.tmp -xy
cat > o92_disp.tmp << EOF
  -2.0000000000000000       -7.0000000000000000        0.0000000000000000        9.2363264330245567E-002   1.1997320862263442E-002  -2.6759027063015094E-002
  -2.0000000000000000       -3.5000000000000000        0.0000000000000000       0.26563363177654353       -1.4098549560357256E-002  -1.9478349656136415E-002
  -2.0000000000000000        0.0000000000000000        0.0000000000000000       0.31359590532656739       -3.8874234337665070E-002  -2.2185088079829423E-002
  -2.0000000000000000        3.5000000000000000        0.0000000000000000       0.34257976047121147       -6.3892513570741480E-002  -3.4986197245589932E-003
  -2.0000000000000000        7.0000000000000000        0.0000000000000000       0.15858608694944923       -5.9003792733126677E-002  -6.4783443229684424E-003
   0.0000000000000000       -7.0000000000000000        0.0000000000000000        4.5771912399021186E-003  -5.5434008422655182E-002  -4.3851328742015777E-002
   0.0000000000000000       -3.5000000000000000        0.0000000000000000       -1.0978868052022694       0.18698321990661329       0.67922745926762684
   0.0000000000000000        0.0000000000000000        0.0000000000000000       -1.4164390793141619       0.29266260947764339       0.41134534913089860
   0.0000000000000000        3.5000000000000000        0.0000000000000000       -1.3304241268406902       0.72222663304432344       0.85806259330791590
   0.0000000000000000        7.0000000000000000        0.0000000000000000        4.3408212130634949E-002   8.2590750564958446E-002  -1.5277874870628167E-002
   2.0000000000000000       -7.0000000000000000        0.0000000000000000      -0.14807635651179163       -6.6910377700777152E-002   3.9175764679269615E-002
   2.0000000000000000       -3.5000000000000000        0.0000000000000000      -0.64956529903464455       0.10283404062936950       0.52582267290168272
   2.0000000000000000        0.0000000000000000        0.0000000000000000      -0.82202619691162360       0.29712592204123217       0.48615630689334566
   2.0000000000000000        3.5000000000000000        0.0000000000000000      -0.77571088711527281       0.52102669093855447       0.69404692889219433
   2.0000000000000000        7.0000000000000000        0.0000000000000000       -7.5960556910848756E-002  0.35273440070796258       0.15415135035315283
   4.0000000000000000       -7.0000000000000000        0.0000000000000000      -0.18523450564968194       -3.8557437048598389E-002   7.5709618597133668E-002
   4.0000000000000000       -3.5000000000000000        0.0000000000000000      -0.48500864302383345        8.4887554085128492E-003  0.23625675889031938
   4.0000000000000000        0.0000000000000000        0.0000000000000000      -0.62918571152584934        7.6439095752301123E-002  0.26090308831074188
   4.0000000000000000        3.5000000000000000        0.0000000000000000      -0.52837881264198572       0.18546542643434208       0.37919877618154901
   4.0000000000000000        7.0000000000000000        0.0000000000000000       -7.5033653108544082E-002  0.18351370478673845       0.10986768440522109
   6.0000000000000000       -7.0000000000000000        0.0000000000000000      -0.14395015301240430       -1.5037950166512673E-002   5.9034073156255043E-002
   6.0000000000000000       -3.5000000000000000        0.0000000000000000      -0.35616870327823746       -1.0196082753230890E-002   8.6092984192907188E-002
   6.0000000000000000        0.0000000000000000        0.0000000000000000      -0.45997116328091892       -2.3076039178626405E-002   7.9303836036420855E-002
   6.0000000000000000        3.5000000000000000        0.0000000000000000      -0.37977949011823853       -4.9292723514823000E-003  0.10568601091035872
   6.0000000000000000        7.0000000000000000        0.0000000000000000      -0.11876250050031774        1.3870294839854368E-002   3.8811230560353068E-004
EOF

awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$6}' o92_disp.tmp > fltinv_disp.tmp

# Covariance matrix: iobs jobs icmp jcmp cov
cat > cov.tmp << EOF
0.01 0.01 0.04
0.02 0.01 0.05
0.01 0.02 0.03
0.02 0.02 0.08
0.01 0.02 0.05
0.01 0.01 0.04
0.02 0.01 0.03
0.02 0.01 0.05
0.02 0.02 0.06
0.01 0.01 0.04
0.02 0.01 0.05
0.01 0.02 0.03
0.02 0.02 0.08
0.01 0.02 0.05
0.01 0.01 0.04
0.02 0.01 0.03
0.02 0.01 0.05
0.02 0.02 0.06
0.02 0.01 0.03
0.02 0.01 0.05
0.02 0.02 0.06
0.01 0.01 0.04
0.02 0.01 0.05
0.01 0.02 0.03
0.02 0.02 0.08
EOF
awk '{print NR,NR,1,1,$1;print NR,NR,2,2,$2;print NR,NR,3,3,$3}' cov.tmp > j
mv j cov.tmp

cat > noise.tmp << EOF
-1.4152496444935702E-002 -1.2167028474564456E-003 -3.4635782241821289E-002
-8.4763468832385787E-004 7.5453025650004950E-004 -1.2063829950532134E-003
9.6945130095189933E-003 5.7152114352401430E-003 2.0121647387134782E-002
-2.6528798803991203E-003 1.9361163888658797E-002 2.7268521639765526E-002
-6.6578228558812820E-004 -1.6482822749079494E-002 -3.5189259417203009E-002
-1.0586362712237301E-003 -5.6500325397569303E-003 2.9778957366943363E-002
1.5334482399784789E-003 -2.7342840116851186E-003 -1.5969733802639708E-002
-5.4498953478676936E-003 8.5568838581746940E-004 -2.7436671816572850E-002
-7.9106420886759857E-003 -7.7488142616894788E-003 1.6351620153504973E-002
-3.3051888553463685E-003 6.3877203026596382E-003 8.3482303485578422E-004
-1.4367604742244799E-002 -7.8806167050283788E-004 -4.6090146108549473E-002
-2.3938993714293657E-003 4.3311392014123957E-004 9.2591178052279424E-003
-9.3238487535593460E-003 9.9525214457998473E-003 9.8477113915949465E-004
-9.0432519815406027E-003 1.0101622464705487E-002 4.2777417265639016E-002
-6.9264781718351405E-003 1.2579245530829138E-003 -2.0203347108802016E-002
-1.5681831508266683E-003 -4.1007816183323761E-003 -1.1826992643122770E-002
-8.2280915610644292E-003 2.1531838847666370E-003 5.1082676770735767E-002
-5.6215311799730572E-003 -1.5717176150302499E-003 3.1956130144547443E-002
1.4659825636416067E-002 -3.4483120757706313E-003 -1.6482395176984826E-002
-1.0679529637706523E-002 -3.3939806174258800E-003 -2.6680434844931780E-002
1.5619089408796661E-002 -1.0788891996656147E-002 -2.8875541930295984E-002
2.6976806776864187E-003 -1.3598107865878514E-003 1.4433164985812442E-002
6.9846425737653461E-003 5.9607180253583563E-004 9.9866441926177666E-003
-3.8827776300663854E-003 -9.1281533241271973E-003 6.0975697575783236E-003
9.3867809188609222E-003 -6.6163114139011934E-003 1.8654646922130976E-002
EOF
paste fltinv_disp.tmp noise.tmp |\
    awk '{print $1,$2,$3,$4+$7,$5+$8,$6+$9}' > j
mv j fltinv_disp.tmp

cat > slip.tmp << EOF
0 5
0 5
0 5
0 5
EOF
cat > rake.tmp << EOF
0 180
0 180
0 180
0 180
EOF
cat > step.tmp << EOF
0.05 2
0.05 2
0.05 2
0.05 2
EOF
$BIN_DIR/fltinv \
    -mode anneal \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -flt:slip slip.tmp \
    -flt:rake rake.tmp \
    -anneal:init_mode mean \
    -anneal:step step.tmp \
    -anneal:max_it 2000 \
    -anneal:reset_it 1000 \
    -anneal:temp_start 0.5 \
    -anneal:temp_min 0.0 \
    -anneal:cool 0.99 \
    -anneal:log_file anneal.log \
    -anneal:seed 1 \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
  8.66128699E-01  1.81351839E+00
 -3.67871230E-01  1.45673734E+00
  2.47535983E-01  1.45810405E+00
 -3.11359215E-01  1.03791256E+00
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: anneal, 4 dip-slip faults, 25 observations, covariance matrix, annealing" || exit 1

#echo ----------
#echo Finished Test #9
#echo ----------
#echo



echo "----------------------------------------------------------"
echo "Test #10: 9 strike-slip faults, 25 three-component displacements, simulated annealing with pseudo-coupling"
echo "----------"
# Compute displacements corresponding to pseudo-coupling solution
# Input faults
#  X  Y  Z STR DIP RAK SLIP WID LEN
cat > o92_flt.tmp << EOF
  -3  0  3  90  90 180    0   3   3
   0  0  3  90  90 180    0   3   3
   3  0  3  90  90 180    0   3   3
  -3  0  6  90  90 180    0   3   3
   0  0  6  90  90 180    1   3   3
   3  0  6  90  90 180    0   3   3
  -3  0  9  90  90 180    0   3   3
   0  0  9  90  90 180    0   3   3
   3  0  9  90  90 180    0   3   3
EOF
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$8*1e3,$9*1e3}' o92_flt.tmp > fltinv_flt.tmp

# Slip constraints
#   SS    DS
cat > fltinv_slip.tmp << EOF
 99999 99999
 99999 99999
 99999 99999
 99999 99999
    -1     0
 99999 99999
 99999 99999
 99999 99999
 99999 99999
EOF

# Pre-stresses
awk '{print 0,0,0,0,0,0}' o92_flt.tmp > fltinv_sts.tmp

# Calculate fault slip surrounding central fault
echo vp 6800 vs 3926 dens 3000 > haf.tmp
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -flt:slip fltinv_slip.tmp \
    -haf haf.tmp \
    -prests fltinv_sts.tmp \
    -gf:model okada_rect \
    -lsqr:mode gesv \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
 -9.36180539E-02  2.56503115E-02
 -1.44914388E-01 -2.87198170E-11
 -9.36180539E-02 -2.56503116E-02
 -2.05038818E-01  4.53406048E-03
 -1.00000000E+00  0.00000000E+00
 -2.05038818E-01 -4.53406048E-03
 -8.70859154E-02 -1.45560101E-02
 -1.35046459E-01  2.16533838E-11
 -8.70859154E-02  1.45560101E-02
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, pre-stresses, fix central fault slip, 8 surrounding unlocked faults" || exit 1


# Calculate displacements at grid of stations around faults
paste o92_flt.tmp inversion.tmp |\
    awk '{print $1,$2,$3,$4,$5,atan2($11,$10)/0.0174533,sqrt($10*$10+$11*$11),$8,$9}' > o92_flt_psc.tmp
# $BIN_DIR/grid -x -4.5 5.2 -nx 5 -y -3.6 2.9 -ny 5 -z 0.0 -o o92_sta_psc.tmp
# $BIN_DIR/o92util -flt o92_flt_psc.tmp -sta o92_sta_psc.tmp -haf haf.tmp -disp o92_disp_psc.tmp -xy
cat > o92_disp_psc.tmp << EOF
  -4.5000000000000000       -3.6000000000000001        0.0000000000000000       -1.9225765032057085E-002  -1.8571483875364930E-002   1.4815552271441086E-002
  -4.5000000000000000       -1.9750000000000001        0.0000000000000000       -1.7719660115527523E-002  -1.4804150678294343E-002   1.5588704164818942E-002
  -4.5000000000000000      -0.35000000000000009        0.0000000000000000       -4.7144156696321176E-003  -6.8112451873538786E-003   4.8946645304318061E-003
  -4.5000000000000000        1.2749999999999999        0.0000000000000000        1.4097369998547785E-002  -1.1492895234051349E-002  -1.3347654591952782E-002
  -4.5000000000000000        2.8999999999999999        0.0000000000000000        1.9391327049255818E-002  -1.7556134398886532E-002  -1.5755714492428137E-002
  -2.0750000000000002       -3.6000000000000001        0.0000000000000000       -1.5054406314088040E-002  -1.3459124501625595E-002   1.1593021224209651E-002
  -2.0750000000000002       -1.9750000000000001        0.0000000000000000       -1.5462138193981697E-002  -1.1098796314641501E-002   1.3200921433491230E-002
  -2.0750000000000002      -0.35000000000000009        0.0000000000000000       -4.6108536661924483E-003  -4.4003678000469060E-003   4.3440716684123972E-003
  -2.0750000000000002        1.2749999999999999        0.0000000000000000        1.3033780860622965E-002  -8.4289114898991983E-003  -1.1590684768884399E-002
  -2.0750000000000002        2.8999999999999999        0.0000000000000000        1.5805484069779427E-002  -1.3014192815644570E-002  -1.2788462710735629E-002
  0.34999999999999964       -3.6000000000000001        0.0000000000000000       -1.2482633382155514E-002   2.5894489879697909E-003  -2.2790548231487306E-003
  0.34999999999999964       -1.9750000000000001        0.0000000000000000       -1.3376729395398771E-002   2.1716014696139116E-003  -2.6653992750765090E-003
  0.34999999999999964      -0.35000000000000009        0.0000000000000000       -4.2743318218549774E-003   8.2458077643462037E-004  -8.4955399361278847E-004
  0.34999999999999964        1.2749999999999999        0.0000000000000000        1.1631501154010524E-002   1.6303687116957270E-003   2.3366709002309799E-003
  0.34999999999999964        2.8999999999999999        0.0000000000000000        1.3277772089753448E-002   2.5283505502148459E-003   2.5465921557052268E-003
   2.7749999999999995       -3.6000000000000001        0.0000000000000000       -1.6533374140964379E-002   1.6273614817967699E-002  -1.3775680329480322E-002
   2.7749999999999995       -1.9750000000000001        0.0000000000000000       -1.6438695381172833E-002   1.3275727952619106E-002  -1.5384756422338868E-002
   2.7749999999999995      -0.35000000000000009        0.0000000000000000       -4.7284000013612763E-003   5.4193394620114044E-003  -4.9576374759494621E-003
   2.7749999999999995        1.2749999999999999        0.0000000000000000        1.3600927333566821E-002   1.0086815315434967E-002   1.3387375345697727E-002
   2.7749999999999995        2.8999999999999999        0.0000000000000000        1.7164973328333222E-002   1.5647932244086237E-002   1.5067106626164822E-002
   5.1999999999999993       -3.6000000000000001        0.0000000000000000       -1.9520066389670378E-002   1.7860395149456035E-002  -1.3837391226908928E-002
   5.1999999999999993       -1.9750000000000001        0.0000000000000000       -1.7437417282997258E-002   1.3774210325020604E-002  -1.3795641888845436E-002
   5.1999999999999993      -0.35000000000000009        0.0000000000000000       -4.4801051773024131E-003   6.7733875018869665E-003  -3.9472019429309477E-003
   5.1999999999999993        1.2749999999999999        0.0000000000000000        1.3656182514119789E-002   1.0665410560028841E-002   1.1379000708997049E-002
   5.1999999999999993        2.8999999999999999        0.0000000000000000        1.9440689477291561E-002   1.6638724067361787E-002   1.4428744143992013E-002
EOF
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$6}' o92_disp_psc.tmp > fltinv_disp_psc.tmp

# Search for fault slip using simulated annealing with pseudo-coupling
cat > fltinv_slip_psc.tmp << EOF
-1 0
-1 0
-1 0
-1 0
-1 0
-1 0
-1 0
-1 0
-1 0
EOF
rm inversion.tmp
$BIN_DIR/fltinv \
    -mode anneal-psc \
    -xy \
    -flt fltinv_flt.tmp \
    -flt:slip fltinv_slip_psc.tmp \
    -disp fltinv_disp_psc.tmp \
    -haf haf.tmp \
    -gf:model okada_rect \
    -anneal:max_it 1000 \
    -anneal:temp_start 0.5 \
    -anneal:init_mode rand0.5 \
    -anneal:log_file anneal.log \
    -anneal-psc:min_flip 1 \
    -anneal-psc:max_flip 2 \
    -anneal:seed 12345 \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
 -9.36180539E-02  2.56503115E-02
 -1.44914388E-01 -2.87198170E-11
 -9.36180539E-02 -2.56503116E-02
 -2.05038818E-01  4.53406048E-03
 -1.00000000E+00  0.00000000E+00
 -2.05038818E-01 -4.53406048E-03
 -8.70859154E-02 -1.45560101E-02
 -1.35046459E-01  2.16533838E-11
 -8.70859154E-02  1.45560101E-02
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: simulated annealing + pseudo-coupling" || exit 1

ANNEAL_POST=$(which $BIN_DIR/anneal_post)
if [ "$ANNEAL_POST" != "" ]
then
    echo "Found anneal_post; testing"
    grep Iteration anneal.log | awk '{print $6}' | head -20 > fit.tmp
    $BIN_DIR/anneal_post -f anneal.log -obj obj.tmp || exit 1
    awk '{print $2}' obj.tmp | head -20 > fit.tmp
    cat > answer_old.tmp << EOF
    -4.2736531432366512E-002
    -5.9343060813405443E-002
    -4.9868570712277106E-002
    -7.4047685432829843E-002
    -4.7958495759038808E-002
    -5.4701619915554729E-002
    -5.4749118857992350E-002
    -5.7047801792104898E-002
    -8.3495981770522393E-002
    -7.6468065589400114E-002
    -5.7047801792104898E-002
    -4.7958495759038808E-002
    -2.4825797010700944E-002
    -4.0353960414947730E-003
    -2.7057133608867173E-002
    -2.8900615810068201E-002
    -2.4102062485012607E-002
    -2.2437408725433487E-002
    -1.9555864370221204E-002
    -4.6156206640372330E-002
EOF
    cat > answer.tmp << EOF
    5.9347999999999998E-002
    4.9862999999999998E-002
    7.4039999999999995E-002
    4.7960000000000003E-002
    5.4713999999999999E-002
    5.4771000000000000E-002
    5.7063999999999997E-002
    8.3502000000000007E-002
    7.6480000000000006E-002
    5.7063999999999997E-002
    4.7960000000000003E-002
    2.4826000000000001E-002
    4.0322999999999999E-003
    2.7050000000000001E-002
    2.8895000000000001E-002
    2.4101000000000001E-002
    2.2432000000000001E-002
    1.9557000000000001E-002
    4.6156999999999997E-002
    5.2974000000000000E-002
EOF
    $TEST_BIN_DIR/test_values.sh fit.tmp answer.tmp 1 "fltinv: simulated annealing + pseudo-coupling, first 20 fits" || exit 1

    echo "Testing anneal_post with resampling"
    $BIN_DIR/anneal_post -f anneal.log -plocked plock.tmp -resample 1000 -seed 37 || exit 1
    cat > answer.tmp << EOF
   2.1000000000000001E-002
   9.1999999999999998E-002
  0.10400000000000000
  0.13500000000000001
  0.81699999999999995
   5.7000000000000002E-002
  0.18500000000000000
  0.18200000000000000
  0.20899999999999999
EOF
    $TEST_BIN_DIR/test_values.sh plock.tmp answer.tmp 1 "fltinv: simulated annealing + pseudo-coupling, resampled, probability locked" || exit 1
else
    echo "Could not locate anneal_post; not testing"
fi

#echo ----------
#echo Finished Test #10
#echo ----------
#echo



echo "----------------------------------------------------------"
echo "Test #11: 4 triangular strike-slip faults, 9 3-component displacements"
echo "----------"
# Input fault slip
# x1 y1 z1 x2 y2 z2 x3 y3 z3 ss ds ts
cat > tri_flt.tmp << EOF
76.3 35.8 2.0 76.4 35.7 5.0 72.7 32.5 2.0 1.8 0.6 0.0
76.9 35.3 6.0 76.4 35.7 5.0 72.7 32.5 2.0 1.3 0.4 0.0
76.9 35.3 6.0 72.4 29.1 4.8 72.7 32.5 2.0 0.7 1.1 0.0
71.0 29.9 1.5 72.4 29.1 4.8 72.7 32.5 2.0 0.2 1.4 0.0
EOF

# Station locations
# $BIN_DIR/grid -x 70.0 80.0 -nx 3 -y 30.0 40.0 -ny 3 -z 0.0 -o tri_sta.tmp

# Calculate "observed" displacements with no noise
# $BIN_DIR/triutil -flt tri_flt.tmp -sta tri_sta.tmp -disp tri_disp.tmp -xy
cat > tri_disp.tmp << EOF
   70.000000000000000        30.000000000000000        0.0000000000000000       -9.9634470672579842E-003   1.1158951018305340E-004  -7.3615236784321560E-003
   70.000000000000000        35.000000000000000        0.0000000000000000        1.3372107126126978E-002  -6.1252358203966913E-002  -4.0286377341635274E-002
   70.000000000000000        40.000000000000000        0.0000000000000000        2.1347297140414685E-002  -5.3316116651404670E-002  -1.7980506136206904E-002
   75.000000000000000        30.000000000000000        0.0000000000000000        7.1635071168204473E-002   2.3912374086484018E-002  0.10851703175324244
   75.000000000000000        35.000000000000000        0.0000000000000000        3.5569888615788470E-002   5.7455271348374692E-002   8.8318945992800574E-002
   75.000000000000000        40.000000000000000        0.0000000000000000        7.7017116313884918E-003  -7.1598610073826674E-002  -2.6569542442407012E-002
   80.000000000000000        30.000000000000000        0.0000000000000000        3.6591925791223792E-002   9.1695112123304882E-004   1.9638848360642325E-002
   80.000000000000000        35.000000000000000        0.0000000000000000       0.13428832189744305        2.4502373833784227E-002   7.4083584178240805E-002
   80.000000000000000        40.000000000000000        0.0000000000000000        2.2193803899185068E-002   1.0075141310700794E-002   8.7468694825564158E-003
EOF

# Prepare displacement observation and fault geometry files for fltinv
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$6}' tri_disp.tmp > fltinv_disp.tmp
awk '{c=1e3;print $1*c,$2*c,$3*c,$4*c,$5*c,$6*c,$7*c,$8*c,$9*c}' tri_flt.tmp > fltinv_flt.tmp

# Actual solution
#echo "Actual solution (strike-slip dip-slip):"
#awk '{printf("%16.8e%16.8e\n"),$10,$11}' tri_flt.tmp

# Linear least squares solutions
#echo "Least squares solution:"
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model triangle \
    -o inversion.tmp || exit 1
cat > answer.tmp << EOF
  1.80000000E+00  6.00000000E-01
  1.30000000E+00  4.00000000E-01
  7.00000000E-01  1.10000000E+00
  2.00000000E-01  1.40000000E+00
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, three-component disp" || exit 1

#echo ----------
#echo Finished Test #11
#echo ----------
#echo




echo "----------------------------------------------------------"
echo "Test #12: 1 triangular dip-slip fault, pre-stresses"
echo "----------"
# One triangular fault
#
#   *(0,10,0)
#   |   --__
#   |       *(15,0,5)
#   |   --
#   *(0,-10,0)
#
echo "0 -10 0 15 0 5 0 10 0" > tri.tmp

# Triangle center
awk '{print ($1+$4+$7)/3,($2+$5+$8)/3,($3+$6+$9)/3}' tri.tmp > center.tmp

# One meter normal slip
echo "0 -1 0" > slip.tmp

# Compute stress at triangle center
paste tri.tmp slip.tmp > triutil_flt.tmp
cp center.tmp triutil_sta.tmp
echo shearmod 40e9 lame 40e9 > haf.tmp
# $BIN_DIR/triutil \
#     -flt triutil_flt.tmp \
#     -sta triutil_sta.tmp \
#     -strain triutil_stn.tmp \
#     -stress triutil_sts.tmp \
#     -xy \
#     -haf haf.tmp
cat > triutil_sts.tmp << EOF
   5.0000000000000000        0.0000000000000000        1.6666700000000001       -3271536.3033847655       -1495673.3460660102       -204487.87081140187        5.3332305132811976E-010  -1009349.3938727777        4.9625255669716852E-011
EOF

# Prepare fltinv files
awk '{print $1*1e3,$2*1e3,$3*1e3,$4*1e3,$5*1e3,$6*1e3,$7*1e3,$8*1e3,$9*1e3}' tri.tmp > fltinv_flt.tmp
awk '{print $4,$5,$6,$7,$8,$9}' triutil_sts.tmp > fltinv_sts.tmp
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -gf:model triangle \
    -prests fltinv_sts.tmp \
    -haf haf.tmp \
    -o inversion.tmp || exit 1
echo   "7.97433674E-17  9.99999930E-01" > answer.tmp
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: triangular dip-slip fault with pre-stresses" || exit 1

#echo ----------
#echo Finished Test #12
#echo ----------
#echo



echo "----------------------------------------------------------"
echo "Test #13: 1 triangular dip-slip fault, pre-stresses"
echo "----------"
# One locked triangular fault, resolved onto adjacent unlocked fault
#
#   *(0,10,0)
#   |   --__
#   |       *(15,0,5)
#   |   --  |
#   *(0,-10,0)
#    -      |
#      --   |
#         --*(15,-20,5)
#
#
echo "0 -10 0 15 0 5 0 10 0" > tri1.tmp
echo "0 -10 0 15 0 5 15 -20 5" > tri2.tmp

# Triangle centers
awk '{print ($1+$4+$7)/3,($2+$5+$8)/3,($3+$6+$9)/3}' tri1.tmp > center1.tmp
awk '{print ($1+$4+$7)/3,($2+$5+$8)/3,($3+$6+$9)/3}' tri2.tmp > center2.tmp

# One meter normal slip on locked fault
echo "0 -1 0" > slip1.tmp

# Compute stress at triangle 2 center
paste tri1.tmp slip1.tmp > triutil_flt1.tmp
echo lame 40e9 shear 40e9 > haf.tmp
cp center2.tmp triutil_sta2.tmp
# $BIN_DIR/triutil \
#     -flt triutil_flt1.tmp \
#     -sta triutil_sta2.tmp \
#     -strain triutil_stn2.tmp \
#     -stress triutil_sts2.tmp \
#     -haf haf.tmp \
#     -xy
cat > triutil_sts2.tmp << EOF
   10.000000000000000       -10.000000000000000        3.3333300000000001        189358.49723258859       -90092.941388008534       -219101.55215006717        220279.70328027982        367749.96302672598       -221143.21426432792
EOF

# Prepare fltinv files
awk '{print $1*1e3,$2*1e3,$3*1e3,$4*1e3,$5*1e3,$6*1e3,$7*1e3,$8*1e3,$9*1e3}' tri2.tmp > fltinv_flt2.tmp
awk '{print $4,$5,$6,$7,$8,$9}' triutil_sts2.tmp > fltinv_sts2.tmp
echo vp 6800 vs 3926 dens 3000 > haf.tmp
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -lsqr:mode gesv \
    -flt fltinv_flt2.tmp \
    -gf:model triangle \
    -prests fltinv_sts2.tmp \
    -haf haf.tmp \
    -o inversion.tmp || exit 1
echo -4.57250048E-02 -1.53447831E-01 > answer.tmp
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, minimize shear traction on 1 other triangular fault" || exit 1

#echo ----------
#echo Finished Test #13
#echo ----------
#echo



echo "----------------------------------------------------------"
echo "Test #14: 1 triangular dip-slip fault, pre-stresses on itself, geographic coordinates"
echo "----------"
# One triangular fault, using geographic coordinates
#
#       *(3.00, 42.00, 2.0)
#      /    --__
#     /         *(3.05, 41.92, 5.0)
#    /      --
#    *(2.97, 41.85, 1.8)
#
echo "3.00 42.00 2.0" > pt1.tmp
echo "3.05 41.92 5.0" > pt2.tmp
echo "2.97 41.85 1.8" > pt3.tmp
paste pt1.tmp pt2.tmp pt3.tmp > tri.tmp

# Triangle center
awk '{printf("%.8f %.8f %.8f\n"),($1+$4+$7)/3,($2+$5+$8)/3,($3+$6+$9)/3}' tri.tmp > center.tmp

# One meter normal slip
echo "0 -1 0" > slip.tmp

# Compute stress at triangle center
paste tri.tmp slip.tmp > triutil_flt.tmp
cp center.tmp triutil_sta.tmp
echo lame 40e9 shear 40e9 > haf.tmp
# $BIN_DIR/triutil \
#     -flt triutil_flt.tmp \
#     -sta triutil_sta.tmp \
#     -strain triutil_stn.tmp \
#     -stress triutil_sts.tmp \
#     -haf haf.tmp \
#     -geo
cat > triutil_sts.tmp << EOF
   3.0066666700000000        41.923333329999998        2.9333333300000000       -5943270.3094363473       -605514.94146627421        4017984.6113139139        703467.33757494017       -2786515.3232462630        378545.85600292875
EOF

# Prepare fltinv files
awk '{print $1,$2,$3*1e3,$4,$5,$6*1e3,$7,$8,$9*1e3}' tri.tmp > fltinv_flt.tmp
awk '{print $4,$5,$6,$7,$8,$9}' triutil_sts.tmp > fltinv_sts.tmp
$BIN_DIR/fltinv \
    -mode lsqr \
    -flt fltinv_flt.tmp \
    -gf:model triangle \
    -prests fltinv_sts.tmp \
    -geo \
    -haf haf.tmp \
    -o inversion.tmp || exit 1
echo "5.06302972E-09  1.00000001E+00" > answer.tmp
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: triangular fault with self-prestress, geographic coordinates" || exit 1

rm *.tmp
#echo ----------
#echo Finished Test #14
#echo ----------
#echo




echo "----------------------------------------------------------"
echo "Test #15: 1 triangular dip-slip fault, pre-stresses on adjacent fault, geographic coordinates"
echo "----------"
# Two triangular faults, using geographic coordinates
#
# -73.0 -72.8
#   1*---* -31.7
#    |1 /|3
#    | / |
#    |/ 2|4
#   2*---* -32.0
#
echo "-73.0 -31.7 6.0" > pt1.tmp
echo "-73.0 -32.0 6.0" > pt2.tmp
echo "-72.8 -31.7 8.3" > pt3.tmp
echo "-72.8 -32.0 8.3" > pt4.tmp
paste pt1.tmp pt2.tmp pt3.tmp > tri1.tmp
paste pt2.tmp pt3.tmp pt4.tmp > tri2.tmp

# One meter normal slip
echo "0 -1 0" > slip.tmp
paste tri1.tmp slip.tmp > triutil_flt.tmp

# Compute stress at adjacent triangle center
awk '{printf("%.8f %.8f %.8f\n"),($1+$4+$7)/3,($2+$5+$8)/3,($3+$6+$9)/3}' tri2.tmp > center.tmp
cp center.tmp triutil_sta.tmp
echo lame 40e9 shear 40e9 > haf.tmp
# $BIN_DIR/triutil \
#     -flt triutil_flt.tmp \
#     -sta triutil_sta.tmp \
#     -strain triutil_stn.tmp \
#     -stress triutil_sts.tmp \
#     -haf haf.tmp \
#     -geo

# Prepare fltinv files
awk '{print $1,$2,$3*1e3,$4,$5,$6*1e3,$7,$8,$9*1e3}' tri2.tmp > fltinv_flt.tmp
# awk '{print $4,$5,$6,$7,$8,$9}' triutil_sts.tmp > fltinv_sts.tmp
# $BIN_DIR/fltinv \
#     -mode lsqr \
#     -flt fltinv_flt.tmp \
#     -gf:model triangle \
#     -prests fltinv_sts.tmp \
#     -geo \
#     -haf haf.tmp \
#     -o inversion.tmp
# echo inversion with pre-stresses
# cat inversion.tmp
cat > fltinv_sts.tmp << EOF
119598.06719676481 80436.466337221820 142008.94233743200 50158.102151294443 736087.74462110875 -94954.489191024870
EOF

cat tri1.tmp tri2.tmp | awk '{print $1,$2,$3*1e3,$4,$5,$6*1e3,$7,$8,$9*1e3}' > fltinv_flt.tmp
cat tri1.tmp tri2.tmp | awk '{print 0,0,0,0,0,0}' > fltinv_sts.tmp
cat > fltinv_slip.tmp << EOF
0 -1
99999 99999
EOF
$BIN_DIR/fltinv \
    -mode lsqr \
    -flt fltinv_flt.tmp \
    -flt:slip fltinv_slip.tmp \
    -gf:model triangle \
    -prests fltinv_sts.tmp \
    -geo \
    -haf haf.tmp \
    -o inversion.tmp || exit 1

cat > answer.tmp << EOF
  0.00000000E+00 -1.00000000E+00
 -2.89142383E-02 -2.93023235E-01
EOF
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: 2 triangular faults, one locked, geographic coordinates" || exit 1

rm *.tmp
#echo ----------
#echo Finished Test #15
#echo ----------
#echo




echo "----------------------------------------------------------"
echo "Test #16: Euler pole, least-squares"
echo "----------"
# Make a couple of points in North America
cat > coords.tmp << EOF
-124 41
-123 44
-122 46
-124 48
-120 43
-121 45
-121 47
EOF

# Calculate motion with respect to pole
# $BIN_DIR/platemotion -f coords.tmp -pole -101/24/0.42 |\
#     awk '{print $1,$2,0,$3/1e3,$4/1e3,0}' > vel.tmp
cat > vel.tmp << EOF
-124.00000000000000 41.000000000000000 0 -0.0114292 -0.0166703 0
-123.00000000000000 44.000000000000000 0 -0.0138149 -0.0159823 0
-122.00000000000000 46.000000000000000 0 -0.0154564 -0.0152895 0
-124.00000000000000 48.000000000000000 0 -0.0164749 -0.0166703 0
-120.00000000000000 43.000000000000000 0 -0.0136194 -0.0138901 0
-121.00000000000000 45.000000000000000 0 -0.0149171 -0.014592 0
-121.00000000000000 47.000000000000000 0 -0.0163661 -0.014592 0
EOF

# Build Euler pole input file
cat > euler.tmp << EOF
1 # One pole
# DOES NOT MATTER WHAT IS IN THIS LINE FOR LEAST-SQUARES
EOF
# Points to model with rigid-body rotations
awk '{print 1,3,NR}' vel.tmp >> euler.tmp

$BIN_DIR/fltinv \
    -mode lsqr \
    -geo \
    -disp vel.tmp \
    -disp:unit m/yr \
    -euler euler.tmp euler_out.tmp || exit 1

awk '{print $1}' euler_out.tmp > inversion.tmp
echo -101 > answer.tmp
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 1 "fltinv: Euler pole lon" -zero 100 || exit 1
awk '{print $2}' euler_out.tmp > inversion.tmp
echo 24 > answer.tmp
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 1 "fltinv: Euler pole lat" -zero 100 || exit 1
awk '{print $3}' euler_out.tmp > inversion.tmp
echo 0.42 > answer.tmp
$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 1 "fltinv: Euler pole vel" -zero 1.0 || exit 1

rm *.tmp
#echo ----------
#echo Finished Test #16
#echo ----------
#echo



echo "----------------------------------------------------------"
echo "Test #17: Annealing with pseudo-coupling, plus an Euler pole"
echo "----------"
# Stations near North Island, New Zealand
# $BIN_DIR/grid -x 173 177 -dx 1 -y -42 -37 -dy 1 -z 0 > sta.tmp

# Put 8 triangles near the Hikurangi subduction zone
cat > tri.tmp << EOF
177.10 -38.1 30 177.95 -38.4 15 176.3 -39.3 30
177.95 -38.4 15 177.20 -39.6 15 176.3 -39.3 30
177.95 -38.4 15 178.80 -38.7  0 177.2 -39.6 15
178.80 -38.7  0 178.10 -39.9  0 177.2 -39.6 15
175.50 -40.5 30 176.45 -40.8 15 176.3 -39.3 30
176.45 -40.8 15 176.30 -39.3 30 177.2 -39.6 15
176.45 -40.8 15 177.20 -39.6 15 177.4 -41.1  0
177.20 -39.6 15 178.10 -39.9  0 177.4 -41.1  0
EOF

# Calculate pseudo-coupling slip
awk '{print $1,$2,$3*1e3,$4,$5,$6*1e3,$7,$8,$9*1e3}' tri.tmp > fltinv_tri.tmp

# Slip constraints
#   SS    DS
cat > fltinv_slip.tmp << EOF
 99999 99999
 99999 99999
 99999 99999
     0   -30
 99999 99999
     0   -30
 99999 99999
 99999 99999
EOF

# Pre-stresses
awk '{print 0,0,0,0,0,0}' fltinv_tri.tmp > fltinv_sts.tmp

# Calculate fault slip surrounding central fault
echo vp 6800 vs 3926 dens 3000 > haf.tmp
# $BIN_DIR/fltinv \
#     -mode lsqr \
#     -geo \
#     -flt fltinv_tri.tmp \
#     -flt:slip fltinv_slip.tmp \
#     -haf haf.tmp \
#     -prests fltinv_sts.tmp \
#     -gf:model triangle \
#     -lsqr:mode gesv \
#     -o inversion.tmp

# Calculate locking generated velocity
# paste tri.tmp inversion.tmp > tri_slip.tmp
# $BIN_DIR/triutil -flt tri_slip.tmp -sta sta.tmp -haf haf.tmp -disp locking_vel.tmp

# Set Euler pole far north of points to give them E velocity
# $BIN_DIR/platemotion -f sta.tmp -pole 175.0/3.0/0.13 > euler_vel.tmp

# Superimpose Euler velocity and locking signal
# paste locking_vel.tmp euler_vel.tmp |\
#     awk '{print $1,$2,$3*1e3,$4+$9,$5+$10,$6}' > fltinv_vel.tmp
cat > fltinv_vel.tmp << EOF
173.00000000000000 -42.000000000000000 0 10.1015 -0.520803 3.8538993612826977E-002
173.00000000000000 -41.000000000000000 0 9.75289 -0.521407 2.0285413283549381E-002
173.00000000000000 -40.000000000000000 0 9.38642 -0.446009 -1.3514960981347135E-002
173.00000000000000 -39.000000000000000 0 9.18584 -0.345235 -3.8375537368887547E-002
173.00000000000000 -38.000000000000000 0 9.11329 -0.314503 -4.1351171582214852E-002
173.00000000000000 -37.000000000000000 0 9.04627 -0.339216 -3.2222137123447098E-002
174.00000000000000 -42.000000000000000 0 10.1108 -0.271607 6.8629726779668720E-002
174.00000000000000 -41.000000000000000 0 9.59931 -0.311408 6.0172816967284048E-002
174.00000000000000 -40.000000000000000 0 8.93802 -0.132101 1.2909641225707924E-002
174.00000000000000 -39.000000000000000 0 8.8282 0.0835566 -3.9408928249684783E-002
174.00000000000000 -38.000000000000000 0 8.959 0.0745014 -4.6587535042090789E-002
174.00000000000000 -37.000000000000000 0 9.00081 -0.015961 -3.2421977109257899E-002
175.00000000000000 -42.000000000000000 0 10.1648 0.00480915 0.11467680673651726
175.00000000000000 -41.000000000000000 0 9.34826 -0.178907 0.19103502503146674
175.00000000000000 -40.000000000000000 0 7.28245 0.481017 0.35906844551456119
175.00000000000000 -39.000000000000000 0 7.99185 0.866093 4.4469352456502798E-002
175.00000000000000 -38.000000000000000 0 8.74353 0.573843 -3.1620889986092249E-002
175.00000000000000 -37.000000000000000 0 8.96904 0.322094 -2.2435956100920418E-002
176.00000000000000 -42.000000000000000 0 10.2672 0.325449 0.15372385666880739
176.00000000000000 -41.000000000000000 0 9.28766 0.554412 0.27233781272856461
176.00000000000000 -40.000000000000000 0 1.80047 2.38877 2.4040853139460570
176.00000000000000 -39.000000000000000 0 5.94941 3.02275 0.98605937819115053
176.00000000000000 -38.000000000000000 0 8.46863 1.22477 6.4278684544963549E-002
176.00000000000000 -37.000000000000000 0 8.98579 0.641412 6.9803250494727132E-003
177.00000000000000 -42.000000000000000 0 10.3558 0.518946 0.11747858916929399
177.00000000000000 -41.000000000000000 0 6.4673 5.68822 -2.1930284425009554
177.00000000000000 -40.000000000000000 0 -5.25102 9.54512 -6.4581624770331363
177.00000000000000 -39.000000000000000 0 3.60728 3.91816 0.25849941437411417
177.00000000000000 -38.000000000000000 0 8.28458 2.17667 0.57110663602468292
177.00000000000000 -37.000000000000000 0 9.08963 0.848983 5.6447270858779486E-002
EOF

# Fault slip range
awk '{print 0,-30}' fltinv_tri.tmp > fltinv_slip_constraint.tmp
# Euler pole range
cat > euler.tmp << EOF
1
175.0 0.0 5000 0 0.2
EOF
awk '{print 1,3,NR}' fltinv_vel.tmp >> euler.tmp

$BIN_DIR/fltinv \
    -mode anneal-psc \
    -geo \
    -flt fltinv_tri.tmp \
    -flt:slip fltinv_slip_constraint.tmp \
    -disp fltinv_vel.tmp \
    -disp:unit mm/yr \
    -gf:model triangle \
    -anneal:init_mode rand \
    -anneal:seed 5267146 \
    -anneal:max_it 1000 \
    -anneal:temp_start 1 \
    -anneal:temp_min 1 \
    -anneal:log_file anneal.log \
    -euler euler.tmp euler_pole.tmp \
    -o inversion.tmp || exit 1
# $BIN_DIR/anneal_post -f anneal.log -obj obj.tmp

rm *.tmp
#echo ----------
#echo Finished Test #17
#echo ----------
#echo



exit









# Triangular faults, using geographic coordinates
#
#   *-*-*-*-*-*
#   |\|\|\|\|\|
#   *-*-*-*-*-*
#   |\|\|\|\|\|
#   *-*-*-*-*-*
#   |\|\|\|\|\|
#   *-*-*-*-*-*
#   |\|\|\|\|\|
#   *-*-*-*-*-*
#   |\|\|\|\|\|
#   *-*-*-*-*-*
#
# WID=4
# WID2=2
# LEN=6
# grid -x 2.1 18.1 -dx 4 -y 3 27 -dy 6 |\
#     awk '{print $1*cos(30*0.0174533),$2,$1*sin(30*0.0174533)}' |\
#     awk '{
#         print $1+2*cos(30*0.0174533),$2+3,$3+2*cos(30*0.0174533)
#         print $1+2*cos(30*0.0174533),$2-3,$3+2*cos(30*0.0174533)
#         print $1-2*cos(30*0.0174533),$2-3,$3-2*cos(30*0.0174533)
#         print $1-2*cos(30*0.0174533),$2+3,$3-2*cos(30*0.0174533)
#     }'
# exit
#echo -73.0 -31.7 6.0 > pt1.tmp
#echo -73.0 -32.0 6.0 > pt2.tmp
#echo -72.8 -31.7 8.3 > pt3.tmp
#echo -72.8 -32.0 8.3 > pt4.tmp
#paste pt1.tmp pt2.tmp pt3.tmp > tri1.tmp
#paste pt2.tmp pt3.tmp pt4.tmp > tri2.tmp
#
## One meter normal slip
#echo 0 -1 0 > slip.tmp
#paste tri1.tmp slip.tmp > triutil_flt.tmp
#
## Compute stress at adjacent triangle center
#awk '{printf("%.8f %.8f %.8f\n"),($1+$4+$7)/3,($2+$5+$8)/3,($3+$6+$9)/3}' tri2.tmp > center.tmp
#cp center.tmp triutil_sta.tmp
#echo lame 40e9 shear 40e9 > haf.tmp
#$BIN_DIR/triutil \
#    -flt triutil_flt.tmp \
#    -sta triutil_sta.tmp \
#    -strain triutil_stn.tmp \
#    -stress triutil_sts.tmp \
#    -haf haf.tmp \
#    -geo
#
## Prepare fltinv files
#awk '{print $1,$2,$3*1e3,$4,$5,$6*1e3,$7,$8,$9*1e3}' tri2.tmp > fltinv_flt.tmp
#awk '{print $4,$5,$6,$7,$8,$9}' triutil_sts.tmp > fltinv_sts.tmp
#$BIN_DIR/fltinv \
#    -mode lsqr \
#    -flt fltinv_flt.tmp \
#    -gf:model triangle \
#    -prests fltinv_sts.tmp \
#    -geo \
#    -haf haf.tmp \
#    -o inversion.tmp -v 1
#cat inversion.tmp
#
#cat tri1.tmp tri2.tmp | awk '{print $1,$2,$3*1e3,$4,$5,$6*1e3,$7,$8,$9*1e3}' > fltinv_flt.tmp
#cat tri1.tmp tri2.tmp | awk '{print 0,0,0,0,0,0}' > fltinv_sts.tmp
#cat > fltinv_slip.tmp << EOF
#0 -1
#99999 99999
#EOF
#$BIN_DIR/fltinv \
#    -mode lsqr \
#    -flt fltinv_flt.tmp \
#    -flt:slip fltinv_slip.tmp \
#    -gf:model triangle \
#    -prests fltinv_sts.tmp \
#    -geo \
#    -haf haf.tmp \
#    -o inversion.tmp -v 1
#cat inversion.tmp
#
#echo 5.06302972E-09  1.00000001E+00 > answer.tmp
#$TEST_BIN_DIR/test_values.sh inversion.tmp answer.tmp 2 "fltinv: triangular fault with self-prestress, geographic coordinates" || exit 1

#echo ----------
#echo Finished Test #16
#echo ----------
#echo
