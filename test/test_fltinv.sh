#!/bin/bash

BIN_DIR=`./define_bin_dir.sh`

echo ---------------------------------------------------------
echo Test \#1: 1 strike-slip fault, 4 3-component displacements
echo ----------
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
$BIN_DIR/o92util -flt o92_flt.tmp -sta o92_sta.tmp -disp o92_disp.tmp -xy

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
    -o inversion.tmp
echo 1.00000000E+00 -1.29382228E-16 > answer.tmp
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 1 fault, three-component disp" || exit 1

#echo ----------
#echo Finished Test \#1
#echo ----------
#echo



echo ----------------------------------------------------------
echo Test \#2: 4 strike-slip faults, 9 3-component displacements, inversion constraints
echo ----------
# Input fault slip
#  X  Y  Z STR DIP RAK SLIP WID LEN
cat > o92_flt.tmp << EOF
   4  0  4  90  90 170  2.0   6   8
   4  0 10  90  90 180  1.5   6   8
  -4  0  4  90  90 180  1.5   6   8
  -4  0 10  90  90 190  1.0   6   8
EOF

# Station locations
$BIN_DIR/grid -x -6 6 -dx 6 -y -6 6 -dy 6 -z 0.0 -o o92_sta.tmp

# Calculate "observed" displacements with no noise
$BIN_DIR/o92util -flt o92_flt.tmp -sta o92_sta.tmp -disp o92_disp.tmp -xy

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
    -o inversion.tmp
cat > answer.tmp << EOF
 -1.96961551E+00  3.47296355E-01
 -1.50000000E+00 -4.38361161E-10
 -1.50000000E+00 -1.04931971E-10
 -9.84807752E-01 -1.73648177E-01
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, three-component disp" || exit 1

#echo Least squares solution + only horizontal displacements:
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -disp:components 12 \
    -o inversion.tmp
cat > answer.tmp << EOF
 -1.96961560E+00  3.47295266E-01
 -1.49999504E+00  3.26290455E-06
 -1.49999990E+00 -1.08904575E-06
 -9.84812714E-01 -1.73644915E-01
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, horizontal disp" || exit 1

#echo Least squares solution + only vertical displacements:
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -disp:components 3 \
    -o inversion.tmp
cat > answer.tmp << EOF
 -5.92420991E+00  5.83424688E-01
  2.05865282E+00 -5.30728300E-01
 -5.45459440E+00 -2.36128333E-01
  2.57384507E+00  3.57080123E-01
EOF
#./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, vertical disp" || exit 1

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
    -o inversion.tmp
cat > answer.tmp << EOF
 -1.96961551E+00  3.47296355E-01
 -1.50000000E+00  1.83697020E-16
 -1.50000000E+00  1.83697020E-16
 -9.84807753E-01 -1.73648178E-01
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, fixed rake (input value)" || exit 1

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
    -o inversion.tmp
cat > answer.tmp << EOF
 -2.18026313E+00  2.67005226E-16
 -1.25371674E+00  1.53536020E-16
 -1.72480941E+00  2.11228232E-16
  7.35844003E-01 -9.01149003E-17
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, fixed rake (180)" || exit 1

#echo "Least squares solution + fixed rake (all 180) + nnls:"
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -flt:rake rake.tmp \
    -lsqr:mode nnls \
    -o inversion.tmp
cat > answer.tmp << EOF
 -2.13317688E+00  2.61238824E-16
 -1.37064417E+00  1.67855500E-16
 -1.55525804E+00  1.90464178E-16
 -0.00000000E+00  0.00000000E+00
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, fixed rake (180), nnls" || exit 1

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
    -o inversion.tmp
cat > answer.tmp << EOF
 -1.96961551E+00  3.47296355E-01
 -1.50000000E+00 -4.38365233E-10
 -1.50000000E+00 -1.04934617E-10
 -9.84807752E-01 -1.73648177E-01
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, rotated rakes, nnls" || exit 1

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
    -o inversion.tmp
cat > answer.tmp << EOF
 -2.00000000E+00  3.59056487E-01
 -1.37164031E+00 -3.09765641E-02
 -1.50515200E+00 -1.60242221E-02
 -1.00000000E+00 -1.17142278E-01
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, 2 fixed slip components" || exit 1

#echo Least squares solution + damping = 0.1:
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -damp 0.1 \
    -o inversion.tmp
cat > answer.tmp << EOF
 -2.08579403E+00  4.00222223E-01
 -1.02542077E+00 -2.01476020E-01
 -1.59764763E+00 -5.74055238E-02
 -6.61180279E-01  5.77337231E-02
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, damping=0.1" || exit 1

#echo Least squares solution + damping = 1.0:
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -damp 1.0 \
    -o inversion.tmp
cat > answer.tmp << EOF
 -9.42573327E-02  3.08437760E-03
 -2.42517153E-02  8.46777249E-03
 -7.77603451E-02  2.83108205E-02
 -1.96980984E-02  7.18217238E-03
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, damping=1.0" || exit 1

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
    -o inversion.tmp
cat > answer.tmp << EOF
 -1.91156100E+00  3.24048634E-01
 -1.64834091E+00  1.07707954E-01
 -1.43189167E+00  1.10122449E-02
 -1.30272323E+00 -2.30797720E-01
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, smoothing=0.1" || exit 1

#echo Least squares solution + smoothing = 1.0:
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -smooth 1.0 smooth.tmp \
    -o inversion.tmp
cat > answer.tmp << EOF
 -1.59988382E+00  1.33773112E-01
 -1.59933238E+00  1.33078461E-01
 -1.59830262E+00  1.31169871E-01
 -1.59827693E+00  1.30915893E-01
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, smoothing=1.0" || exit 1

#echo Least squares solution + damping + smoothing:
$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -disp fltinv_disp.tmp \
    -gf:model okada_rect \
    -damp 0.05 \
    -smooth 0.05 smooth.tmp \
    -o inversion.tmp
cat > answer.tmp << EOF
 -1.95967192E+00  3.35790450E-01
 -1.51533042E+00  3.93285795E-02
 -1.48272842E+00  9.86849084E-03
 -1.05891374E+00 -2.01470950E-01
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, damping=0.05 smoothing=0.05" || exit 1

#echo ----------
#echo Finished Test \#2
#echo ----------
#echo


echo ---------------------------------------------------------
echo Test \#3: 4 strike-slip faults, 16 line-of-sight displacements
echo ----------
# Input fault slip
#  X  Y  Z STR DIP RAK SLIP WID LEN
cat > o92_flt.tmp << EOF
   4  0  4  90  90 170  2.0   6   8
   4  0 10  90  90 180  1.5   6   8
  -4  0  4  90  90 180  1.5   6   8
  -4  0 10  90  90 190  1.0   6   8
EOF

# Station locations
$BIN_DIR/grid -x -6 6 -dx 4 -y -6 6 -dy 4 -z 0.0 -o o92_sta.tmp

# Calculate "observed" displacements with no noise
$BIN_DIR/o92util -flt o92_flt.tmp -sta o92_sta.tmp -disp o92_disp.tmp -xy

# Prepare displacement observation and fault geometry files for fltinv
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$6}' o92_disp.tmp > fltinv_disp.tmp
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$8*1e3,$9*1e3}' o92_flt.tmp > fltinv_flt.tmp

# Calculate line-of-sight displacements
AZ="45"
INC="35"
$BIN_DIR/vec2los -f o92_disp.tmp -o o92_los.tmp -a $AZ -i $INC

# Prepare displacement and fault files for fltinv
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,'"$AZ"','"$INC"'}' o92_los.tmp > fltinv_los.tmp

$BIN_DIR/fltinv \
    -mode lsqr \
    -xy \
    -flt fltinv_flt.tmp \
    -los fltinv_los.tmp \
    -gf:model okada_rect \
    -o inversion.tmp
cat > answer.tmp << EOF
 -1.96961559E+00  3.47296349E-01
 -1.49999974E+00 -4.50134408E-07
 -1.50000043E+00 -1.39591787E-07
 -9.84805251E-01 -1.73647574E-01
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, 16 los disp" || exit 1

#echo ----------
#echo Finished Test \#3
#echo ----------
#echo


echo ---------------------------------------------------------
echo Test \#4: 4 strike-slip faults, 9 three-component, 16 line-of-sight displacements
echo ----------
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
$BIN_DIR/grid -x -6 6 -dx 6 -y -6 6 -dy 6 -z 0.0 -o o92_sta_disp.tmp
$BIN_DIR/grid -x -6 6 -dx 4 -y -6 6 -dy 4 -z 0.0 -o o92_sta_los.tmp

# Calculate "observed" displacements with no noise
$BIN_DIR/o92util -flt o92_flt.tmp -sta o92_sta_disp.tmp -disp o92_disp.tmp -xy
$BIN_DIR/o92util -flt o92_flt.tmp -sta o92_sta_los.tmp -disp o92_disp_los.tmp -xy
AZ="45"
INC="35"
$BIN_DIR/vec2los -f o92_disp_los.tmp -o o92_los.tmp -a $AZ -i $INC

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
    -o inversion.tmp
cat > answer.tmp << EOF
 -1.96961540E+00  3.47296424E-01
 -1.49999953E+00 -3.83966638E-07
 -1.50000054E+00 -1.46434557E-07
 -9.84806772E-01 -1.73647802E-01
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, 9 3-comp disp, 16 los disp" || exit 1

#echo ----------
#echo Finished Test \#4
#echo ----------
#echo


echo ----------------------------------------------------------
echo Test \#5: 4 dip-slip faults, 25 3-component displacements, covariance
echo ----------
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
$BIN_DIR/grid -x -2 6 -nx 5 -y -7 7 -ny 5 -z 0.0 -o o92_sta.tmp

# Calculate "observed" displacements with small noise
$BIN_DIR/o92util -flt o92_flt.tmp -sta o92_sta.tmp -disp o92_disp.tmp -xy
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
    -o inversion.tmp
cat > answer.tmp << EOF
  8.54029764E-01  1.82691412E+00
 -3.25610007E-01  1.46402442E+00
  2.56686931E-01  1.46217560E+00
 -3.69135504E-01  1.00378397E+00
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: 4 dip-slip faults, 25 observations, covariance matrix" || exit 1

#echo ----------
#echo Finished Test \#5
#echo ----------
#echo



echo ----------------------------------------------------------------
echo Test \#6: 1 strike-slip fault, pre-stresses from coincident fault
echo ----------
# Fault generating pre-stresses
#    X  Y  Z  STR DIP RAK SLIP WID LEN
echo 0  1 10    0  90   0    1   2   2 > flt.tmp

# Target fault
#  X  Y  Z STR DIP RAK SLIP WID LEN
echo 0  1 10 > sta.tmp

# Elastic half-space properties
echo shearmod 40e9 lame 40e9 > haf.tmp

# Calculate pre-stresses
$BIN_DIR/o92util -flt flt.tmp -sta sta.tmp -stress stress.tmp -xy -haf haf.tmp

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
    -o inversion.tmp
echo -1.00000000E+00  0.00000000E+00 > answer.tmp
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 1 ss fault, pre-stresses from coincident fault" || exit 1

#echo ----------
#echo Finished Test \#6
#echo ----------
#echo



echo --------------------------------------------------------------
echo Test \#7: 2 strike-slip faults, pre-stresses from central fault
echo ----------
# Fault generating pre-stresses
#    X  Y  Z STR DIP RAK SLIP WID LEN
echo 0  0 10   0  90   0    1   2   2 > flt.tmp

# Target faults
#    X   Y  Z STR DIP RAK SLIP WID LEN
echo 0   2 10   0  90   0    0   2   2 > sta.tmp
echo 0  -2 10   0  90   0    0   2   2 >> sta.tmp

# Elastic half-space properties
echo lame 40e9 shearmod 40e9 > haf.tmp

# Calculate pre-stresses
$BIN_DIR/o92util -flt flt.tmp -sta sta.tmp -stress stress.tmp -xy -haf haf.tmp

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
    -o inversion.tmp
cat > answer.tmp << EOF
  1.73043144E-01  6.59943338E-05
  1.73043144E-01 -6.59943338E-05
EOF
./test_values.sh inversion.tmp answer.tmp 2 "lsqr, minimize shear traction on 2 strike-slip faults" || exit 1

#echo ----------
#echo Finished Test \#7
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


echo --------------------------------------------------------------
echo Side Test: Are integer and double precision annealing algorithms the same?
echo ----------
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


echo --------------------------------------------------------------
echo Test \#8: 1 strike-slip fault, 4 3-component displacements, simulated annealing
echo ----------
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
$BIN_DIR/o92util -flt o92_flt.tmp -sta o92_sta.tmp -disp o92_disp.tmp -xy

# Prepare displacement observation and fault geometry files for fltinv
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$6}' o92_disp.tmp > fltinv_disp.tmp
awk '{print $1*1e3,$2*1e3,$3*1e3,$4,$5,$8*1e3,$9*1e3}' o92_flt.tmp > fltinv_flt.tmp

# Actual solution
#echo "Actual solution (strike-slip dip-slip):"
#awk '{printf("%16.8e%16.8e\n"),$7*cos($6*0.0174533),$7*sin($6*0.0174533)}' o92_flt.tmp

# Slip and rake constraints
echo 0 10 > slip.tmp
echo -30 60 > rake.tmp
echo 0.25 1 > step.tmp

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
    -o inversion.tmp
echo 1.18729953E+00 -1.69601249E-01 > answer.tmp
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: anneal, 1 fault, three-component disp" || exit 1

#echo ----------
#echo Finished Test \#8
#echo ----------
#echo




echo ----------------------------------------------------------
echo Test \#9: 4 dip-slip faults, 25 3-component displacements, covariance, annealing
echo ----------
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
$BIN_DIR/grid -x -2 6 -nx 5 -y -7 7 -ny 5 -z 0.0 -o o92_sta.tmp

# Calculate "observed" displacements with small noise
$BIN_DIR/o92util -flt o92_flt.tmp -sta o92_sta.tmp -disp o92_disp.tmp -xy
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
    -o inversion.tmp
cat > answer.tmp << EOF
  8.66128699E-01  1.81351839E+00
 -3.67871230E-01  1.45673734E+00
  2.47535983E-01  1.45810405E+00
 -3.11359215E-01  1.03791256E+00
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: anneal, 4 dip-slip faults, 25 observations, covariance matrix, annealing" || exit 1

#echo ----------
#echo Finished Test \#9
#echo ----------
#echo



echo ----------------------------------------------------------
echo Test \#10: 9 strike-slip faults, 25 three-component displacements, simulated annealing with pseudo-coupling
echo ----------
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
    -o inversion.tmp
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
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, pre-stresses, fix central fault slip, 8 surrounding unlocked faults" || exit 1


# Calculate displacements at grid of stations around faults
paste o92_flt.tmp inversion.tmp |\
    awk '{print $1,$2,$3,$4,$5,atan2($11,$10)/0.0174533,sqrt($10*$10+$11*$11),$8,$9}' > o92_flt_psc.tmp
$BIN_DIR/grid -x -4.5 5.2 -nx 5 -y -3.6 2.9 -ny 5 -z 0.0 -o o92_sta_psc.tmp
$BIN_DIR/o92util -flt o92_flt_psc.tmp -sta o92_sta_psc.tmp -haf haf.tmp -disp o92_disp_psc.tmp -xy
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
    -o inversion.tmp
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
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: simulated annealing + pseudo-coupling" || exit 1

grep Iteration anneal.log | awk '{print $6}' | head -20 > fit.tmp
$BIN_DIR/anneal_post -f anneal.log -obj obj.tmp
awk '{print $2}' obj.tmp | head -20 > fit.tmp
cat > answer.tmp << EOF
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
./test_values.sh fit.tmp answer.tmp 1 "fltinv: simulated annealing + pseudo-coupling, first 20 fits" || exit 1

#echo ----------
#echo Finished Test \#10
#echo ----------
#echo


echo ----------------------------------------------------------
echo Test \#11: 4 triangular strike-slip faults, 9 3-component displacements
echo ----------
# Input fault slip
# x1 y1 z1 x2 y2 z2 x3 y3 z3 ss ds ts
cat > tri_flt.tmp << EOF
76.3 35.8 2.0 76.4 35.7 5.0 72.7 32.5 2.0 1.8 0.6 0.0
76.9 35.3 6.0 76.4 35.7 5.0 72.7 32.5 2.0 1.3 0.4 0.0
76.9 35.3 6.0 72.4 29.1 4.8 72.7 32.5 2.0 0.7 1.1 0.0
71.0 29.9 1.5 72.4 29.1 4.8 72.7 32.5 2.0 0.2 1.4 0.0
EOF

# Station locations
$BIN_DIR/grid -x 70.0 80.0 -nx 3 -y 30.0 40.0 -ny 3 -z 0.0 -o tri_sta.tmp

# Calculate "observed" displacements with no noise
$BIN_DIR/triutil -flt tri_flt.tmp -sta tri_sta.tmp -disp tri_disp.tmp -xy

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
    -o inversion.tmp
cat > answer.tmp << EOF
  1.80000000E+00  6.00000000E-01
  1.30000000E+00  4.00000000E-01
  7.00000000E-01  1.10000000E+00
  2.00000000E-01  1.40000000E+00
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, 4 faults, three-component disp" || exit 1

#echo ----------
#echo Finished Test \#11
#echo ----------
#echo




echo ----------------------------------------------------------
echo Test \#12: 1 triangular dip-slip fault, pre-stresses
echo ----------
# One triangular fault
#
#   *(0,10,0)
#   |   --__
#   |       *(15,0,5)
#   |   --
#   *(0,-10,0)
#
echo 0 -10 0 15 0 5 0 10 0 > tri.tmp

# Triangle center
awk '{print ($1+$4+$7)/3,($2+$5+$8)/3,($3+$6+$9)/3}' tri.tmp > center.tmp

# One meter normal slip
echo 0 -1 0 > slip.tmp

# Compute stress at triangle center
paste tri.tmp slip.tmp > triutil_flt.tmp
cp center.tmp triutil_sta.tmp
echo shearmod 40e9 lame 40e9 > haf.tmp
$BIN_DIR/triutil \
    -flt triutil_flt.tmp \
    -sta triutil_sta.tmp \
    -strain triutil_stn.tmp \
    -stress triutil_sts.tmp \
    -xy \
    -haf haf.tmp

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
    -o inversion.tmp
echo   7.97433674E-17  9.99999930E-01 > answer.tmp
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: triangular dip-slip fault with pre-stresses" || exit 1

#echo ----------
#echo Finished Test \#12
#echo ----------
#echo



echo ----------------------------------------------------------
echo Test \#13: 1 triangular dip-slip fault, pre-stresses
echo ----------
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
echo 0 -10 0 15 0 5 0 10 0 > tri1.tmp
echo 0 -10 0 15 0 5 15 -20 5 > tri2.tmp

# Triangle centers
awk '{print ($1+$4+$7)/3,($2+$5+$8)/3,($3+$6+$9)/3}' tri1.tmp > center1.tmp
awk '{print ($1+$4+$7)/3,($2+$5+$8)/3,($3+$6+$9)/3}' tri2.tmp > center2.tmp

# One meter normal slip on locked fault
echo 0 -1 0 > slip1.tmp

# Compute stress at triangle 2 center
paste tri1.tmp slip1.tmp > triutil_flt1.tmp
echo lame 40e9 shear 40e9 > haf.tmp
cp center2.tmp triutil_sta2.tmp
$BIN_DIR/triutil \
    -flt triutil_flt1.tmp \
    -sta triutil_sta2.tmp \
    -strain triutil_stn2.tmp \
    -stress triutil_sts2.tmp \
    -haf haf.tmp \
    -xy

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
    -o inversion.tmp
echo -4.57250048E-02 -1.53447831E-01 > answer.tmp
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: lsqr, minimize shear traction on 1 other triangular fault" || exit 1

#echo ----------
#echo Finished Test \#13
#echo ----------
#echo


echo ----------------------------------------------------------
echo Test \#14: 1 triangular dip-slip fault, pre-stresses on itself, geographic coordinates
echo ----------
# One triangular fault, using geographic coordinates
#
#       *(3.00, 42.00, 2.0)
#      /    --__
#     /         *(3.05, 41.92, 5.0)
#    /      --
#    *(2.97, 41.85, 1.8)
#
echo 3.00 42.00 2.0 > pt1.tmp
echo 3.05 41.92 5.0 > pt2.tmp
echo 2.97 41.85 1.8 > pt3.tmp
paste pt1.tmp pt2.tmp pt3.tmp > tri.tmp

# Triangle center
awk '{printf("%.8f %.8f %.8f\n"),($1+$4+$7)/3,($2+$5+$8)/3,($3+$6+$9)/3}' tri.tmp > center.tmp

# One meter normal slip
echo 0 -1 0 > slip.tmp

# Compute stress at triangle center
paste tri.tmp slip.tmp > triutil_flt.tmp
cp center.tmp triutil_sta.tmp
echo lame 40e9 shear 40e9 > haf.tmp
$BIN_DIR/triutil \
    -flt triutil_flt.tmp \
    -sta triutil_sta.tmp \
    -strain triutil_stn.tmp \
    -stress triutil_sts.tmp \
    -haf haf.tmp \
    -geo

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
    -o inversion.tmp
echo 5.06302972E-09  1.00000001E+00 > answer.tmp
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: triangular fault with self-prestress, geographic coordinates" || exit 1

rm *.tmp
#echo ----------
#echo Finished Test \#14
#echo ----------
#echo




echo ----------------------------------------------------------
echo Test \#15: 1 triangular dip-slip fault, pre-stresses on adjacent fault, geographic coordinates
echo ----------
# Two triangular faults, using geographic coordinates
#
# -73.0 -72.8
#   1*---* -31.7
#    |1 /|3
#    | / |
#    |/ 2|4
#   2*---* -32.0
#
echo -73.0 -31.7 6.0 > pt1.tmp
echo -73.0 -32.0 6.0 > pt2.tmp
echo -72.8 -31.7 8.3 > pt3.tmp
echo -72.8 -32.0 8.3 > pt4.tmp
paste pt1.tmp pt2.tmp pt3.tmp > tri1.tmp
paste pt2.tmp pt3.tmp pt4.tmp > tri2.tmp

# One meter normal slip
echo 0 -1 0 > slip.tmp
paste tri1.tmp slip.tmp > triutil_flt.tmp

# Compute stress at adjacent triangle center
awk '{printf("%.8f %.8f %.8f\n"),($1+$4+$7)/3,($2+$5+$8)/3,($3+$6+$9)/3}' tri2.tmp > center.tmp
cp center.tmp triutil_sta.tmp
echo lame 40e9 shear 40e9 > haf.tmp
$BIN_DIR/triutil \
    -flt triutil_flt.tmp \
    -sta triutil_sta.tmp \
    -strain triutil_stn.tmp \
    -stress triutil_sts.tmp \
    -haf haf.tmp \
    -geo

# Prepare fltinv files
awk '{print $1,$2,$3*1e3,$4,$5,$6*1e3,$7,$8,$9*1e3}' tri2.tmp > fltinv_flt.tmp
awk '{print $4,$5,$6,$7,$8,$9}' triutil_sts.tmp > fltinv_sts.tmp
$BIN_DIR/fltinv \
    -mode lsqr \
    -flt fltinv_flt.tmp \
    -gf:model triangle \
    -prests fltinv_sts.tmp \
    -geo \
    -haf haf.tmp \
    -o inversion.tmp
# echo inversion with pre-stresses
# cat inversion.tmp

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
    -o inversion.tmp

cat > answer.tmp << EOF
  0.00000000E+00 -1.00000000E+00
 -2.89142383E-02 -2.93023235E-01
EOF
./test_values.sh inversion.tmp answer.tmp 2 "fltinv: 2 triangular faults, one locked, geographic coordinates" || exit 1

rm *.tmp
#echo ----------
#echo Finished Test \#15
#echo ----------
#echo




echo ----------------------------------------------------------
echo Test \#16: Euler pole, least-squares
echo ----------
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
$BIN_DIR/platemotion -f coords.tmp -pole -101/24/0.42 |\
    awk '{print $1,$2,0,$3/1e3,$4/1e3,0}' > vel.tmp

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
    -euler euler.tmp euler_out.tmp

awk '{print $1}' euler_out.tmp > inversion.tmp
echo -101 > answer.tmp
./test_values.sh inversion.tmp answer.tmp 1 "fltinv: Euler pole lon" -zero 100 || exit 1
awk '{print $2}' euler_out.tmp > inversion.tmp
echo 24 > answer.tmp
./test_values.sh inversion.tmp answer.tmp 1 "fltinv: Euler pole lat" -zero 100 || exit 1
awk '{print $3}' euler_out.tmp > inversion.tmp
echo 0.42 > answer.tmp
./test_values.sh inversion.tmp answer.tmp 1 "fltinv: Euler pole vel" -zero 1.0 || exit 1

rm *.tmp
#echo ----------
#echo Finished Test \#16
#echo ----------
#echo


echo ----------------------------------------------------------
echo Test \#17: Annealing with pseudo-coupling, plus an Euler pole
echo ----------
# Stations near North Island, New Zealand
grid -x 173 177 -dx 1 -y -42 -37 -dy 1 -z 0 > sta.tmp

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
$BIN_DIR/fltinv \
    -mode lsqr \
    -geo \
    -flt fltinv_tri.tmp \
    -flt:slip fltinv_slip.tmp \
    -haf haf.tmp \
    -prests fltinv_sts.tmp \
    -gf:model triangle \
    -lsqr:mode gesv \
    -o inversion.tmp

# Calculate locking generated velocity
paste tri.tmp inversion.tmp > tri_slip.tmp
$BIN_DIR/triutil -flt tri_slip.tmp -sta sta.tmp -haf haf.tmp -disp locking_vel.tmp

# Set Euler pole far north of points to give them E velocity
platemotion -f sta.tmp -pole 175.0/3.0/0.13 > euler_vel.tmp

# Superimpose Euler velocity and locking signal
paste locking_vel.tmp euler_vel.tmp |\
    awk '{print $1,$2,$3*1e3,$4+$9,$5+$10,$6}' > fltinv_vel.tmp

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
    -o inversion.tmp
$BIN_DIR/anneal_post -f anneal.log -obj obj.tmp

rm *.tmp
#echo ----------
#echo Finished Test \#17
#echo ----------
#echo



#####
#	CLEAN UP
#####
rm anneal.log

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
WID=4
WID2=2
LEN=6
grid -x 2.1 18.1 -dx 4 -y 3 27 -dy 6 |\
    awk '{print $1*cos(30*0.0174533),$2,$1*sin(30*0.0174533)}' |\
    awk '{
        print $1+2*cos(30*0.0174533),$2+3,$3+2*cos(30*0.0174533)
        print $1+2*cos(30*0.0174533),$2-3,$3+2*cos(30*0.0174533)
        print $1-2*cos(30*0.0174533),$2-3,$3-2*cos(30*0.0174533)
        print $1-2*cos(30*0.0174533),$2+3,$3-2*cos(30*0.0174533)
    }'
exit
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
#./test_values.sh inversion.tmp answer.tmp 2 "fltinv: triangular fault with self-prestress, geographic coordinates" || exit 1

#echo ----------
#echo Finished Test \#16
#echo ----------
#echo

#####
#	CLEAN UP
#####
rm *.tmp
