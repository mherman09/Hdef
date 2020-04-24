#!/bin/bash

#####
#	SET PATH TO HDEF EXECUTABLE
#####
# Check if mtutil is set in PATH
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which mtutil | xargs dirname)
fi

# Check for mtutil in same directory as script
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/mtutil | xargs dirname)
fi

# Check for mtutil in relative directory ../bin (assumes script is in Hdef/dir)
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/../bin/mtutil | xargs dirname)
fi

# Check for mtutil in relative directory ../build (assumes script is in Hdef/dir)
if [ "$BIN_DIR" == "" ]
then
    BIN_DIR=$(which $(dirname $0)/../build/mtutil | xargs dirname)
fi

# Hdef executables are required!
if [ "$BIN_DIR" == "" ]
then
    echo "$0: unable to find Hdef executable mtutil; exiting" 1>&2
    exit 1
fi


#####
#	SET PATH TO TEST_VALUES SCRIPT
#####
TEST_BIN_DIR=`echo $0 | xargs dirname`


#####
#	RUN TEST
#####
trap "rm -f *.tmp" 0 1 2 3 8 9


cat > input.tmp << EOF
7.0
8.0
EOF
#$BIN_DIR/mtutil -mag input.tmp -mom
#cat input.tmp | $BIN_DIR/mtutil -mag -mom
#$BIN_DIR/mtutil -mag 7.0 -mom

# Magnitude to seismic moment
echo 1.9952623149688665E+020 > answer.tmp
$BIN_DIR/mtutil -mag 7.5 -mom mom.tmp
$TEST_BIN_DIR/test_values.sh mom.tmp answer.tmp 1 "mtutil: mag2mom" || exit 1

# Seismic moment to magnitude
echo 8.6333333333333329 > answer.tmp
$BIN_DIR/mtutil -mom 1.0e22 -mag mag.tmp
$TEST_BIN_DIR/test_values.sh mag.tmp answer.tmp 1 "mtutil: mom2mag" || exit 1

# Moment tensor to double couple percentage
echo 0.58784927087317074 > answer.tmp
$BIN_DIR/mtutil -mij 4.173e+17,7.304e+17,-1.1477e+18,-3.5502e+18,-3.843e+18,1.1714e+18 -dcp dcp.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp dcp.tmp 1 "mtutil: mij2dcp" || exit 1

# Moment tensor to moment magnitude
echo 6.4522434199195224 > answer.tmp
$BIN_DIR/mtutil -mij 4.173e+17,7.304e+17,-1.1477e+18,-3.5502e+18,-3.843e+18,1.1714e+18 -mag mag.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp mag.tmp 1 "mtutil: mij2mag" || exit 1

# Moment tensor to seismic moment
echo 5.3501397941634130E+018 > answer.tmp
$BIN_DIR/mtutil -mij 4.173e+17,7.304e+17,-1.1477e+18,-3.5502e+18,-3.843e+18,1.1714e+18 -mom mom.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp mom.tmp 1 "mtutil: mij2mom" || exit 1

# Moment tensor to PNT axes
echo 0.55173391074654210       0.61662524753526782       0.56157189729757373       0.76332257667899917      -0.64464366270352003       -4.2109287198004358E-002  0.33604811510326971       0.45189362934020239      -0.82635574185533400       -5.3220383822829574E+021   4.7233928809523143E+019   5.2738044534734329E+021 > answer.tmp
$BIN_DIR/mtutil -mij 1.923e21,-9.270e20,-9.970e20,3.811e21,-3.115e21,1.033d21 -pnt pnt.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp pnt.tmp 12 "mtutil: mij2pnt" || exit 1

# Moment tensor to strike, dip, rake
echo    213.51387373513245        40.404418834012674        65.379290596519624        64.553919649886069        53.896206591202564        109.52541253455682 > answer.tmp
$BIN_DIR/mtutil -mij 1.3892d18,-8.1410d17,-5.7510d17,5.3130d17,-7.6200d16,-8.2620d17 -sdr sdr.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp sdr.tmp 6 "mtutil: mij2sdr" || exit 1

# Moment tensor to slip vectors
echo   0.34714538840760611      -0.72957158138707823       0.58924985103068972      -0.54042032771493020       0.35788462890633410       0.76148832018952395 > answer.tmp
$BIN_DIR/mtutil -mij 1.3892d18,-8.1410d17,-5.7510d17,5.3130d17,-7.6200d16,-8.2620d17 -sv sv.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp sv.tmp 6 "mtutil mij2sv" || exit 1

# Moment tensor to ternary representation
echo   0.78904401692544923        2.4806642679335161E-003  0.20847531880661724 > answer.tmp
$BIN_DIR/mtutil -mij 1.4336e+18,-1.3870e+17,-1.2949e+18,4.7140e+17,-1.8573e+18,1.5080e+17 -ternary ternary.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp ternary.tmp 3 "mtutil: mij2ternary" || exit 1

# PNT axes to double couple percentage
echo   0.58784927087317074 > answer.tmp
$BIN_DIR/mtutil -pnt -0.64401587537772065,0.31483723501120958,-0.69722382899034951,0.60972258020181702,0.76168884485343713,-0.21924524810785431,0.46204104520432626,-0.56631053237277174,-0.68250307945837219,-4.7355476386355784d+018,-1.2291843110556710d+018,5.9647319496912486d+018 -dcp dcp.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp dcp.tmp 1 "mtutil: pnt2dcp" || exit 1

# PNT axes to moment magnitude
echo 6.4522434199195224 > answer.tmp
$BIN_DIR/mtutil -pnt -0.64401587537772065,0.31483723501120958,-0.69722382899034951,0.60972258020181702,0.76168884485343713,-0.21924524810785431,0.46204104520432626,-0.56631053237277174,-0.68250307945837219,-4.7355476386355784d+018,-1.2291843110556710d+018,5.9647319496912486d+018 -mag mag.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp mag.tmp 1 "mtutil: pnt2mag" || exit 1

# PNT axes to moment tensor
echo 1.9229999999999987E+021  -9.2700000000000092E+020  -9.9700000000000039E+020   3.8109999999999997E+021  -3.1150000000000003E+021   1.0330000000000004E+021 > answer.tmp
$BIN_DIR/mtutil -pnt 0.55173391074654210d0,0.61662524753526782d0,0.56157189729757373d0,0.76332257667899917d0,-0.64464366270352003d0,-4.2109287198004358d-002,0.33604811510326971d0,0.45189362934020239d0,-0.82635574185533400d0,-5.3220383822829574d+021,4.7233928809523143d+019,5.2738044534734329d+021 -mij mij.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp mij.tmp 6 "mtutil: pnt2mij" || exit 1

# PNT axes to seismic moment
echo 5.3501397941634130E+018 > answer.tmp
$BIN_DIR/mtutil -pnt -0.64401587537772065,0.31483723501120958,-0.69722382899034951,0.60972258020181702,0.76168884485343713,-0.21924524810785431,0.46204104520432626,-0.56631053237277174,-0.68250307945837219,-4.7355476386355784d+018,-1.2291843110556710d+018,5.9647319496912486d+018 -mom mom.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp mom.tmp 1 "mtutil: pnt2mom" || exit 1

# PNT axes to strike, dip, and rake
echo 213.51387373513245        40.404418834012674        65.379290596519624        64.553919649886069        53.896206591202564        109.52541253455682 > answer.tmp
$BIN_DIR/mtutil -pnt 0.62760373661893987d0,-0.76894766054189612d0,-0.12179098952340160d0,0.76644370223769165d0,0.58278975628002727d0,0.27003731459790342d0,-0.13666602021762725d0,-0.26282236457769625d0,0.95511612047732908d0,-1.5725830409879032d18,2.6280182158839164d16,1.5463028588290642d18 -sdr sdr.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp sdr.tmp 6 "mtutil: pnt2sdr" || exit 1

# PNT axes to slip vectors
echo 0.34714538840760611      -0.72957158138707823       0.58924985103068972      -0.54042032771493020       0.35788462890633410       0.76148832018952395 > answer.tmp
$BIN_DIR/mtutil -mij 1.3892d18,-8.1410d17,-5.7510d17,5.3130d17,-7.6200d16,-8.2620d17 -pnt | $BIN_DIR/mtutil -pnt -sv sv.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp sv.tmp 6 "mtutil: pnt2sv" || exit 1

# PNT axes to ternary representation
echo 0.91224680359566357        7.2920151275247064E-002   1.4833045129089315E-002 > answer.tmp
$BIN_DIR/mtutil -pnt 0.62760373661893987d0,-0.76894766054189612d0,-0.12179098952340160d0,0.76644370223769165d0,0.58278975628002727d0,0.27003731459790342d0,-0.13666602021762725d0,-0.26282236457769625d0,0.95511612047732908d0,-1.5725830409879032d18,2.6280182158839164d16,1.5463028588290642d18 -ternary ternary.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp ternary.tmp 3 "mtutil: pnt2ter" || exit 1

# Strike, dip, rake to moment tensor
echo 0.22003533408899165       0.18410177651383547      -0.40413711060282714       0.74058859778607167      -0.16465065799140619       0.54918192012069555 > answer.tmp
$BIN_DIR/mtutil -sdr 9,39,167 -mij mij.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp mij.tmp 6 "mtutil: sdr2mij" || exit 1

# Strike, dip, rake to PNT axes
echo 0.61670145162508538      -0.77643201628987624      -0.12974067844569803       0.77519713030606596       0.57032765296921284       0.27165378227418457      -0.13692599727134791      -0.26810388348300168       0.95360976239370576       -1.0000000000000002        1.6132928326584306E-016   1.0000000000000000 > answer.tmp
$BIN_DIR/mtutil -sdr 214,40,65 -pnt pnt.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp pnt.tmp 12 "mtutil: sdr2pnt" || exit 1

# Strike, dip, rake to second nodal plane
echo 83.998199510808121        46.432414978165248       -58.179239464578913 > answer.tmp
$BIN_DIR/mtutil -sdr 222,52,-119 -sdr sdr.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp sdr.tmp 3 "mtutil: sdr2sdr" || exit 1

# Strike, dip, rake to slip vector
echo 0.29134327529638415      -0.77164718283107014       0.56540226491273315 > answer.tmp
$BIN_DIR/mtutil -sdr 292,43,124 -sv sv.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp sv.tmp 3 "mtutil: sdr2sv" || exit 1

# Strike, dip, rake to ternary representation
echo 0.78910589852855950        2.3368434059280944E-003  0.20855725806551242 > answer.tmp
$BIN_DIR/mtutil -sdr 357,18,99 -ternary ternary.tmp
$TEST_BIN_DIR/test_values.sh answer.tmp ternary.tmp 3 "mtutil: sdr2ter" || exit 1
