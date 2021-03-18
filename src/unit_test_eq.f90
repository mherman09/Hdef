program main

use test, only: test_value
use io, only: stdout

use eq

implicit none

double precision :: str1, dip1, rak1, str2, dip2, rak2
double precision :: mrr, mpp, mtt, mrp, mrt, mtp
double precision :: pnt(12)
double precision :: slip_vec(3), slip_vec2(3)
double precision :: fth, fss, fno
double precision :: dcp
double precision :: mom, mag
double precision :: wid, len

! mij2dcp
! Moment tensor to double couple percentage
mrr = 4.173d+17
mtt = 7.304d+17
mpp = -1.1477d+18
mrt = -3.5502d+18
mrp = -3.843d+18
mtp = 1.1714d+18
call mij2dcp(mrr,mtt,mpp,mrt,mrp,mtp,dcp)
call test_value(dcp,0.58784927087317074d0,'mij2dcp: dcp')

! mij2mag
call mij2mag(mrr,mtt,mpp,mrt,mrp,mtp,mag,'dc')
call test_value(mag,6.4522434199195224d0,'mij2mag: mag (dc)')
call mij2mag(mrr,mtt,mpp,mrt,mrp,mtp,mag,'nondc')
call test_value(mag,6.4578638576982641d0,'mij2mag: mag (non-dc)')

! mij2mom
! Moment tensor to scalar moment
mrr = 4.173d+17
mtt = 7.304d+17
mpp = -1.1477d+18
mrt = -3.5502d+18
mrp = -3.843d+18
mtp = 1.1714d+18
call mij2mom(mrr,mtt,mpp,mrt,mrp,mtp,mom,'dc')
call test_value(mom,5.3501397941634130d+018,'mij2mom: mom (dc)')
call mij2mom(mrr,mtt,mpp,mrt,mrp,mtp,mom,'nondc')
call test_value(mom,5.4550129578214554d+018,'mij2mom: mom (non-dc)')

! mij2pnt
! Moment tensor to P, N, and T axes
mrr =  1.923d21
mtt = -9.270d20
mpp = -9.970d20
mrt =  3.811d21
mrp = -3.115d21
mtp =  1.033d21
call mij2pnt(mrr,mtt,mpp,mrt,mrp,mtp,pnt)
call test_value(abs(pnt(1)),  0.55173391074654210d0,  'mij2pnt: P axis e-component')
call test_value(abs(pnt(2)),  0.61662524753526782d0,  'mij2pnt: P axis n-component')
call test_value(abs(pnt(3)),  0.56157189729757373d0,  'mij2pnt: P axis z-component')
call test_value(abs(pnt(4)),  0.76332257667899917d0,  'mij2pnt: N axis e-component')
call test_value(abs(pnt(5)),  0.64464366270352003d0,  'mij2pnt: N axis n-component')
call test_value(abs(pnt(6)),  4.2109287198004358d-2,  'mij2pnt: N axis z-component')
call test_value(abs(pnt(7)),  0.33604811510326971d0,  'mij2pnt: T axis e-component')
call test_value(abs(pnt(8)),  0.45189362934020239d0,  'mij2pnt: T axis n-component')
call test_value(abs(pnt(9)),  0.82635574185533400d0,  'mij2pnt: T axis z-component')
call test_value(pnt(10),-5.3220383822829574d+021,'mij2pnt: P magnitude')
call test_value(pnt(11), 4.7233928809523143d+019,'mij2pnt: N magnitude')
call test_value(pnt(12), 5.2738044534734329d+021,'mij2pnt: T magnitude')

! mij2sdr
! Moment tensor to best-fitting double couple strike, dip, and rake
mrr =  1.3892d18
mtt = -8.1410d17
mpp = -5.7510d17
mrt =  5.3130d17
mrp = -7.6200d16
mtp = -8.2620d17
call mij2sdr(mrr,mtt,mpp,mrt,mrp,mtp,str1,dip1,rak1,str2,dip2,rak2)
if (str1.gt.200.0d0) then
    call test_value(str1,213.51387373513245d0,'mij2sdr: str1')
    call test_value(dip1,40.404418834012674d0,'mij2sdr: dip1')
    call test_value(rak1,65.379290596519624d0,'mij2sdr: rak1')
    call test_value(str2,64.553919649886069d0,'mij2sdr: str2')
    call test_value(dip2,53.896206591202564d0,'mij2sdr: dip2')
    call test_value(rak2,109.52541253455682d0,'mij2sdr: rak2')
else
    call test_value(str2,213.51387373513245d0,'mij2sdr: str2')
    call test_value(dip2,40.404418834012674d0,'mij2sdr: dip2')
    call test_value(rak2,65.379290596519624d0,'mij2sdr: rak2')
    call test_value(str1,64.553919649886069d0,'mij2sdr: str1')
    call test_value(dip1,53.896206591202564d0,'mij2sdr: dip1')
    call test_value(rak1,109.52541253455682d0,'mij2sdr: rak1')
endif

! mij2sv
! Moment tensor to slip vectors of best-fitting double couple source
mrr =  1.3892d18
mtt = -8.1410d17
mpp = -5.7510d17
mrt =  5.3130d17
mrp = -7.6200d16
mtp = -8.2620d17
call mij2sv(mrr,mtt,mpp,mrt,mrp,mtp,slip_vec,slip_vec2)
if (slip_vec(1).gt.0.0d0) then
    call test_value(slip_vec(1), 0.34714538840760611d0,'mij2sv: slip_vec_1(1)')
    call test_value(slip_vec(2),-0.72957158138707823d0,'mij2sv: slip_vec_1(2)')
    call test_value(slip_vec(3), 0.58924985103068972d0,'mij2sv: slip_vec_1(3)')
    call test_value(slip_vec2(1),-0.54042032771493020d0,'mij2sv: slip_vec_2(1)')
    call test_value(slip_vec2(2), 0.35788462890633410d0,'mij2sv: slip_vec_2(2)')
    call test_value(slip_vec2(3), 0.76148832018952395d0,'mij2sv: slip_vec_2(3)')
else
    call test_value(slip_vec2(1), 0.34714538840760611d0,'mij2sv: slip_vec_2(1)')
    call test_value(slip_vec2(2),-0.72957158138707823d0,'mij2sv: slip_vec_2(2)')
    call test_value(slip_vec2(3), 0.58924985103068972d0,'mij2sv: slip_vec_2(3)')
    call test_value(slip_vec(1),-0.54042032771493020d0,'mij2sv: slip_vec_1(1)')
    call test_value(slip_vec(2), 0.35788462890633410d0,'mij2sv: slip_vec_1(2)')
    call test_value(slip_vec(3), 0.76148832018952395d0,'mij2sv: slip_vec_1(3)')
endif

! mij2ter
! Moment tensor to fraction thrust, strike-slip, and normal
mrr =  1.4336e+18
mtt = -1.3870e+17
mpp = -1.2949e+18
mrt =  4.7140e+17
mrp = -1.8573e+18
mtp =  1.5080e+17
call mij2ter(mrr,mtt,mpp,mrt,mrp,mtp,fth,fss,fno)
call test_value(fth,0.78904402299330922d0,'mij2ter: fth')
call test_value(fss,2.4806651462488932d-003,'mij2ter: fss')
call test_value(fno,0.20847531186044141d0,'mij2ter: fno')

! pnt2dcp
! P, N, and T axes to double couple percentage
! mrr = 4.173d+17
! mtt = 7.304d+17
! mpp = -1.1477d+18
! mrt = -3.5502d+18
! mrp = -3.843d+18
! mtp = 1.1714d+18
! call mij2pnt(mrr,mtt,mpp,mrt,mrp,mtp,pnt)
pnt(1) = -0.64401587537772065
pnt(2) =  0.31483723501120958
pnt(3) = -0.69722382899034951
pnt(4) =  0.60972258020181702
pnt(5) =  0.76168884485343713
pnt(6) = -0.21924524810785431
pnt(7) =  0.46204104520432626
pnt(8) = -0.56631053237277174
pnt(9) = -0.68250307945837219
pnt(10) = -4.7355476386355784d+018
pnt(11) = -1.2291843110556710d+018
pnt(12) =  5.9647319496912486d+018
call pnt2dcp(pnt,dcp)
call test_value(dcp,0.58784927087317074d0,'pnt2dcp: dcp')

! pnt2mij
! P, N, and T axes to moment tensor
pnt(1) =   0.55173391074654210d0
pnt(2) =   0.61662524753526782d0
pnt(3) =   0.56157189729757373d0
pnt(4) =   0.76332257667899917d0
pnt(5) =  -0.64464366270352003d0
pnt(6) =   -4.2109287198004358d-002
pnt(7) =   0.33604811510326971d0
pnt(8) =   0.45189362934020239d0
pnt(9) =  -0.82635574185533400d0
pnt(10) =   -5.3220383822829574d+021
pnt(11) =    4.7233928809523143d+019
pnt(12) =    5.2738044534734329d+021
call pnt2mij(pnt,mrr,mtt,mpp,mrt,mrp,mtp)
call test_value(mrr, 1.9229999999999987d+021,'pnt2mij: mrr')
call test_value(mtt,-9.2700000000000092d+020,'pnt2mij: mtt')
call test_value(mpp,-9.9700000000000039d+020,'pnt2mij: mpp')
call test_value(mrt, 3.8109999999999997d+021,'pnt2mij: mrt')
call test_value(mrp,-3.1150000000000003d+021,'pnt2mij: mrp')
call test_value(mtp, 1.0330000000000004d+021,'pnt2mij: mtp')

! pnt2mag
call pnt2mag(pnt,mag)
call test_value(mag,8.4494036748523378d0,'pnt2mag: mag')
call pnt2mag_nondc(pnt,mag)
call test_value(mag,8.4494124278966609d0,'pnt2mag_nondc: mag')

! pnt2mom
! P, N, and T axes to scalar moment
pnt(1) = -0.64401587537772065
pnt(2) =  0.31483723501120958
pnt(3) = -0.69722382899034951
pnt(4) =  0.60972258020181702
pnt(5) =  0.76168884485343713
pnt(6) = -0.21924524810785431
pnt(7) =  0.46204104520432626
pnt(8) = -0.56631053237277174
pnt(9) = -0.68250307945837219
pnt(10) = -4.7355476386355784d+018
pnt(11) = -1.2291843110556710d+018
pnt(12) =  5.9647319496912486d+018
call pnt2mom(pnt,mom)
call test_value(mom,5.3501397941634130d+018,'pnt2mom: mom')
call pnt2mom_nondc(pnt,mom)
call test_value(mom,5.4550129578214554d+018,'pnt2mom_nondc: mom')

! pnt2sdr
! P, N, and T axes to best-fitting double couple strike, dip, and rake
pnt(1) =  0.62760373661893987d0
pnt(2) = -0.76894766054189612d0
pnt(3) = -0.12179098952340160d0
pnt(4) =  0.76644370223769165d0
pnt(5) =  0.58278975628002727d0
pnt(6) =  0.27003731459790342d0
pnt(7) = -0.13666602021762725d0
pnt(8) = -0.26282236457769625d0
pnt(9) =  0.95511612047732908d0
pnt(10) = -1.5725830409879032d18
pnt(11) =  2.6280182158839164d16
pnt(12) =  1.5463028588290642d18
call pnt2sdr(pnt,str1,dip1,rak1,str2,dip2,rak2)
call test_value(str1,213.51387373513245d0,'pnt2sdr: str1')
call test_value(dip1,40.404418834012674d0,'pnt2sdr: dip1')
call test_value(rak1,65.379290596519624d0,'pnt2sdr: rak1')
call test_value(str2,64.553919649886069d0,'pnt2sdr: str2')
call test_value(dip2,53.896206591202564d0,'pnt2sdr: dip2')
call test_value(rak2,109.52541253455682d0,'pnt2sdr: rak2')

! pnt2sv
mrr =  1.3892d18
mtt = -8.1410d17
mpp = -5.7510d17
mrt =  5.3130d17
mrp = -7.6200d16
mtp = -8.2620d17
call mij2pnt(mrr,mtt,mpp,mrt,mrp,mtp,pnt)
call pnt2sv(pnt,slip_vec,slip_vec2)
if (slip_vec(1).gt.0.0d0) then
    call test_value(slip_vec(1), 0.34714538840760611d0,'pnt2sv: slip_vec_1(1)')
    call test_value(slip_vec(2),-0.72957158138707823d0,'pnt2sv: slip_vec_1(2)')
    call test_value(slip_vec(3), 0.58924985103068972d0,'pnt2sv: slip_vec_1(3)')
    call test_value(slip_vec2(1),-0.54042032771493020d0,'pnt2sv: slip_vec_2(1)')
    call test_value(slip_vec2(2), 0.35788462890633410d0,'pnt2sv: slip_vec_2(2)')
    call test_value(slip_vec2(3), 0.76148832018952395d0,'pnt2sv: slip_vec_2(3)')
else
    call test_value(slip_vec2(1), 0.34714538840760611d0,'pnt2sv: slip_vec_2(1)')
    call test_value(slip_vec2(2),-0.72957158138707823d0,'pnt2sv: slip_vec_2(2)')
    call test_value(slip_vec2(3), 0.58924985103068972d0,'pnt2sv: slip_vec_2(3)')
    call test_value(slip_vec(1),-0.54042032771493020d0,'pnt2sv: slip_vec_1(1)')
    call test_value(slip_vec(2), 0.35788462890633410d0,'pnt2sv: slip_vec_1(2)')
    call test_value(slip_vec(3), 0.76148832018952395d0,'pnt2sv: slip_vec_1(3)')
endif

! pnt2ter
! P, N, and T axes to fraction thrust, strike-slip, and normal
pnt(1) =  0.62760373661893987d0
pnt(2) = -0.76894766054189612d0
pnt(3) = -0.12179098952340160d0
pnt(4) =  0.76644370223769165d0
pnt(5) =  0.58278975628002727d0
pnt(6) =  0.27003731459790342d0
pnt(7) = -0.13666602021762725d0
pnt(8) = -0.26282236457769625d0
pnt(9) =  0.95511612047732908d0
pnt(10) = -1.5725830409879032d18
pnt(11) =  2.6280182158839164d16
pnt(12) =  1.5463028588290642d18
call pnt2ter(pnt,fth,fss,fno)
call test_value(fth,0.91224680359566357d0,'pnt2ter: fth')
call test_value(fss,7.2920151275247064d-002,'pnt2ter: fss')
call test_value(fno,1.4833045129089315d-002,'pnt2ter: fno')

! sdr2mij
! Strike, dip, and rake to moment tensor in spherical coordinates
str1 = 9.0d0
dip1 = 39.0d0
rak1 = 167.0d0
call sdr2mij(str1,dip1,rak1,mrr,mtt,mpp,mrt,mrp,mtp)
call test_value(mrr, 0.22003533408899165d0,'sdr2mij: mrr')
call test_value(mtt, 0.18410177651383547d0,'sdr2mij: mtt')
call test_value(mpp,-0.40413711060282714d0,'sdr2mij: mpp')
call test_value(mrt, 0.74058859778607167d0,'sdr2mij: mrt')
call test_value(mrp,-0.16465065799140619d0,'sdr2mij: mrp')
call test_value(mtp, 0.54918192012069555d0,'sdr2mij: mtp')

! sdr2pnt
! Strike, dip, and rake to P, N, and T axes
str1 = 214.0d0
dip1 = 40.0d0
rak1 = 65.0d0
call sdr2pnt(str1,dip1,rak1,pnt)
call test_value(abs(pnt(1)), 0.61670145162508538d0,'sdr2pnt: P axis e-component')
call test_value(abs(pnt(2)), 0.77643201628987624d0,'sdr2pnt: P axis n-component')
call test_value(abs(pnt(3)), 0.12974067844569803d0,'sdr2pnt: P axis z-component')
call test_value(abs(pnt(4)), 0.77519713030606596d0,'sdr2pnt: N axis e-component')
call test_value(abs(pnt(5)), 0.57032765296921284d0,'sdr2pnt: N axis n-component')
call test_value(abs(pnt(6)), 0.27165378227418457d0,'sdr2pnt: N axis z-component')
call test_value(abs(pnt(7)), 0.13692599727134791d0,'sdr2pnt: T axis e-component')
call test_value(abs(pnt(8)), 0.26810388348300168d0,'sdr2pnt: T axis n-component')
call test_value(abs(pnt(9)), 0.95360976239370576d0,'sdr2pnt: T axis z-component')
call test_value(pnt(10),-1.0d0,'sdr2pnt: P magnitude')
call test_value(pnt(11), 0.0d0,'sdr2pnt: N magnitude')
call test_value(pnt(12), 1.0d0,'sdr2pnt: T magnitude')

! sdr2sdr
! Strike, dip, and rake to strike, dip, and rake of second nodal plane
str1 = 222.0d0
dip1 = 52.0d0
rak1 = -119.0d0
call sdr2sdr(str1,dip1,rak1,str2,dip2,rak2)
call test_value(str2,83.998199510808121d0,'sdr2sdr: normal str')
call test_value(dip2,46.432414978165248d0,'sdr2sdr: normal dip')
call test_value(rak2,-58.179239464578913d0,'sdr2sdr: normal rak')
str1 = 109.0d0
dip1 = 82.0d0
rak1 = 52.0d0
call sdr2sdr(str1,dip1,rak1,str2,dip2,rak2)
call test_value(str2,8.8996492587401121d0,'sdr2sdr: oblique str')
call test_value(dip2,38.708111069705964d0,'sdr2sdr: oblique dip')
call test_value(rak2,167.14107738093375d0,'sdr2sdr: oblique rak')

! sdr2sv
! Strike, dip, and rake to slip vector
str1 = 292.0d0
dip1 = 43.0d0
rak1 = 124.0d0
call sdr2sv(str1,dip1,rak1,slip_vec)
call test_value(slip_vec(1), 0.29134327529638415d0,'sdr2sv: slip vector e-component')
call test_value(slip_vec(2),-0.77164718283107014d0,'sdr2sv: slip vector e-component')
call test_value(slip_vec(3), 0.56540226491273315d0,'sdr2sv: slip vector e-component')

! sdr2ter
! Strike, dip, and rake to fraction thrust, strike-slip, and normal
str1 = 357.0d0
dip1 = 18.0d0
rak1 = 99.0d0
call sdr2ter(str1,dip1,rak1,fth,fss,fno)
call test_value(fth,0.78910589852855950d0,'sdr2ter: fth')
call test_value(fss,2.3368434059280944d-3,'sdr2ter: fss')
call test_value(fno,0.20855725806551242d0,'sdr2ter: fno')

! mom2mag
mag = 7.0d0
call mag2mom(mag,mom)
call test_value(mom,3.5481338923357311d+019,'mag2mom')

! mag2mom
call mom2mag(mom,mag)
call test_value(mag,7.0d0,'mom2mag')

! empirical
mag = 7.1d0
! Wells and Coppersmith (1994)
call empirical(mag,wid,len,'WC','ss')
call test_value(wid,14.354894333536560d0,'empirical: WC-ss wid')
call test_value(len,67.920363261718506d0,'empirical: WC-ss len')
call empirical(mag,wid,len,'WC','th')
call test_value(wid,19.998618696327419d0,'empirical: WC-th wid')
call test_value(len,49.888448746001167d0,'empirical: WC-th len')
call empirical(mag,wid,len,'WC','no')
call test_value(wid,22.130947096056374d0,'empirical: WC-no wid')
call test_value(len,46.773514128719810d0,'empirical: WC-no len')
call empirical(mag,wid,len,'WC','ot')
call test_value(wid,18.281002161427416d0,'empirical: WC-ot wid')
call test_value(len,56.104797603246951d0,'empirical: WC-ot len')
call empirical(mag,wid,len,'WC:print','')
! Mai and Beroza (2000)
call empirical(mag,wid,len,'MB','ss')
call test_value(wid,14.757065332758916d0,'empirical: MB-ss wid')
call test_value(len,37.153522909717204d0,'empirical: MB-ss len')
call empirical(mag,wid,len,'MB','th')
call test_value(wid,24.266100950824043d0,'empirical: MB-th wid')
call test_value(len,30.902954325135855d0,'empirical: MB-th len')
call empirical(mag,wid,len,'MB','no')
call test_value(wid,24.266100950824043d0,'empirical: MB-no wid')
call test_value(len,30.902954325135855d0,'empirical: MB-no len')
call empirical(mag,wid,len,'MB','ot')
call test_value(wid,17.947336268325209d0,'empirical: MB-ot wid')
call test_value(len,35.727283815192813d0,'empirical: MB-ot len')
call empirical(mag,wid,len,'MB:print','')
! Blaser et al. (2010)
call empirical(mag,wid,len,'B','ss')
call test_value(wid,16.710906143107071d0,'empirical: B-ss wid')
call test_value(len,71.449632607551280d0,'empirical: B-ss len')
call empirical(mag,wid,len,'B','th')
call test_value(wid,25.468302525850405d0,'empirical: B-th wid')
call test_value(len,47.533522594280484d0,'empirical: B-th len')
call empirical(mag,wid,len,'B','no')
call test_value(wid,22.698648518838201d0,'empirical: B-no wid')
call test_value(len,60.534087475391331d0,'empirical: B-no len')
call empirical(mag,wid,len,'B','ot')
call test_value(wid,22.438819237827637d0,'empirical: B-ot wid')
call test_value(len,54.575786109127051d0,'empirical: B-ot len')
call empirical(mag,wid,len,'B:print','')
! Yen and Ma (2011)
call empirical(mag,wid,len,'YM','ss')
call test_value(wid,40.179081084893753d0,'empirical: YM-ss wid')
call test_value(len,54.954087385762257d0,'empirical: YM-ss len')
call empirical(mag,wid,len,'YM','th')
call test_value(wid,33.806483620598016d0,'empirical: YM-th wid')
call test_value(len,41.114972110451959d0,'empirical: YM-th len')
call empirical(mag,wid,len,'YM','no')
call test_value(wid,33.806483620598016d0,'empirical: YM-no wid')
call test_value(len,41.114972110451959d0,'empirical: YM-no len')
call empirical(mag,wid,len,'YM','ot')
call test_value(wid,38.018939632056046d0,'empirical: YM-ot wid')
call test_value(len,62.950618285719294d0,'empirical: YM-ot len')
call empirical(mag,wid,len,'YM:print','')
! Allen and Hayes (2017)
call empirical(mag,wid,len,'AH','ss')
call test_value(wid,12.445146117713850d0,'empirical: AH-ss wid')
call test_value(len,46.025657358135582d0,'empirical: AH-ss len')
call empirical(mag,wid,len,'AH','th')
call test_value(wid,31.477483141013156d0,'empirical: AH-th wid')
call test_value(len,37.411058827205331d0,'empirical: AH-th len')
call empirical(mag,wid,len,'AH','no')
call test_value(wid,20.183663636815606d0,'empirical: AH-no wid')
call test_value(len,40.086671762730262d0,'empirical: AH-no len')
call empirical(mag,wid,len,'AH','ot')
call test_value(wid,29.853826189179586d0,'empirical: AH-ot wid')
call test_value(len,27.733201046518410d0,'empirical: AH-ot len')
call empirical(mag,wid,len,'AH:print','')

write(stdout,*) 'eq_module unit test passed'
end
