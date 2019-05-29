program main

use test, only: test_value
use io, only: stdout
use trig, only: d2r

use okada92

implicit none

double precision :: disp(3), strain(3,3), sta_coord(3), evdp, dip, moment(4), lambda, shear_modulus
double precision :: slip_vec(3), wid, len
! double precision :: ux,uy,uz,xin,yin,stdp,dipin,rakin,slip,vp,vs,dens,area
double precision :: uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz

#ifdef UNIT_TEST
  lambda = 40.0d0
  shear_modulus = 40.0d0
  call halfspace_vars(lambda,shear_modulus)
  call test_value(a,0.66666666666666663d0,'halfspace_vars(): a')
  call test_value(CA1,0.16666666666666669d0,'halfspace_vars(): CA1')
  call test_value(CA2,0.33333333333333331d0,'halfspace_vars(): CA2')
  call test_value(CB,0.50000000000000011d0,'halfspace_vars(): CB')
  call test_value(CC,0.33333333333333337d0,'halfspace_vars(): CC')
  write(stdout,*) 'subroutine halfspace_vars() passed unit test'
  write(stdout,*)


  dip = 70.0d0
  call dip_vars(dip)
  call test_value(sd,0.93969262078590832d0,'dip_vars(): sd')
  call test_value(cd,0.34202014332566882d0,'dip_vars(): cd')
  call test_value(s2d,0.64278760968653947d0,'dip_vars(): s2d')
  call test_value(c2d,-0.76604444311897790d0,'dip_vars(): c2d')
  call test_value(sdsd,0.88302222155948884d0,'dip_vars(): sdsd')
  call test_value(cdcd,0.11697777844051105d0,'dip_vars(): cdcd')
  call test_value(cdsd,0.32139380484326974d0,'dip_vars(): cdsd')
  write(stdout,*) 'subroutine dip_vars() passed unit test'
  write(stdout,*)


  dip = 70.0d0
  call dip_vars(dip)
  sta_coord(1) = 2.0d0
  sta_coord(2) = 3.0d0
  sta_coord(3) = 0.0d0
  evdp = 4.0d0
  call pt_src_vars(sta_coord(1),sta_coord(2),sta_coord(3),evdp)
  call test_value(p,4.7848309131206399d0,'pt_src_vars(): p')
  call test_value(q,1.4509972890550495d0,'pt_src_vars(): q')
  call test_value(s,4.9925406015355298d0,'pt_src_vars(): s')
  call test_value(t,0.27301710938922441d0,'pt_src_vars(): t')
  call test_value(xx,4.0d0,'pt_src_vars(): xx')
  call test_value(xy,6.0d0,'pt_src_vars(): xy')
  call test_value(yy,9.0d0,'pt_src_vars(): yy')
  call test_value(dd,16.0d0,'pt_src_vars(): dd')
  call test_value(xd,8.0d0,'pt_src_vars(): xd')
  call test_value(xc,8.0d0,'pt_src_vars(): xc')
  call test_value(xq,2.9019945781100991d0,'pt_src_vars(): xq')
  call test_value(yq,4.3529918671651489d0,'pt_src_vars(): yq')
  call test_value(dq,5.8039891562201982d0,'pt_src_vars(): dq')
  call test_value(pq,6.9427766835248459d0,'pt_src_vars(): pq')
  call test_value(xs,9.9850812030710596d0,'pt_src_vars(): xs')
  call test_value(xt,0.54603421877844882d0,'pt_src_vars(): xt')
  call test_value(ys,14.977621804606589d0,'pt_src_vars(): ys')
  call test_value(yt,0.81905132816767323d0,'pt_src_vars(): yt')
  call test_value(qq,2.1053931328451032d0,'pt_src_vars(): qq')
  call test_value(R,5.3851648071345037d0,'pt_src_vars(): R')
  call test_value(R2,28.999999999999996d0,'pt_src_vars(): R2')
  call test_value(R3,156.16977940690060d0,'pt_src_vars(): R3')
  call test_value(R4,840.99999999999989d0,'pt_src_vars(): R4')
  call test_value(R5,4528.9236028001169d0,'pt_src_vars(): R5')
  call test_value(R7,131338.78448120339d0,'pt_src_vars(): R7')
  call test_value(Rd,9.3851648071345046d0,'pt_src_vars(): Rd')
  call test_value(A3,0.58620689655172409d0,'pt_src_vars(): A3')
  call test_value(A5,0.31034482758620685d0,'pt_src_vars(): A5')
  call test_value(A7,3.4482758620689502d-2,'pt_src_vars(): A7')
  call test_value(B3,6.8965517241379226d-2,'pt_src_vars(): B3')
  call test_value(B5,-0.55172413793103470d0,'pt_src_vars(): B5')
  call test_value(B7,-1.1724137931034484d0,'pt_src_vars(): B7')
  call test_value(C3,-0.65517241379310365d0,'pt_src_vars(): C3')
  call test_value(C5,-1.7586206896551726d0,'pt_src_vars(): C5')
  call test_value(C7,-2.8620689655172420d0,'pt_src_vars(): C7')
  call test_value(I1,4.4511857695225003d-003,'pt_src_vars(): I1')
  call test_value(I2,1.4062132997373915d-003,'pt_src_vars(): I2')
  call test_value(I3,1.1400361746955842d-002,'pt_src_vars(): I3')
  call test_value(I4,-6.4425920722921085d-003,'pt_src_vars(): I4')
  call test_value(I5,1.5490988410148277d-002,'pt_src_vars(): I5')
  call test_value(J1,-2.3037701706801377d-003,'pt_src_vars(): J1')
  call test_value(J2,8.3831231772032186d-004,'pt_src_vars(): J2')
  call test_value(J3,2.9153389890690740d-003,'pt_src_vars(): J3')
  call test_value(J4,-1.6706841541556938d-003,'pt_src_vars(): J4')
  call test_value(K1,-1.7024420377114214d-003,'pt_src_vars(): K1')
  call test_value(K2,1.3075030688791331d-004,'pt_src_vars(): K2')
  call test_value(K3,-5.4300227400023550d-003,'pt_src_vars(): K3')
  call test_value(U2,0.18917678161950324d0,'pt_src_vars(): U2')
  call test_value(V2,1.4014492135054364d0,'pt_src_vars(): V2')
  call test_value(W2,1.1288694024054116d0,'pt_src_vars(): W2')
  call test_value(U3,1.3427079288808756d0,'pt_src_vars(): U3')
  call test_value(V3,5.0611389600960157d0,'pt_src_vars(): V3')
  call test_value(W3,1.6847280722065445d0,'pt_src_vars(): W3')
  write(stdout,*) 'subroutine pt_src_vars() passed unit test'
  write(stdout,*)


  dip = 70.0d0
  call dip_vars(dip)
  sta_coord(1) = 2.0d0
  sta_coord(2) = 3.0d0
  sta_coord(3) = 0.0d0
  evdp = 4.0d0
  len = 3.0d0
  wid = 2.0d0
  call rect_src_coords(sta_coord(1),sta_coord(2),sta_coord(3),evdp,wid,len)
  call test_value(p,4.7848309131206399d0,'rect_src_coords(): p')
  call test_value(q,1.4509972890550495d0,'rect_src_coords(): q')
  call test_value(eta_vec(1),3.7848309131206399d0,'rect_src_coords(): eta_vec(1)')
  call test_value(eta_vec(2),5.7848309131206399d0,'rect_src_coords(): eta_vec(2)')
  call test_value(ksi_vec(1),0.5000000000000000d0,'rect_src_coords(): ksi_vec(1)')
  call test_value(ksi_vec(2),3.5000000000000000d0,'rect_src_coords(): ksi_vec(2)')
  write(stdout,*) 'subroutine rect_src_coords() passed unit test'
  call rect_src_vars(ksi_vec(1),eta_vec(1))
  call test_value(ybar,2.6579798566743307d0,'rect_src_vars(): ybar')
  call test_value(dbar,3.0603073792140916d0,'rect_src_vars(): cbar')
  call test_value(cbar,3.0603073792140916d0,'rect_src_vars(): dbar')
  call test_value(R,4.0841569722231199d0,'rect_src_vars(): R')
  call test_value(R2,16.680338173758720d0,'rect_src_vars(): R2')
  call test_value(R3,68.125119451396145d0,'rect_src_vars(): R3')
  ! call test_value(R4,840.99999999999989d0,'rect_src_vars(): R4')
  call test_value(R5,1136.3500305769958d0,'rect_src_vars(): R5')
  ! call test_value(R7,131338.78448120339d0,'rect_src_vars(): R7')
  call test_value(Rd,7.1444643514372110d0,'rect_src_vars(): Rd')
  call test_value(qq,2.1053931328451032d0,'rect_src_vars(): qq')
  call test_value(TH,0.30910022531517090d0,'rect_src_vars(): TH')
  call test_value(Rk,4.5841569722231199d0,'rect_src_vars(): Rk')
  call test_value(X11,5.3411908204318771d-002,'rect_src_vars(): X11')
  call test_value(X32,6.0549197882804537d-003,'rect_src_vars(): X32')
  call test_value(X53,1.3937426671487758d-003,'rect_src_vars(): X53')
  call test_value(logRk,1.5226062223303141d0,'rect_src_vars(): logRk')
  call test_value(Re,7.8689878853437598d0,'rect_src_vars(): Re')
  call test_value(Y11,3.1115637101260175d-002,'rect_src_vars(): Y11')
  call test_value(Y32,2.8335909219300561d-003,'rect_src_vars(): Y32')
  call test_value(Y53,5.6987957651285224d-004,'rect_src_vars(): Y53')
  call test_value(logRe,2.0629294500096305d0,'rect_src_vars(): logRe')
  call test_value(Y0,3.0407239370777661d-002,'rect_src_vars(): Y0')
  call test_value(I1,-24.464842377362366d0,'rect_src_vars(): I1')
  call test_value(I2,1.2599616181546565d0,'rect_src_vars(): I2')
  call test_value(I3,-0.75170987457349314d0,'rect_src_vars(): I3')
  call test_value(I4,26.009469280193130d0,'rect_src_vars(): I4')
  call test_value(Z32,1.2387401944116128d-002,'rect_src_vars(): Z53')
  call test_value(Z53,2.1980039131027692d-003,'rect_src_vars(): Z32')
  call test_value(Z0,1.1837900965840436d-002,'rect_src_vars(): Z0')
  call test_value(D11,3.4271088684950102d-002,'rect_src_vars(): D11')
  call test_value(K1,7.3562247838898460d-003,'rect_src_vars(): K1')
  call test_value(K2,0.11862078136010840d0,'rect_src_vars(): K2')
  call test_value(K3,-0.13432880842342962d0,'rect_src_vars(): K3')
  call test_value(K4,-1.5915028157424500d-003,'rect_src_vars(): K4')
  call test_value(J1,-3.6668885560592841d-002,'rect_src_vars(): J1')
  call test_value(J2,6.3749960046040766d-003,'rect_src_vars(): J2')
  call test_value(J3,3.9930048199627510d-003,'rect_src_vars(): J3')
  call test_value(J4,-1.3985998433743778d-002,'rect_src_vars(): J4')
  call test_value(J6,-1.1485677136015674d-002,'rect_src_vars(): J6')
  call test_value(E2,0.17347006835163784d0,'rect_src_vars(): E2')
  call test_value(F2,4.5587540534369766d-002,'rect_src_vars(): F2')
  call test_value(G2,7.7029412273117825d-002,'rect_src_vars(): G2')
  call test_value(H2,2.8818645950405460d-002,'rect_src_vars(): H2')
  call test_value(P2,8.8840473466618580d-003,'rect_src_vars(): P2')
  call test_value(Q2,1.9608281698289438d-003,'rect_src_vars(): Q2')
  call test_value(E3,0.14892464696612687d0,'rect_src_vars(): E3')
  call test_value(F3,3.9258435961432761d-002,'rect_src_vars(): F3')
  call test_value(G3,6.3422754458309291d-002,'rect_src_vars(): G3')
  call test_value(H3,2.4055253239807024d-002,'rect_src_vars(): H3')
  call test_value(P3,1.2387401944116128d-002,'rect_src_vars(): P3')
  call test_value(Q3,1.7300629665110766d-002,'rect_src_vars(): Q3')
  write(stdout,*) 'subroutine rect_src_vars() passed unit test'
  write(stdout,*)


  ! Half-space
  lambda = 40.0d0
  shear_modulus = 40.0d0
  call halfspace_vars(lambda,shear_modulus)
  ! Fault dip
  dip = 70.0d0
  call dip_vars(dip)
  ! Geometry
  sta_coord(1) = -12.0d0
  sta_coord(2) = 6.2d0
  sta_coord(3) = 13.0d0
  evdp = 7.0d0
  wid = 3.2d0
  len = 6.5d0
  call rect_src_coords(sta_coord(1),sta_coord(2),-sta_coord(3),evdp,wid,len)
  call rect_src_vars(ksi_vec(2),eta_vec(1))
  disp = uA_ss_rect()
  call test_value(disp(1), 0.72547327384003113d0,'uA_ss_rect()(1)')
  call test_value(disp(2),-1.5927062399294267d-2,'uA_ss_rect()(2)')
  call test_value(disp(3), 0.61666031833564772d0,'uA_ss_rect()(3)')
  disp = uB_ss_rect()
  call test_value(disp(1),-11.993002387482605d0,'uB_ss_rect()(1)')
  call test_value(disp(2),0.11463965563116107d0,'uB_ss_rect()(2)')
  call test_value(disp(3),-0.99745107154736778d0,'uB_ss_rect()(3)')
  disp = uC_ss_rect()
  call test_value(disp(1),-1.4461853668621687d-003,'uC_ss_rect()(1)')
  call test_value(disp(2),5.0207492407645738d-003,'uC_ss_rect()(2)')
  call test_value(disp(3),-2.0079399191054680d-002,'uC_ss_rect()(3)')
  call rect_src_coords(sta_coord(1),sta_coord(2),-sta_coord(3),evdp,wid,len)
  call rect_src_vars(ksi_vec(1),eta_vec(2))
  disp = uA_ds_rect()
  call test_value(disp(1),-1.2424808371844004d-2,'uA_ds_rect()(1)')
  call test_value(disp(2), 0.72190403451379292d0,'uA_ds_rect()(2)')
  call test_value(disp(3), 0.41256779800532334d0,'uA_ds_rect()(3)')
  disp = uB_ds_rect()
  call test_value(disp(1),-0.24780594882338475d0,'uB_ds_rect()(1)')
  call test_value(disp(2),-1.3701174902198523d0,'uB_ds_rect()(2)')
  call test_value(disp(3),-3.8188429813460623d0,'uB_ds_rect()(3)')
  disp = uC_ds_rect()
  call test_value(disp(1), 5.1793376586471750d-003,'uC_ds_rect()(1)')
  call test_value(disp(2), 8.6684162903643731d-003,'uC_ds_rect()(2)')
  call test_value(disp(3),-7.2807977939472246d-002,'uC_ds_rect()(3)')
  write(stdout,*) 'displacement functions for rectangular source passed unit tests'
  write(stdout,*)


  ! Half-space
  lambda = 40.0d0
  shear_modulus = 40.0d0
  call halfspace_vars(lambda,shear_modulus)
  ! Fault dip
  dip = 70.0d0
  call dip_vars(dip)
  ! Geometry
  sta_coord(1) = -12.0d0
  sta_coord(2) = 6.2d0
  sta_coord(3) = 13.0d0
  evdp = 7.0d0
  wid = 3.2d0
  len = 6.5d0
  call rect_src_coords(sta_coord(1),sta_coord(2),-sta_coord(3),evdp,wid,len)
  call rect_src_vars(ksi_vec(1),eta_vec(1))
  ! X-derivatives of components of displacement
  disp = duAdx_ss_rect()
  call test_value(disp(1), 3.4306426504990750d-004,'duAdx_ss_rect(): (1)')
  call test_value(disp(2),-3.4508486032947076d-004,'duAdx_ss_rect(): (2)')
  call test_value(disp(3),-2.3607160020503188d-003,'duAdx_ss_rect(): (3)')
  disp = duAdx_ds_rect()
  call test_value(disp(1),-3.4508486032947076d-004,'duAdx_ds_rect(): (1)')
  call test_value(disp(2), 9.0562367901465455d-004,'duAdx_ds_rect(): (2)')
  call test_value(disp(3), 6.7897701068703959d-003,'duAdx_ds_rect(): (3)')
  ! disp = duAdx_ts_rect()
  ! write(0,*) 'duAdx_ts_rect()',disp
  disp = duBdx_ss_rect()
  call test_value(disp(1), 6.9599894428799694d-004,'duBdx_ss_rect(): (1)')
  call test_value(disp(2), 1.9194218761293762d-003,'duBdx_ss_rect(): (2)')
  call test_value(disp(3), 2.0942277171944572d-003,'duBdx_ss_rect(): (3)')
  disp = duBdx_ds_rect()
  call test_value(disp(1), 2.7417691868664293d-003,'duBdx_ds_rect(): (1)')
  call test_value(disp(2),-5.1586589356422364d-003,'duBdx_ds_rect(): (2)')
  call test_value(disp(3),-6.7076617173327598d-004,'duBdx_ds_rect(): (3)')
  ! disp = duBdx_ts_rect()
  ! write(0,*) 'duBdx_ts_rect()',disp
  disp = duCdx_ss_rect()
  call test_value(disp(1), 4.2126939364670769d-005,'duCdx_ss_rect(): (1)')
  call test_value(disp(2), 1.1209090419535319d-004,'duCdx_ss_rect(): (2)')
  call test_value(disp(3),-3.4475338947556833d-004,'duCdx_ss_rect(): (3)')
  disp = duCdx_ds_rect()
  call test_value(disp(1), 1.6966631277932049d-004,'duCdx_ds_rect(): (1)')
  call test_value(disp(2),-1.4986921869860253d-004,'duCdx_ds_rect(): (2)')
  call test_value(disp(3), 1.1331156239589844d-003,'duCdx_ds_rect(): (3)')
  ! disp = duCdx_ts_rect()
  ! write(0,*) 'duCdx_ts_rect()',disp
  ! Y-derivatives of components of displacement
  disp = duAdy_ss_rect()
  call test_value(disp(1),2.8890966553217295d-002,'duAdy_ss_rect(): (1)')
  call test_value(disp(2),1.2845371264737372d-002,'duAdy_ss_rect(): (2)')
  call test_value(disp(3),2.7617711810937568d-003,'duAdy_ss_rect(): (3)')
  disp = duAdy_ds_rect()
  call test_value(disp(1), 1.2845371264737372d-002,'duAdy_ds_rect(): (1)')
  call test_value(disp(2), 6.0051670894042106d-002,'duAdy_ds_rect(): (2)')
  call test_value(disp(3), 6.8785272909704029d-003,'duAdy_ds_rect(): (3)')
  ! disp = duAdy_ts_rect()
  disp = duBdy_ss_rect()
  call test_value(disp(1),-5.4892739626334362d-002,'duBdy_ss_rect(): (1)')
  call test_value(disp(2),-2.7969238876413290d-002,'duBdy_ss_rect(): (2)')
  call test_value(disp(3),-3.1020145458715599d-003,'duBdy_ss_rect(): (3)')
  disp = duBdy_ds_rect()
  call test_value(disp(1),-3.8965904427687062d-002,'duBdy_ds_rect(): (1)')
  call test_value(disp(2),-0.14704676790658636d0,'duBdy_ds_rect(): (2)')
  call test_value(disp(3),-9.1052106567444321d-003,'duBdy_ds_rect(): (3)')
  ! disp = duBdy_ts_rect()
  disp = duCdy_ss_rect()
  call test_value(disp(1),3.6229956319952023d-004,'duCdy_ss_rect(): (1)')
  call test_value(disp(2),2.7645676601602590d-004,'duCdy_ss_rect(): (2)')
  call test_value(disp(3),5.2251433363535108d-004,'duCdy_ss_rect(): (3)')
  disp = duCdy_ds_rect()
  call test_value(disp(1),-1.1160386819908798d-003,'duCdy_ds_rect(): (1)')
  call test_value(disp(2),-6.3814177065406037d-004,'duCdy_ds_rect(): (2)')
  call test_value(disp(3), 2.7583611996049840d-003,'duCdy_ds_rect(): (3)')
  ! disp = duCdy_ts_rect()
  ! Z-derivatives of components of displacement
  disp = duAdz_ss_rect()
  call test_value(disp(1), 8.5466041996219398d-003,'duAdz_ss_rect(): (1)')
  call test_value(disp(2), 4.2102278432794181d-003,'duAdz_ss_rect(): (2)')
  call test_value(disp(3),-6.2203199241805318d-003,'duAdz_ss_rect(): (3)')
  disp = duAdz_ds_rect()
  call test_value(disp(1), 4.2102278432794181d-003,'duAdz_ds_rect(): (1)')
  call test_value(disp(2), 1.7605900194629188d-002,'duAdz_ds_rect(): (2)')
  call test_value(disp(3),-1.2506584807178457d-002,'duAdz_ds_rect(): (3)')
  ! disp = duAdy_ts_rect()
  disp = duBdz_ss_rect()
  call test_value(disp(1),-1.7349770792071152d-002,'duBdz_ss_rect(): (1)')
  call test_value(disp(2),-1.0130288362773653d-002,'duBdz_ss_rect(): (2)')
  call test_value(disp(3), 1.0421179348291247d-002,'duBdz_ss_rect(): (3)')
  disp = duBdz_ds_rect()
  call test_value(disp(1),-9.6899784763487366d-003,'duBdz_ds_rect(): (1)')
  call test_value(disp(2),-4.0686184682977758d-002,'duBdz_ds_rect(): (2)')
  call test_value(disp(3),-2.2296367147699551d-003,'duBdz_ds_rect(): (3)')
  ! disp = duBdy_ts_rect()
  disp = duCdz_ss_rect()
  call test_value(disp(1),-3.6434767894171044d-005,'duCdz_ss_rect(): (1)')
  call test_value(disp(2), 2.3748100977727307d-004,'duCdz_ss_rect(): (2)')
  call test_value(disp(3),-7.0849983200360817d-004,'duCdz_ss_rect(): (3)')
  disp = duCdz_ds_rect()
  call test_value(disp(1),-1.5620330749121738d-004,'duCdz_ds_rect(): (1)')
  call test_value(disp(2), 5.7029847969910998d-004,'duCdz_ds_rect(): (2)')
  call test_value(disp(3),-5.4083583036719013d-003,'duCdz_ds_rect(): (3)')
  ! disp = duCdz_ts_rect()
  write(stdout,*) 'displacement derivative functions for rectangular source passed unit tests'
  write(stdout,*)

#endif

!----
! Point source tests
!----
! Point source example from Okada (1985)
lambda = 40.0d0            ! lambda = shear_modulus implies
shear_modulus = 40.0d0     !     a Poisson's ratio of 0.25
dip = 70.0d0
sta_coord(1) = 2.0d0
sta_coord(2) = 3.0d0
sta_coord(3) = 0.0d0
evdp = 4.0d0

! Strike-slip
moment(1) = 40.0d0
moment(2) =  0.0d0
moment(3) =  0.0d0
call o92_pt_disp(disp,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(disp(1),-9.4474186064630631d-4,'o92_pt_disp(): disp(1) (strike-slip)')
call test_value(disp(2),-1.0229885427654646d-3,'o92_pt_disp(): disp(2) (strike-slip)')
call test_value(disp(3),-7.4201246532591128d-4,'o92_pt_disp(): disp(3) (strike-slip)')
write(stdout,*) 'subroutine o92_pt_disp() passed strike-slip unit test'
call o92_pt_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz,&
                     sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(uxx,-2.2862094355883986d-004,'o92_pt_partials(): uxx (strike-slip)')
call test_value(uxy,-1.4246398245262621d-004,'o92_pt_partials(): uxy (strike-slip)')
call test_value(uxz, 6.2590905699873601d-005,'o92_pt_partials(): uxz (strike-slip)')
call test_value(uyx,-2.0511017314976633d-004,'o92_pt_partials(): uyx (strike-slip)')
call test_value(uyy,-3.0067825210386353d-004,'o92_pt_partials(): uyy (strike-slip)')
call test_value(uyz, 1.6933001353308135d-004,'o92_pt_partials(): uyz (strike-slip)')
call test_value(uzx,-6.2590905699873859d-005,'o92_pt_partials(): uzx (strike-slip)')
call test_value(uzy,-1.6933001353308129d-004,'o92_pt_partials(): uzy (strike-slip)')
call test_value(uzz, 1.7643306522090115d-004,'o92_pt_partials(): uzz (strike-slip)')
write(stdout,*) 'subroutine o92_pt_partials() passed strike-slip unit test'
call o92_pt_strain(strain,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(strain(1,1),-2.2862094355883986d-004,'o92_pt_strain(): strain(1,1) (strike-slip)')
call test_value(strain(1,2),-1.7378707780119627d-004,'o92_pt_strain(): strain(1,2) (strike-slip)')
call test_value(strain(1,3),-1.2874900798265365d-019,'o92_pt_strain(): strain(1,3) (strike-slip)')
call test_value(strain(2,1),-1.7378707780119627d-004,'o92_pt_strain(): strain(2,1) (strike-slip)')
call test_value(strain(2,2),-3.0067825210386353d-004,'o92_pt_strain(): strain(2,2) (strike-slip)')
call test_value(strain(2,3), 2.7105054312137611d-020,'o92_pt_strain(): strain(2,3) (strike-slip)')
call test_value(strain(3,1),-1.2874900798265365d-019,'o92_pt_strain(): strain(3,1) (strike-slip)')
call test_value(strain(3,2), 2.7105054312137611d-020,'o92_pt_strain(): strain(3,2) (strike-slip)')
call test_value(strain(3,3), 1.7643306522090115d-004,'o92_pt_strain(): strain(3,3) (strike-slip)')
write(stdout,*) 'subroutine o92_pt_strain() passed strike-slip unit test'

! Dip-slip
moment(1) =  0.0d0
moment(2) = 40.0d0
moment(3) =  0.0d0
call o92_pt_disp(disp,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(disp(1),-1.1723215364075840d-3,'o92_pt_disp(): disp(1) (dip-slip)')
call test_value(disp(2),-2.0819985410726201d-3,'o92_pt_disp(): disp(2) (dip-slip)')
call test_value(disp(3),-2.5315947108865925d-3,'o92_pt_disp(): disp(3) (dip-slip)')
write(stdout,*) 'subroutine o92_pt_disp() passed dip-slip unit test'
call o92_pt_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz,&
                     sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(uxx,-1.5259408689661085d-004,'o92_pt_partials(): uxx (dip-slip)')
call test_value(uxy,-3.5441800891339358d-004,'o92_pt_partials(): uxy (dip-slip)')
call test_value(uxz,-8.7070529207727333d-004,'o92_pt_partials(): uxz (dip-slip)')
call test_value(uyx, 6.9826591999645748d-004,'o92_pt_partials(): uyx (dip-slip)')
call test_value(uyy,-1.1537526800074350d-003,'o92_pt_partials(): uyy (dip-slip)')
call test_value(uyz, 6.3453607666584869d-004,'o92_pt_partials(): uyz (dip-slip)')
call test_value(uzx, 8.7070529207727322d-004,'o92_pt_partials(): uzx (dip-slip)')
call test_value(uzy,-6.3453607666584847d-004,'o92_pt_partials(): uzy (dip-slip)')
call test_value(uzz, 4.3544892230134886d-004,'o92_pt_partials(): uzz (dip-slip)')
write(stdout,*) 'subroutine o92_pt_partials() passed dip-slip unit test'
call o92_pt_strain(strain,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(strain(1,1),-1.5259408689661085d-004,'o92_pt_strain(): strain(1,1) (dip-slip)')
call test_value(strain(1,2), 1.7192395554153195d-004,'o92_pt_strain(): strain(1,2) (dip-slip)')
call test_value(strain(1,3),-5.4210108624275222d-020,'o92_pt_strain(): strain(1,3) (dip-slip)')
call test_value(strain(2,1), 1.7192395554153195d-004,'o92_pt_strain(): strain(2,1) (dip-slip)')
call test_value(strain(2,2),-1.1537526800074350d-003,'o92_pt_strain(): strain(2,2) (dip-slip)')
call test_value(strain(2,3), 1.0842021724855044d-019,'o92_pt_strain(): strain(2,3) (dip-slip)')
call test_value(strain(3,1),-5.4210108624275222d-020,'o92_pt_strain(): strain(3,1) (dip-slip)')
call test_value(strain(3,2), 1.0842021724855044d-019,'o92_pt_strain(): strain(3,2) (dip-slip)')
call test_value(strain(3,3), 4.3544892230134886d-004,'o92_pt_strain(): strain(3,3) (dip-slip)')
write(stdout,*) 'subroutine o92_pt_strain() passed dip-slip unit test'

! Tensile-slip
moment(1) =  0.0d0
moment(2) =  0.0d0
moment(3) = 40.0d0
call o92_pt_disp(disp,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(disp(1),-3.5716333691232548d-4,'o92_pt_disp(): disp(1) (tensile-slip)')
call test_value(disp(2), 3.5310854913572190d-4,'o92_pt_disp(): disp(2) (tensile-slip)')
call test_value(disp(3),-2.0068126969335952d-4,'o92_pt_disp(): disp(3) (tensile-slip)')
write(stdout,*) 'subroutine o92_pt_disp() passed tensile-slip unit test'
call o92_pt_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz,&
                     sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(uxx,-1.3597209844017221d-004,'o92_pt_partials(): uxx (tensile-slip)')
call test_value(uxy, 5.0725453276876352d-004,'o92_pt_partials(): uxy (tensile-slip)')
call test_value(uxz,-7.5405344488748576d-005,'o92_pt_partials(): uxz (tensile-slip)')
call test_value(uyx,-6.7733492861951009d-005,'o92_pt_partials(): uyx (tensile-slip)')
call test_value(uyy, 6.8111287703614769d-004,'o92_pt_partials(): uyy (tensile-slip)')
call test_value(uyz,-8.1037165493104131d-004,'o92_pt_partials(): uyz (tensile-slip)')
call test_value(uzx, 7.5405344488748590d-005,'o92_pt_partials(): uzx (tensile-slip)')
call test_value(uzy, 8.1037165493104120d-004,'o92_pt_partials(): uzy (tensile-slip)')
call test_value(uzz,-1.8171359286532499d-004,'o92_pt_partials(): uzz (tensile-slip)')
write(stdout,*) 'subroutine o92_pt_partials() passed tensile-slip unit test'
call o92_pt_strain(strain,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(strain(1,1),-1.3597209844017221d-004,'o92_pt_strain(): strain(1,1) (tensile-slip)')
call test_value(strain(1,2), 2.1976051995340625d-004,'o92_pt_strain(): strain(1,2) (tensile-slip)')
call test_value(strain(1,3), 6.7762635780344027d-021,'o92_pt_strain(): strain(1,3) (tensile-slip)')
call test_value(strain(2,1), 2.1976051995340625d-004,'o92_pt_strain(): strain(2,1) (tensile-slip)')
call test_value(strain(2,2), 6.8111287703614769d-004,'o92_pt_strain(): strain(2,2) (tensile-slip)')
call test_value(strain(2,3),-5.4210108624275222d-020,'o92_pt_strain(): strain(2,3) (tensile-slip)')
call test_value(strain(3,1), 6.7762635780344027d-021,'o92_pt_strain(): strain(3,1) (tensile-slip)')
call test_value(strain(3,2),-5.4210108624275222d-020,'o92_pt_strain(): strain(3,2) (tensile-slip)')
call test_value(strain(3,3),-1.8171359286532499d-004,'o92_pt_strain(): strain(3,3) (tensile-slip)')
write(stdout,*) 'subroutine o92_pt_strain() passed tensile-slip unit test'
write(stdout,*) 'point source example from Okada (1985) passed unit test'
write(stdout,*)

! Point source, station depth greater than zero, oblique slip
lambda = 40.0d0
shear_modulus = 40.0d0
sta_coord(1) = -12.0d0
sta_coord(2) = 6.2d0
sta_coord(3) = 13.0d0
evdp = 7.0d0
moment(1) = -27.0d0
moment(2) = 13.0d0
moment(3) = 0.0d0

! Displacement
call o92_pt_disp(disp,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(disp(1), 1.9245175992850267d-004,'o92_pt_disp(): disp(1) (z>0; oblique-slip)')
call test_value(disp(2),-1.9098044959056914d-004,'o92_pt_disp(): disp(2) (z>0; oblique-slip)')
call test_value(disp(3), 8.1576627770328566d-005,'o92_pt_disp(): disp(3) (z>0; oblique-slip)')

! Partial derivatives of displacement
call o92_pt_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz,&
                     sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(uxx, 2.0874574675642794d-005,'o92_pt_partials(): uxx (z>0; oblique-slip)')
call test_value(uxy, 2.3620493327193696d-006,'o92_pt_partials(): uxy (z>0; oblique-slip)')
call test_value(uxz, 2.2899507819320573d-005,'o92_pt_partials(): uxz (z>0; oblique-slip)')
call test_value(uyx,-2.1572559953257881d-005,'o92_pt_partials(): uyx (z>0; oblique-slip)')
call test_value(uyy,-6.8032602399881741d-006,'o92_pt_partials(): uyy (z>0; oblique-slip)')
call test_value(uyz,-2.1249535452408235d-005,'o92_pt_partials(): uyz (z>0; oblique-slip)')
call test_value(uzx, 1.2032014067111290d-005,'o92_pt_partials(): uzx (z>0; oblique-slip)')
call test_value(uzy,-6.7523163638339241d-006,'o92_pt_partials(): uzy (z>0; oblique-slip)')
call test_value(uzz,-1.7420298917986620d-006,'o92_pt_partials(): uzz (z>0; oblique-slip)')

! Strain tensor
call o92_pt_strain(strain,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(strain(1,1), 2.0874574675642794d-5,'o92_pt_strain(): e11 (z>0; oblique-slip)')
call test_value(strain(1,2),-9.6052553102692557d-6,'o92_pt_strain(): e12 (z>0; oblique-slip)')
call test_value(strain(1,3), 1.7465760943215932d-5,'o92_pt_strain(): e13 (z>0; oblique-slip)')
call test_value(strain(2,1),-9.6052553102692557d-6,'o92_pt_strain(): e21 (z>0; oblique-slip)')
call test_value(strain(2,2),-6.8032602399881741d-6,'o92_pt_strain(): e22 (z>0; oblique-slip)')
call test_value(strain(2,3),-1.4000925908121080d-5,'o92_pt_strain(): e23 (z>0; oblique-slip)')
call test_value(strain(3,1), 1.7465760943215932d-5,'o92_pt_strain(): e31 (z>0; oblique-slip)')
call test_value(strain(3,2),-1.4000925908121080d-5,'o92_pt_strain(): e32 (z>0; oblique-slip)')
call test_value(strain(3,3),-1.7420298917986620d-6,'o92_pt_strain(): e33 (z>0; oblique-slip)')
write(stdout,*) 'point source example with depth>0, oblique slip passed unit test'
write(stdout,*)

!----
! Rectangular source
!----
! Example from Okada (1985)
lambda = 40.0d0            ! lambda = shear_modulus implies
shear_modulus = 40.0d0     !     a Poisson's ratio of 0.25
dip = 70.0d0
sta_coord(1) = 2.0d0
sta_coord(2) = 3.0d0
sta_coord(3) = 0.0d0
evdp = 4.0d0
len = 3.0d0
wid = 2.0d0

! Original Okada derivation has fault origin at the bottom left corner, but I put the fault
! center at the coordinates. To compare with Okada (1985), shift fault.
sta_coord(1) = sta_coord(1)-len/2.0d0
sta_coord(2) = sta_coord(2)-wid/2.0d0*dcos(dip*d2r)
evdp = evdp-wid/2.0d0*dsin(dip*d2r)

! Strike-slip
slip_vec(1) = 1.0d0
slip_vec(2) = 0.0d0
slip_vec(3) = 0.0d0
call o92_rect_disp(disp,sta_coord,evdp,dip,slip_vec,wid,len,lambda,shear_modulus)
call test_value(disp(1),-8.6891650042561549d-003,'o92_rect_disp(): disp(1) (strike-slip)')
call test_value(disp(2),-4.2975821897415613d-003,'o92_rect_disp(): disp(2) (strike-slip)')
call test_value(disp(3),-2.7474058276388967d-003,'o92_rect_disp(): disp(3) (strike-slip)')
write(stdout,*) 'subroutine o92_rect_disp() passed strike-slip unit test'

! Dip-slip
slip_vec(1) = 0.0d0
slip_vec(2) = 1.0d0
slip_vec(3) = 0.0d0
call o92_rect_disp(disp,sta_coord,evdp,dip,slip_vec,wid,len,lambda,shear_modulus)
call test_value(disp(1),-4.6823487628354045d-003,'o92_rect_disp(): disp(1) (dip-slip)')
call test_value(disp(2),-3.5267267968718020d-002,'o92_rect_disp(): disp(2) (dip-slip)')
call test_value(disp(3),-3.5638557673268470d-002,'o92_rect_disp(): disp(3) (dip-slip)')
write(stdout,*) 'subroutine o92_rect_disp() passed dip-slip unit test'

! Tensile-slip
slip_vec(1) = 0.0d0
slip_vec(2) = 0.0d0
slip_vec(3) = 1.0d0
call o92_rect_disp(disp,sta_coord,evdp,dip,slip_vec,wid,len,lambda,shear_modulus)
call test_value(disp(1),-2.6599600964428191d-004,'o92_rect_disp(): disp(1) (tensile-slip)')
call test_value(disp(2), 1.0564074876983659d-002,'o92_rect_disp(): disp(2) (tensile-slip)')
call test_value(disp(3), 3.2141931142207676d-003,'o92_rect_disp(): disp(3) (tensile-slip)')
write(stdout,*) 'subroutine o92_rect_disp() passed tensile-slip unit test'
write(stdout,*) 'rectangular source example from Okada (1985) passed unit test'
write(stdout,*)


! Example with rectangular source, station depth greater than zero
lambda = 40.0d0
shear_modulus = 40.0d0
dip = 70.0d0
sta_coord(1) = -12.0d0
sta_coord(2) = 6.2d0
sta_coord(3) = 13.0d0
evdp = 7.0d0
wid = 3.2d0
len = 6.5d0

! Strike-slip displacement
slip_vec(1) = 1.0d0
slip_vec(2) = 0.0d0
slip_vec(3) = 0.0d0
call o92_rect_disp(disp,sta_coord,evdp,dip,slip_vec,wid,len,lambda,shear_modulus)
call test_value(disp(1),-6.6858712769214453d-003,'o92_rect_disp(): disp(1) (z>0; strike-slip)')
call test_value(disp(2), 5.6833739558216096d-003,'o92_rect_disp(): disp(2) (z>0; strike-slip)')
call test_value(disp(3),-3.8337220569204383d-003,'o92_rect_disp(): disp(3) (z>0; strike-slip)')
write(stdout,*)

! Dip-slip displacement
slip_vec(1) = 0.0d0
slip_vec(2) = 1.0d0
slip_vec(3) = 0.0d0
call o92_rect_disp(disp,sta_coord,evdp,dip,slip_vec,wid,len,lambda,shear_modulus)
call test_value(disp(1),-1.4891278409154987d-003,'o92_rect_disp(): disp(1) (z>0; dip-slip)')
call test_value(disp(2),-6.6628749834942423d-004,'o92_rect_disp(): disp(2) (z>0; dip-slip)')
call test_value(disp(3),-2.6808905339339140d-003,'o92_rect_disp(): disp(3) (z>0; dip-slip)')
write(stdout,*)

! Strike-slip partial derivatives of displacement
slip_vec(1) = 1.0d0
slip_vec(2) = 0.0d0
slip_vec(3) = 0.0d0
call o92_rect_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz, &
                       sta_coord,evdp,dip,slip_vec,wid,len,lambda,shear_modulus)
call test_value(uxx,-7.4224535411782154d-004,'o92_rect_partials(): uxx (z>0; strike-slip)')
call test_value(uyx, 7.0705631705157462d-004,'o92_rect_partials(): uyx (z>0; strike-slip)')
call test_value(uzx,-5.7940045023416157d-004,'o92_rect_partials(): uzx (z>0; strike-slip)')
call test_value(uxy, 6.4347351147793649d-005,'o92_rect_partials(): uxy (z>0; strike-slip)')
call test_value(uyy, 1.5409553602864544d-004,'o92_rect_partials(): uyy (z>0; strike-slip)')
call test_value(uzy, 2.2074703722434528d-004,'o92_rect_partials(): uzy (z>0; strike-slip)')
call test_value(uxz,-6.3176025002264938d-004,'o92_rect_partials(): uxz (z>0; strike-slip)')
call test_value(uyz, 5.5337876385493990d-004,'o92_rect_partials(): uyz (z>0; strike-slip)')
call test_value(uzz, 1.0685389278657297d-004,'o92_rect_partials(): uzz (z>0; strike-slip)')

slip_vec(1) = 0.0d0
slip_vec(2) = 1.0d0
slip_vec(3) = 0.0d0
call o92_rect_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz, &
                       sta_coord,evdp,dip,slip_vec,wid,len,lambda,shear_modulus)
call test_value(uxx,-2.7933167653439758d-004,'o92_rect_partials(): uxx (z>0; dip-slip)')
call test_value(uyx, 9.6965877919922243d-005,'o92_rect_partials(): uyx (z>0; dip-slip)')
call test_value(uzx,-4.7400064111388521d-004,'o92_rect_partials(): uzx (z>0; dip-slip)')
call test_value(uxy, 1.8756328884862554d-004,'o92_rect_partials(): uxy (z>0; dip-slip)')
call test_value(uyy,-6.2549348110294185d-005,'o92_rect_partials(): uyy (z>0; dip-slip)')
call test_value(uzy,-2.0932522925848662d-005,'o92_rect_partials(): uzy (z>0; dip-slip)')
call test_value(uxz, 2.4695902664444695d-004,'o92_rect_partials(): uxz (z>0; dip-slip)')
call test_value(uyz,-3.3607226779516899d-004,'o92_rect_partials(): uyz (z>0; dip-slip)')
call test_value(uzz, 1.6913995963324557d-004,'o92_rect_partials(): uzz (z>0; dip-slip)')

! slip_vec(1) = 0.0d0
! slip_vec(2) = 0.0d0
! slip_vec(3) = 1.0d0
! call o92_rect_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz, &
!                        sta_coord,evdp,dip,slip_vec,wid,len,lambda,shear_modulus)
! call test_value(uxx,-7.4224535411782154d-004,'o92_rect_partials(): uxx (z>0; tensile-slip)')
! call test_value(uyx, 7.0705631705157462d-004,'o92_rect_partials(): uyx (z>0; tensile-slip)')
! call test_value(uzx,-5.7940045023416157d-004,'o92_rect_partials(): uzx (z>0; tensile-slip)')
! call test_value(uxy, 6.4347351147793649d-005,'o92_rect_partials(): uxy (z>0; tensile-slip)')
! call test_value(uyy, 1.5409553602864544d-004,'o92_rect_partials(): uyy (z>0; tensile-slip)')
! call test_value(uzy, 2.2074703722434528d-004,'o92_rect_partials(): uzy (z>0; tensile-slip)')
! call test_value(uxz,-6.3176025002264938d-004,'o92_rect_partials(): uxz (z>0; tensile-slip)')
! call test_value(uyz, 5.5337876385493990d-004,'o92_rect_partials(): uyz (z>0; tensile-slip)')
! call test_value(uzz, 1.0685389278657297d-004,'o92_rect_partials(): uzz (z>0; tensile-slip)')

slip_vec(1) = 3.0d0
slip_vec(2) = -8.0d0
slip_vec(3) = 0.0d0
call o92_rect_strain(strain,sta_coord,evdp,dip,slip_vec,wid,len,lambda,shear_modulus)
call test_value(strain(1,1), 7.9173499217153817d-006,'o92_rect_strain(): e11 (z>0; oblique-slip)')
call test_value(strain(1,2), 1.8988835224861660d-005,'o92_rect_strain(): e12 (z>0; oblique-slip)')
call test_value(strain(1,3),-9.0857459250746377d-004,'o92_rect_strain(): e13 (z>0; oblique-slip)')
call test_value(strain(2,1), 1.8988835224861660d-005,'o92_rect_strain(): e21 (z>0; oblique-slip)')
call test_value(strain(2,2), 9.6268139296828981d-004,'o92_rect_strain(): e22 (z>0; oblique-slip)')
call test_value(strain(2,3), 2.5892078645029983d-003,'o92_rect_strain(): e23 (z>0; oblique-slip)')
call test_value(strain(3,1),-9.0857459250746377d-004,'o92_rect_strain(): e31 (z>0; oblique-slip)')
call test_value(strain(3,2), 2.5892078645029983d-003,'o92_rect_strain(): e32 (z>0; oblique-slip)')
call test_value(strain(3,3),-1.0325579987062454d-003,'o92_rect_strain(): e33 (z>0; oblique-slip)')

!----
! Compare with original code results
!----
! sta_coord(1) = 1111.949d0
! sta_coord(2) = -1111.949d0
! sta_coord(3) = 0.0d0
! evdp = 10000.0d0
! dip = 90.0d0
! wid = 4000.0d0
! len = 6000.0d0
! lambda = 4.0d10
! shear_modulus = 4.0d10
! slip_vec(1) = -1.0d0
! slip_vec(2) = 0.0d0
! slip_vec(3) = 0.0d0
! xin = sta_coord(1)
! yin = sta_coord(2)
! stdp = sta_coord(3)
! dipin = dip
! rakin = datan2(moment(2),moment(1))/d2r
! area = 1.0d0
! slip = dsqrt(moment(1)*moment(1)+moment(2)*moment(2))/shear_modulus/area
! dens = 1.0d0
! vp = dsqrt((lambda+2.0d0*shear_modulus)/dens)
! vs = dsqrt(shear_modulus/dens)
! call o92pt(ux,uy,uz,xin,yin,stdp,evdp,dipin,rakin,area,slip,vp,vs,dens)
! write(6,*) ux,uy,uz
! call o92ptstn(strain,xin,yin,stdp,evdp,dipin,rakin,area,slip,vp,vs,dens)
! write(6,*) strain(1,:)
! write(6,*) strain(2,:)
! write(6,*) strain(3,:)

! rakin = datan2(slip_vec(2),slip_vec(1))/d2r
! slip = dsqrt(slip_vec(1)*slip_vec(1)+slip_vec(2)*slip_vec(2))
! write(0,*) 'sta_coord:',sta_coord,'    x,y,z:',xin,yin,stdp
! write(0,*) 'evdp:',evdp
! write(0,*) 'dip:',dip,dipin
! write(0,*) 'slip_vec:',slip_vec,'    slip,rakin:',slip,rakin
! write(0,*) 'wid:',wid
! write(0,*) 'len:',len
! write(0,*) 'shear:',shear_modulus
! write(0,*) 'vp,vs,dens:',vp,vs,dens
! call o92_rect_disp(disp,sta_coord,evdp,dip,slip_vec,wid,len,lambda,shear_modulus)
! write(0,*) disp
! call o92rect(ux,uy,uz,xin,yin,stdp,evdp,dipin,rakin,wid,len,slip,vp,vs,dens)
! write(0,*) ux,uy,uz
! call o92rectstn(strain,xin,yin,stdp,evdp,dipin,rakin,wid,len,slip,vp,vs,dens)

write(stdout,*) 'okada92_module unit test passed'
end
