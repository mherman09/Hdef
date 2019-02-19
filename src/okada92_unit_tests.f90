program main
use okada92_module
implicit none
double precision :: disp(3), strain(3,3), sta_coord(3), evdp, dip, moment(4), lambda, shear_modulus
double precision :: slip_vec(3), wid, len
double precision :: ux,uy,uz,xin,yin,stdp,dipin,rakin,slip,vp,vs,dens
double precision :: uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz

! Example from Okada (1985) for surface deformation
lambda = 40.0d0
shear_modulus = 40.0d0
call halfspace_vars(lambda,shear_modulus)
call test_value(a,0.66666666666666663d0,'halfspace_vars(): a')
call test_value(CA1,0.16666666666666669d0,'halfspace_vars(): CA1')
call test_value(CA2,0.33333333333333331d0,'halfspace_vars(): CA2')
call test_value(CB,0.50000000000000011d0,'halfspace_vars(): CB')
call test_value(CC,0.33333333333333337d0,'halfspace_vars(): CC')
write(6,*) 'subroutine halfspace_vars() passed unit test'
write(6,*)

dip = 70.0d0
call dip_vars(dip)
call test_value(sd,0.93969262078590832d0,'dip_vars(): sd')
call test_value(cd,0.34202014332566882d0,'dip_vars(): cd')
call test_value(s2d,0.64278760968653947d0,'dip_vars(): s2d')
call test_value(c2d,-0.76604444311897790d0,'dip_vars(): c2d')
call test_value(sdsd,0.88302222155948884d0,'dip_vars(): sdsd')
call test_value(cdcd,0.11697777844051105d0,'dip_vars(): cdcd')
call test_value(cdsd,0.32139380484326974d0,'dip_vars(): cdsd')
write(6,*) 'subroutine dip_vars() passed unit test'
write(6,*)

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
write(6,*) 'subroutine pt_src_vars() passed unit test'
write(6,*)

!----
! Point source example from Okada (1985)
!----
moment = 0.0d0
moment(1) = 40.0d0
call o92_pt_disp(disp,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(disp(1),-9.4474186064630631d-4,'o92_pt_disp(): disp(1) (strike-slip)')
call test_value(disp(2),-1.0229885427654646d-3,'o92_pt_disp(): disp(2) (strike-slip)')
call test_value(disp(3),-7.4201246532591128d-4,'o92_pt_disp(): disp(3) (strike-slip)')
write(6,*) 'subroutine o92_pt_disp() passed strike-slip unit test'
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
write(6,*) 'subroutine o92_pt_partials() passed strike-slip unit test'
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
write(6,*) 'subroutine o92_pt_strain() passed strike-slip unit test'
moment = 0.0d0
moment(2) = 40.0d0
call o92_pt_disp(disp,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(disp(1),-1.1723215364075840d-3,'o92_pt_disp(): disp(1) (dip-slip)')
call test_value(disp(2),-2.0819985410726201d-3,'o92_pt_disp(): disp(2) (dip-slip)')
call test_value(disp(3),-2.5315947108865925d-3,'o92_pt_disp(): disp(3) (dip-slip)')
write(6,*) 'subroutine o92_pt_disp() passed dip-slip unit test'
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
write(6,*) 'subroutine o92_pt_partials() passed dip-slip unit test'
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
write(6,*) 'subroutine o92_pt_strain() passed dip-slip unit test'
moment = 0.0d0
moment(3) = 40.0d0
call o92_pt_disp(disp,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(disp(1),-3.5716333691232548d-4,'o92_pt_disp(): disp(1) (tensile-slip)')
call test_value(disp(2), 3.5310854913572190d-4,'o92_pt_disp(): disp(2) (tensile-slip)')
call test_value(disp(3),-2.0068126969335952d-4,'o92_pt_disp(): disp(3) (tensile-slip)')
write(6,*) 'subroutine o92_pt_disp() passed tensile-slip unit test'
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
write(6,*) 'subroutine o92_pt_partials() passed tensile-slip unit test'
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
write(6,*) 'subroutine o92_pt_strain() passed tensile-slip unit test'
write(6,*) 'point source example from Okada (1985) passed unit test'
write(6,*)

!----
! Rectangular fault example from Okada (1985)
!----
len = 3.0d0
wid = 2.0d0
call rect_src_coords(sta_coord(1),sta_coord(2),sta_coord(3),evdp,wid,len)
call test_value(p,4.7848309131206399d0,'rect_src_coords(): p')
call test_value(q,1.4509972890550495d0,'rect_src_coords(): q')
call test_value(eta_vec(1),3.7848309131206399d0,'rect_src_coords(): eta_vec(1)')
call test_value(eta_vec(2),5.7848309131206399d0,'rect_src_coords(): eta_vec(2)')
call test_value(ksi_vec(1),0.5000000000000000d0,'rect_src_coords(): ksi_vec(1)')
call test_value(ksi_vec(2),3.5000000000000000d0,'rect_src_coords(): ksi_vec(2)')
write(6,*) 'subroutine rect_src_coords() passed unit test'
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
write(6,*) 'subroutine rect_src_vars() passed unit test'
write(6,*)

! Original Okada derivation have fault origin at its bottom left corner, but I use the fault
! center as the coordinates. To compare with Okada (1985), shift fault.
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
write(6,*) 'subroutine o92_rect_disp() passed strike-slip unit test'
! Dip-slip
slip_vec(1) = 0.0d0
slip_vec(2) = 1.0d0
slip_vec(3) = 0.0d0
call o92_rect_disp(disp,sta_coord,evdp,dip,slip_vec,wid,len,lambda,shear_modulus)
call test_value(disp(1),-4.6823487628354045d-003,'o92_rect_disp(): disp(1) (dip-slip)')
call test_value(disp(2),-3.5267267968718020d-002,'o92_rect_disp(): disp(2) (dip-slip)')
call test_value(disp(3),-3.5638557673268470d-002,'o92_rect_disp(): disp(3) (dip-slip)')
write(6,*) 'subroutine o92_rect_disp() passed dip-slip unit test'
! Tensile-slip
slip_vec(1) = 0.0d0
slip_vec(2) = 0.0d0
slip_vec(3) = 1.0d0
call o92_rect_disp(disp,sta_coord,evdp,dip,slip_vec,wid,len,lambda,shear_modulus)
call test_value(disp(1),-2.6599600964428191d-004,'o92_rect_disp(): disp(1) (tensile-slip)')
call test_value(disp(2), 1.0564074876983659d-002,'o92_rect_disp(): disp(2) (tensile-slip)')
call test_value(disp(3), 3.2141931142207676d-003,'o92_rect_disp(): disp(3) (tensile-slip)')
write(6,*) 'subroutine o92_rect_disp() passed tensile-slip unit test'
write(6,*) 'rectangular source example from Okada (1985) passed unit test'
write(6,*)

!----
! Example with point source, station depth greater than zero, oblique slip
!----
sta_coord(1) = -12.0d0
sta_coord(2) = 6.2d0
sta_coord(3) = 13.0d0
evdp = 7.0d0
moment(1) = -27.0d0
moment(2) = 13.0d0
moment(3) = 0.0d0
call o92_pt_disp(disp,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(disp(1), 1.9245175992850267d-004,'o92_pt_disp(): disp(1) (z>0; oblique-slip)')
call test_value(disp(2),-1.9098044959056914d-004,'o92_pt_disp(): disp(2) (z>0; oblique-slip)')
call test_value(disp(3), 8.1576627770328566d-005,'o92_pt_disp(): disp(3) (z>0; oblique-slip)')
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
call o92_pt_strain(strain,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(strain(1,1), 2.0874574675642794d-5,'o92_pt_strain(): strain(1,1) (z>0; oblique-slip)')
call test_value(strain(1,2),-9.6052553102692557d-6,'o92_pt_strain(): strain(1,2) (z>0; oblique-slip)')
call test_value(strain(1,3), 1.7465760943215932d-5,'o92_pt_strain(): strain(1,3) (z>0; oblique-slip)')
call test_value(strain(2,1),-9.6052553102692557d-6,'o92_pt_strain(): strain(2,1) (z>0; oblique-slip)')
call test_value(strain(2,2),-6.8032602399881741d-6,'o92_pt_strain(): strain(2,2) (z>0; oblique-slip)')
call test_value(strain(2,3),-1.4000925908121080d-5,'o92_pt_strain(): strain(2,3) (z>0; oblique-slip)')
call test_value(strain(3,1), 1.7465760943215932d-5,'o92_pt_strain(): strain(3,1) (z>0; oblique-slip)')
call test_value(strain(3,2),-1.4000925908121080d-5,'o92_pt_strain(): strain(3,2) (z>0; oblique-slip)')
call test_value(strain(3,3),-1.7420298917986620d-6,'o92_pt_strain(): strain(3,3) (z>0; oblique-slip)')
write(6,*) 'example with depth>0, oblique slip passed unit test'
write(6,*)

!----
! Example with rectangular source, station depth greater than zero
!----
sta_coord(1) = -12.0d0
sta_coord(2) = 6.2d0
sta_coord(3) = 13.0d0
evdp = 7.0d0
wid=3.2d0
len=6.5d0

! Strike-slip
slip_vec(1) = 1.0d0
slip_vec(2) = 0.0d0
slip_vec(3) = 0.0d0
! Test components of displacement
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
! Test total displacement
call o92_rect_disp(disp,sta_coord,evdp,dip,slip_vec,wid,len,lambda,shear_modulus)
call test_value(disp(1),-6.6858712769214453d-003,'o92_rect_disp(): disp(1) (z>0; strike-slip)')
call test_value(disp(2), 5.6833739558216096d-003,'o92_rect_disp(): disp(2) (z>0; strike-slip)')
call test_value(disp(3),-3.8337220569204383d-003,'o92_rect_disp(): disp(3) (z>0; strike-slip)')
write(6,*)

! Dip-slip
slip_vec(1) = 0.0d0
slip_vec(2) = 1.0d0
slip_vec(3) = 0.0d0
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
! Test total displacement
call o92_rect_disp(disp,sta_coord,evdp,dip,slip_vec,wid,len,lambda,shear_modulus)
call test_value(disp(1),-1.4891278409154987d-003,'o92_rect_disp(): disp(1) (z>0; dip-slip)')
call test_value(disp(2),-6.6628749834942423d-004,'o92_rect_disp(): disp(2) (z>0; dip-slip)')
call test_value(disp(3),-2.6808905339339140d-003,'o92_rect_disp(): disp(3) (z>0; dip-slip)')
write(6,*)
! write(6,*) disp

xin = sta_coord(1)
yin = sta_coord(2)
stdp = sta_coord(3)
dipin = dip
! rakin = datan2(moment(2),moment(1))/d2r
! area = 1.0d0
! slip = dsqrt(moment(1)*moment(1)+moment(2)*moment(2))/shear_modulus/area
dens = 1.0d0
vp = dsqrt((lambda+2.0d0*shear_modulus)/dens)
vs = dsqrt(shear_modulus/dens)
! call o92pt(ux,uy,uz,xin,yin,stdp,evdp,dipin,rakin,area,slip,vp,vs,dens)
! write(6,*) ux,uy,uz
! call o92ptstn(strain,xin,yin,stdp,evdp,dipin,rakin,area,slip,vp,vs,dens)
! write(6,*) strain(1,:)
! write(6,*) strain(2,:)
! write(6,*) strain(3,:)

rakin = datan2(slip_vec(2),slip_vec(1))/d2r
slip = dsqrt(slip_vec(1)*slip_vec(1)+slip_vec(2)*slip_vec(2))
call o92rect(ux,uy,uz,xin,yin,stdp,evdp,dipin,rakin,wid,len,slip,vp,vs,dens)
! write(6,*) ux,uy,uz

write(6,*) 'okada92_module unit test passed'
end

!--------------------------------------------------------------------------------------------------!

subroutine test_value(actual_value,expected_value,string)
implicit none
double precision :: actual_value, expected_value, diff
character(len=*) :: string

diff = (actual_value-expected_value)/actual_value

if (dabs(diff).gt.1.0d-10) then
    write(0,*) 'okada92 unit test FAILED for '//trim(string)
    write(0,*) '    expected value: ',expected_value
    write(0,*) '    computed value: ',actual_value
    stop
else
    write(6,*) 'okada92 unit test passed for '//trim(string)
endif
return
end
