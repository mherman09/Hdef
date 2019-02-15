program main
use okada92_module
implicit none
double precision :: disp(3), strain(3,3), sta_coord(3), evdp, dip, moment(4), lambda, shear_modulus
double precision :: slip_vec(3), wid, len
double precision :: ux,uy,uz,xin,yin,stdp,dipin,rakin,area,slip,vp,vs,dens
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
write(0,*) 'subroutine halfspace_vars() passed unit test'
write(0,*)

dip = 70.0d0
call dip_vars(dip)
call test_value(sd,0.93969262078590832d0,'dip_vars(): sd')
call test_value(cd,0.34202014332566882d0,'dip_vars(): cd')
call test_value(s2d,0.64278760968653947d0,'dip_vars(): s2d')
call test_value(c2d,-0.76604444311897790d0,'dip_vars(): c2d')
call test_value(sdsd,0.88302222155948884d0,'dip_vars(): sdsd')
call test_value(cdcd,0.11697777844051105d0,'dip_vars(): cdcd')
call test_value(cdsd,0.32139380484326974d0,'dip_vars(): cdsd')
write(0,*) 'subroutine dip_vars() passed unit test'
write(0,*)

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
write(0,*) 'subroutine pt_src_vars() passed unit test'
write(0,*)

moment = 0.0d0
moment(1) = 40.0d0
call o92_pt_disp(disp,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(disp(1),-9.4474186064630631d-4,'o92_pt_disp(): disp(1)')
call test_value(disp(2),-1.0229885427654646d-3,'o92_pt_disp(): disp(2)')
call test_value(disp(3),-7.4201246532591128d-4,'o92_pt_disp(): disp(3)')
write(0,*) 'subroutine o92_pt_disp() passed strike-slip unit test'
call o92_pt_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz,&
                     sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(uxx,-2.2862094355883986d-004,'o92_pt_partials(): uxx')
call test_value(uxy,-1.4246398245262621d-004,'o92_pt_partials(): uxy')
call test_value(uxz,6.2590905699873601d-005,'o92_pt_partials(): uxz')
call test_value(uyx,-2.0511017314976633d-004,'o92_pt_partials(): uyx')
call test_value(uyy,-3.0067825210386353d-004,'o92_pt_partials(): uyy')
call test_value(uyz,1.6933001353308135d-004,'o92_pt_partials(): uyz')
call test_value(uzx,-6.2590905699873859d-005,'o92_pt_partials(): uzx')
call test_value(uzy,-1.6933001353308129d-004,'o92_pt_partials(): uzy')
call test_value(uzz,1.7643306522090115d-004,'o92_pt_partials(): uzz')
write(0,*) 'subroutine o92_pt_partials() passed strike-slip unit test'
call o92_pt_strain(strain,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(strain(1,1),-2.2862094355883986d-004,'o92_pt_strain(): strain(1,1)')
call test_value(strain(1,2),-1.7378707780119627d-004,'o92_pt_strain(): strain(1,2)')
call test_value(strain(1,3),-1.2874900798265365d-019,'o92_pt_strain(): strain(1,3)')
call test_value(strain(2,1),-1.7378707780119627d-004,'o92_pt_strain(): strain(2,1)')
call test_value(strain(2,2),-3.0067825210386353d-004,'o92_pt_strain(): strain(2,2)')
call test_value(strain(2,3), 2.7105054312137611d-020,'o92_pt_strain(): strain(2,3)')
call test_value(strain(3,1),-1.2874900798265365d-019,'o92_pt_strain(): strain(3,1)')
call test_value(strain(3,2), 2.7105054312137611d-020,'o92_pt_strain(): strain(3,2)')
call test_value(strain(3,3), 1.7643306522090115d-004,'o92_pt_strain(): strain(3,3)')
write(0,*) 'subroutine o92_pt_strain() passed strike-slip unit test'
moment = 0.0d0
moment(2) = 40.0d0
call o92_pt_disp(disp,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(disp(1),-1.1723215364075840d-3,'o92_pt_disp(): disp(1)')
call test_value(disp(2),-2.0819985410726201d-3,'o92_pt_disp(): disp(2)')
call test_value(disp(3),-2.5315947108865925d-3,'o92_pt_disp(): disp(3)')
write(0,*) 'subroutine o92_pt_disp() passed dip-slip unit test'
call o92_pt_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz,&
                     sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(uxx,-1.5259408689661085d-004,'o92_pt_partials(): uxx')
call test_value(uxy,-3.5441800891339358d-004,'o92_pt_partials(): uxy')
call test_value(uxz,-8.7070529207727333d-004,'o92_pt_partials(): uxz')
call test_value(uyx, 6.9826591999645748d-004,'o92_pt_partials(): uyx')
call test_value(uyy,-1.1537526800074350d-003,'o92_pt_partials(): uyy')
call test_value(uyz, 6.3453607666584869d-004,'o92_pt_partials(): uyz')
call test_value(uzx, 8.7070529207727322d-004,'o92_pt_partials(): uzx')
call test_value(uzy,-6.3453607666584847d-004,'o92_pt_partials(): uzy')
call test_value(uzz, 4.3544892230134886d-004,'o92_pt_partials(): uzz')
write(0,*) 'subroutine o92_pt_partials() passed dip-slip unit test'
call o92_pt_strain(strain,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(strain(1,1),-1.5259408689661085d-004,'o92_pt_strain(): strain(1,1)')
call test_value(strain(1,2), 1.7192395554153195d-004,'o92_pt_strain(): strain(1,2)')
call test_value(strain(1,3),-5.4210108624275222d-020,'o92_pt_strain(): strain(1,3)')
call test_value(strain(2,1), 1.7192395554153195d-004,'o92_pt_strain(): strain(2,1)')
call test_value(strain(2,2),-1.1537526800074350d-003,'o92_pt_strain(): strain(2,2)')
call test_value(strain(2,3), 1.0842021724855044d-019,'o92_pt_strain(): strain(2,3)')
call test_value(strain(3,1),-5.4210108624275222d-020,'o92_pt_strain(): strain(3,1)')
call test_value(strain(3,2), 1.0842021724855044d-019,'o92_pt_strain(): strain(3,2)')
call test_value(strain(3,3), 4.3544892230134886d-004,'o92_pt_strain(): strain(3,3)')
write(0,*) 'subroutine o92_pt_strain() passed dip-slip unit test'
moment = 0.0d0
moment(3) = 40.0d0
call o92_pt_disp(disp,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(disp(1),-3.5716333691232548d-4,'o92_pt_disp(): disp(1)')
call test_value(disp(2), 3.5310854913572190d-4,'o92_pt_disp(): disp(2)')
call test_value(disp(3),-2.0068126969335952d-4,'o92_pt_disp(): disp(3)')
write(0,*) 'subroutine o92_pt_disp() passed tensile-slip unit test'
call o92_pt_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz,&
                     sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(uxx,-1.3597209844017221d-004,'o92_pt_partials(): uxx')
call test_value(uxy, 5.0725453276876352d-004,'o92_pt_partials(): uxy')
call test_value(uxz,-7.5405344488748576d-005,'o92_pt_partials(): uxz')
call test_value(uyx,-6.7733492861951009d-005,'o92_pt_partials(): uyx')
call test_value(uyy, 6.8111287703614769d-004,'o92_pt_partials(): uyy')
call test_value(uyz,-8.1037165493104131d-004,'o92_pt_partials(): uyz')
call test_value(uzx, 7.5405344488748590d-005,'o92_pt_partials(): uzx')
call test_value(uzy, 8.1037165493104120d-004,'o92_pt_partials(): uzy')
call test_value(uzz,-1.8171359286532499d-004,'o92_pt_partials(): uzz')
write(0,*) 'subroutine o92_pt_partials() passed tensile-slip unit test'
call o92_pt_strain(strain,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(strain(1,1),-1.3597209844017221d-004,'o92_pt_strain(): strain(1,1)')
call test_value(strain(1,2), 2.1976051995340625d-004,'o92_pt_strain(): strain(1,2)')
call test_value(strain(1,3), 6.7762635780344027d-021,'o92_pt_strain(): strain(1,3)')
call test_value(strain(2,1), 2.1976051995340625d-004,'o92_pt_strain(): strain(2,1)')
call test_value(strain(2,2), 6.8111287703614769d-004,'o92_pt_strain(): strain(2,2)')
call test_value(strain(2,3),-5.4210108624275222d-020,'o92_pt_strain(): strain(2,3)')
call test_value(strain(3,1), 6.7762635780344027d-021,'o92_pt_strain(): strain(3,1)')
call test_value(strain(3,2),-5.4210108624275222d-020,'o92_pt_strain(): strain(3,2)')
call test_value(strain(3,3),-1.8171359286532499d-004,'o92_pt_strain(): strain(3,3)')
write(0,*) 'subroutine o92_pt_strain() passed tensile-slip unit test'
write(0,*) 'point source example from Okada (1985) passed unit test'
write(0,*)

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
write(0,*) 'subroutine rect_src_coords() passed unit test'
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
write(0,*) 'R:',R
write(0,*) 'R2:',R2
write(0,*) 'R3:',R3
write(0,*) 'R5:',R5
write(0,*) 'Rd:',Rd
write(0,*) 'qq:',qq
! TH = datan(ksi*eta/(q*R))
! Rk = R + ksi
! X11 = 1.0d0/(R*Rk)
! X32 = (2.0d0*R+ksi)/(R3*Rk*Rk)
! X53 = (8.0d0*R2+9.0d0*R*ksi+3.0d0*ksi*ksi)/(R5*Rk*Rk*Rk)
! logRk = dlog(R+ksi)
! Re = R+eta
! Y11 = 1.0d0/(R*Re)
! Y32 = (2.0d0*R+eta)/(R3*Re*Re)
! Y53 = (8.0d0*R2+9.0d0*R*eta+3.0d0*eta*eta)/(R5*Re*Re*Re)
! logRe = dlog(R+eta)
! Y0  = Y11 - ksi*ksi*Y32
! I4 = sd*ksi/(Rd*cd) + 2.0d0*datan((eta*(xxx+q*cd)+xxx*(R+xxx)*sd)/(ksi*(R+xxx)*cd))/(cdcd)
! I3 = ybar/(cd*Rd) - (logRe-sd*dlog(Rd))/(cdcd)
! I1 = -ksi*cd/Rd - I4*sd
! I2 = dlog(Rd) + I3*sd
! Z32 = sd/R3 - h*Y32
! Z53 = 3.0d0*sd/R5 - h*Y53
! Z0  = Z32 - ksi*ksi*Z53
! D11 = 1.0d0/(R*Rd)
! K1 =ksi*(D11-Y11*sd)/cd
! K3 = (q*Y11-ybar*D11)/cd
! K2 = 1.0d0/R+K3*sd
! K4 = ksi*Y11*cd-K1*sd
! J2 = ksi*ybar*D11/Rd
! J5 = -(dbar+ybar*ybar/Rd)*D11
! J3 = (K1-J2*sd)/cd
! J6 = (K3-J5*sd)/cd
! J1 = J5*cd - J6*sd
! J4 = -ksi*Y11 - J2*cd + J3*sd
! E2 = sd/R - ybar*q/R3
! F2 = dbar/R3 + ksi*ksi*Y32*sd
! G2 = 2.0d0*X11*sd - ybar*q*X32
! H2 = dbar*q*X32 + ksi*q*Y32*sd
! P2 = cd/R3 + q*Y32*sd
! Q2 = 3.0d0*cbar*dbar/R5 - (z*Y32+Z32+Z0)*sd
! E3 = cd/R + dbar*q/R3
! F3 = ybar/R3 + ksi*ksi*Y32*cd
! G3 = 2.0d0*X11*cd + dbar*q*X32
! H3 = ybar*q*X32 + ksi*q*Y32*cd
! P3 = sd/R3 - q*Y32*cd
! Q3 = 3.0d0*cbar*ybar/R5 + q*Y32 - (z*Y32+Z32+Z0)*cd
write(0,*) 'subroutine rect_src_vars() passed unit test'

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
call test_value(disp(1),-8.6891650042561549d-003,'o92_rect_disp(): disp(1)')
call test_value(disp(2),-4.2975821897415613d-003,'o92_rect_disp(): disp(2)')
call test_value(disp(3),-2.7474058276388967d-003,'o92_rect_disp(): disp(3)')
write(0,*) 'subroutine o92_rect_disp() passed strike-slip unit test'
! Dip-slip
slip_vec(1) = 0.0d0
slip_vec(2) = 1.0d0
slip_vec(3) = 0.0d0
call o92_rect_disp(disp,sta_coord,evdp,dip,slip_vec,wid,len,lambda,shear_modulus)
call test_value(disp(1),-4.6823487628354045d-003,'o92_rect_disp(): disp(1)')
call test_value(disp(2),-3.5267267968718020d-002,'o92_rect_disp(): disp(2)')
call test_value(disp(3),-3.5638557673268470d-002,'o92_rect_disp(): disp(3)')
write(0,*) 'subroutine o92_rect_disp() passed dip-slip unit test'
! Tensile-slip
slip_vec(1) = 0.0d0
slip_vec(2) = 0.0d0
slip_vec(3) = 1.0d0
call o92_rect_disp(disp,sta_coord,evdp,dip,slip_vec,wid,len,lambda,shear_modulus)
write(0,*) 'disp(1):',disp(1)
write(0,*) 'disp(2):',disp(2)
write(0,*) 'disp(3):',disp(3)
call test_value(disp(1),-2.6599600964428191d-004,'o92_rect_disp(): disp(1)')
call test_value(disp(2), 1.0564074876983659d-002,'o92_rect_disp(): disp(2)')
call test_value(disp(3), 3.2141931142207676d-003,'o92_rect_disp(): disp(3)')
write(0,*) 'subroutine o92_rect_disp() passed tensile-slip unit test'
write(0,*) 'rectangular source example from Okada (1985) passed unit test'
write(0,*)

! Example with depth greater than zero, oblique slip
sta_coord(1) = -12.0d0
sta_coord(2) = 6.2d0
sta_coord(3) = 13.0d0
evdp = 7.0d0
moment(1) = -27.0d0
moment(2) = 13.0d0
moment(3) = 0.0d0
call o92_pt_disp(disp,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(disp(1), 1.9245175992850267d-004,'o92_pt_disp(): disp(1)')
call test_value(disp(2),-1.9098044959056914d-004,'o92_pt_disp(): disp(2)')
call test_value(disp(3), 8.1576627770328566d-005,'o92_pt_disp(): disp(3)')
call o92_pt_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz,&
                     sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(uxx, 2.0874574675642794d-005,'o92_pt_partials(): uxx')
call test_value(uxy, 2.3620493327193696d-006,'o92_pt_partials(): uxy')
call test_value(uxz, 2.2899507819320573d-005,'o92_pt_partials(): uxz')
call test_value(uyx,-2.1572559953257881d-005,'o92_pt_partials(): uyx')
call test_value(uyy,-6.8032602399881741d-006,'o92_pt_partials(): uyy')
call test_value(uyz,-2.1249535452408235d-005,'o92_pt_partials(): uyz')
call test_value(uzx, 1.2032014067111290d-005,'o92_pt_partials(): uzx')
call test_value(uzy,-6.7523163638339241d-006,'o92_pt_partials(): uzy')
call test_value(uzz,-1.7420298917986620d-006,'o92_pt_partials(): uzz')
call o92_pt_strain(strain,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value(strain(1,1), 2.0874574675642794d-005,'o92_pt_strain(): strain(1,1)')
call test_value(strain(1,2),-9.6052553102692557d-006,'o92_pt_strain(): strain(1,2)')
call test_value(strain(1,3), 1.7465760943215932d-005,'o92_pt_strain(): strain(1,3)')
call test_value(strain(2,1),-9.6052553102692557d-006,'o92_pt_strain(): strain(2,1)')
call test_value(strain(2,2),-6.8032602399881741d-006,'o92_pt_strain(): strain(2,2)')
call test_value(strain(2,3),-1.4000925908121080d-005,'o92_pt_strain(): strain(2,3)')
call test_value(strain(3,1), 1.7465760943215932d-005,'o92_pt_strain(): strain(3,1)')
call test_value(strain(3,2),-1.4000925908121080d-005,'o92_pt_strain(): strain(3,2)')
call test_value(strain(3,3),-1.7420298917986620d-006,'o92_pt_strain(): strain(3,3)')
write(0,*) 'example with depth>0, oblique slip passed unit test'

xin = sta_coord(1)
yin = sta_coord(2)
stdp = sta_coord(3)
dipin = dip
rakin = datan2(moment(2),moment(1))/d2r
area = 1.0d0
slip = dsqrt(moment(1)*moment(1)+moment(2)*moment(2))/shear_modulus/area
dens = 1.0d0
vp = dsqrt((lambda+2.0d0*shear_modulus)/dens)
vs = dsqrt(shear_modulus/dens)
! call o92pt(ux,uy,uz,xin,yin,stdp,evdp,dipin,rakin,area,slip,vp,vs,dens)
! write(0,*) ux,uy,uz
call o92ptstn(strain,xin,yin,stdp,evdp,dipin,rakin,area,slip,vp,vs,dens)
write(0,*) strain(1,:)
write(0,*) strain(2,:)
write(0,*) strain(3,:)

write(0,*) 'okada92_module unit test passed'
end

!--------------------------------------------------------------------------------------------------!

subroutine test_value(actual_value,expected_value,string)
implicit none
double precision :: actual_value, expected_value, diff
character(len=*) :: string

diff = (actual_value-expected_value)/actual_value

if (dabs(diff).gt.1.0d-10) then
    write(0,*) 'tri_disloc unit test FAILED for '//trim(string)
    write(0,*) '    expected value: ',expected_value
    write(0,*) '    computed value: ',actual_value
    stop
else
    write(0,*) 'tri_disloc unit test passed for '//trim(string)
endif
return
end
