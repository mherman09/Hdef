program main

use test, only: test_value
use okada92_module

implicit none

double precision :: disp(3), strain(3,3), sta_coord(3), evdp, dip, moment(4), lambda, shear_modulus
double precision :: slip_vec(3), wid, len
double precision :: ux,uy,uz,xin,yin,stdp,dipin,rakin,slip,vp,vs,dens
double precision :: uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz

! Half-space, dip, and geometric variables are private to module; use checking subroutines built
! into module for these
call test_halfspace_vars()
call test_dip_vars()
call test_pt_src_vars()
call test_rect_src_vars()

!----
! Point source
!----
! Point source example from Okada (1985)
lambda = 40.0d0            ! lambda = shear_modulus implies
shear_modulus = 40.0d0     !     Poisson's ratio of 0.25
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

! Dip-slip
moment(1) =  0.0d0
moment(2) = 40.0d0
moment(3) =  0.0d0
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

! Tensile-slip
moment(1) =  0.0d0
moment(2) =  0.0d0
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
! Rectangular source
!----
! Example from Okada (1985)
lambda = 40.0d0            ! lambda = shear_modulus implies
shear_modulus = 40.0d0     !     Poisson's ratio of 0.25
dip = 70.0d0
sta_coord(1) = 2.0d0
sta_coord(2) = 3.0d0
sta_coord(3) = 0.0d0
evdp = 4.0d0
len = 3.0d0
wid = 2.0d0

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

! Example with rectangular source, station depth greater than zero
lambda = 40.0d0
shear_modulus = 40.0d0
dip = 70.0d0
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

! Time to test rectangular fault STRAINS! WOOOOOO!
! Strike-slip
slip_vec(1) = 1.0d0
slip_vec(2) = 0.0d0
slip_vec(3) = 0.0d0
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
write(0,*) disp
! call test_value(disp(1),2.8890966553217295d-002,'duAdy_ss_rect(): (1)')
! call test_value(disp(2),1.2845371264737372d-002,'duAdy_ss_rect(): (2)')
! call test_value(disp(3),2.7617711810937568d-003,'duAdy_ss_rect(): (3)')
disp = duAdz_ds_rect()
! call test_value(disp(1), 1.2845371264737372d-002,'duAdy_ds_rect(): (1)')
! call test_value(disp(2), 6.0051670894042106d-002,'duAdy_ds_rect(): (2)')
! call test_value(disp(3), 6.8785272909704029d-003,'duAdy_ds_rect(): (3)')
! disp = duAdy_ts_rect()
disp = duBdz_ss_rect()
! call test_value(disp(1),-5.4892739626334362d-002,'duBdy_ss_rect(): (1)')
! call test_value(disp(2),-2.7969238876413290d-002,'duBdy_ss_rect(): (2)')
! call test_value(disp(3),-3.1020145458715599d-003,'duBdy_ss_rect(): (3)')
disp = duBdz_ds_rect()
! call test_value(disp(1),-3.8965904427687062d-002,'duBdy_ds_rect(): (1)')
! call test_value(disp(2),-0.14704676790658636d0,'duBdy_ds_rect(): (2)')
! call test_value(disp(3),-9.1052106567444321d-003,'duBdy_ds_rect(): (3)')
! disp = duBdy_ts_rect()
disp = duCdz_ss_rect()
! call test_value(disp(1),3.6229956319952023d-004,'duCdy_ss_rect(): (1)')
! call test_value(disp(2),2.7645676601602590d-004,'duCdy_ss_rect(): (2)')
! call test_value(disp(3),5.2251433363535108d-004,'duCdy_ss_rect(): (3)')
disp = duCdz_ds_rect()
! call test_value(disp(1),-1.1160386819908798d-003,'duCdy_ds_rect(): (1)')
! call test_value(disp(2),-6.3814177065406037d-004,'duCdy_ds_rect(): (2)')
! call test_value(disp(3), 2.7583611996049840d-003,'duCdy_ds_rect(): (3)')
disp = duCdz_ts_rect()
write(0,*)

call o92_rect_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz, &
                       sta_coord,evdp,dip,slip_vec,wid,len,lambda,shear_modulus)
write(0,*) uxx,uyy,uzz

!----
! Compare with original code results
!----
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
call o92rectstn(strain,xin,yin,stdp,evdp,dipin,rakin,wid,len,slip,vp,vs,dens)
write(0,*) strain(1,1),strain(2,2),strain(3,3)

write(6,*) 'okada92_module unit test passed'
end
