program main
use okada92_module
implicit none
double precision :: disp(3), strain(3,3), sta_coord(3), evdp, dip, moment(4), lambda, shear_modulus
double precision :: ux,uy,uz,xin,yin,stdp,dipin,rakin,area,slip,vp,vs,dens
double precision :: uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz
! Example from Okada (1985) for surface deformation
lambda = 40.0d0
shear_modulus = 40.0d0
call halfspace_vars(lambda,shear_modulus)
call test_value((a-0.66666666666666663d0)/a ,'halfspace_vars(): a')
call test_value((CA1-0.16666666666666669d0)/CA1 ,'halfspace_vars(): CA1')
call test_value((CA2-0.33333333333333331d0)/CA2 ,'halfspace_vars(): CA2')
call test_value((CB-0.50000000000000011d0)/CB ,'halfspace_vars(): CB')
call test_value((CC-0.33333333333333337d0)/CC ,'halfspace_vars(): CC')
write(0,*) 'subroutine halfspace_vars() passed unit test'
write(0,*)

dip = 70.0d0
call dip_vars(dip)
call test_value((sd-0.93969262078590832d0)/sd ,'dip_vars(): sd')
call test_value((cd-0.34202014332566882d0)/cd ,'dip_vars(): cd')
call test_value((s2d-0.64278760968653947d0)/s2d ,'dip_vars(): s2d')
call test_value((c2d+0.76604444311897790d0)/c2d ,'dip_vars(): c2d')
call test_value((sdsd-0.88302222155948884d0)/sdsd ,'dip_vars(): sdsd')
call test_value((cdcd-0.11697777844051105d0)/cdcd ,'dip_vars(): cdcd')
call test_value((cdsd-0.32139380484326974d0)/cdsd ,'dip_vars(): cdsd')
write(0,*) 'subroutine dip_vars() passed unit test'
write(0,*)

sta_coord(1) = 2.0d0
sta_coord(2) = 3.0d0
sta_coord(3) = 0.0d0
evdp = 4.0d0
call pt_src_vars(sta_coord(1),sta_coord(2),sta_coord(3),evdp)
call test_value((p-4.7848309131206399d0)/p ,'pt_src_vars(): p')
call test_value((q-1.4509972890550495d0)/q ,'pt_src_vars(): q')
call test_value((s-4.9925406015355298d0)/s ,'pt_src_vars(): s')
call test_value((t-0.27301710938922441d0)/t ,'pt_src_vars(): t')
call test_value((xx-4.0d0)/xx ,'pt_src_vars(): xx')
call test_value((xy-6.0d0)/xy ,'pt_src_vars(): xy')
call test_value((yy-9.0d0)/yy ,'pt_src_vars(): yy')
call test_value((dd-16.0d0)/dd ,'pt_src_vars(): dd')
call test_value((xd-8.0d0)/xd ,'pt_src_vars(): xd')
call test_value((xc-8.0d0)/xc ,'pt_src_vars(): xc')
call test_value((xq-2.9019945781100991d0)/xq ,'pt_src_vars(): xq')
call test_value((yq-4.3529918671651489d0)/yq ,'pt_src_vars(): yq')
call test_value((dq-5.8039891562201982d0)/dq ,'pt_src_vars(): dq')
call test_value((pq-6.9427766835248459d0)/pq ,'pt_src_vars(): pq')
call test_value((R-5.3851648071345037d0)/R ,'pt_src_vars(): R')
call test_value((R2-28.999999999999996d0)/R2 ,'pt_src_vars(): R2')
call test_value((R3-156.16977940690060d0)/R3 ,'pt_src_vars(): R3')
call test_value((R4-840.99999999999989d0)/R4 ,'pt_src_vars(): R4')
call test_value((R5-4528.9236028001169d0)/R5 ,'pt_src_vars(): R5')
call test_value((R7-131338.78448120339d0)/R7 ,'pt_src_vars(): R7')
call test_value((Rd-9.3851648071345046d0)/Rd ,'pt_src_vars(): Rd')
call test_value((A3-0.58620689655172409d0)/A3 ,'pt_src_vars(): A3')
call test_value((A5-0.31034482758620685d0)/A5 ,'pt_src_vars(): A5')
call test_value((A7-3.4482758620689502d-2)/A7 ,'pt_src_vars(): A7')
call test_value((B3-6.8965517241379226d-2)/B3 ,'pt_src_vars(): B3')
call test_value((B5+0.55172413793103470d0)/B5 ,'pt_src_vars(): B5')
call test_value((B7+1.1724137931034484d0)/B7 ,'pt_src_vars(): B7')
call test_value((C3+0.65517241379310365d0)/C3 ,'pt_src_vars(): C3')
call test_value((C5+1.7586206896551726d0)/C5 ,'pt_src_vars(): C5')
call test_value((C7+2.8620689655172420d0)/C7 ,'pt_src_vars(): C7')
call test_value((I1-4.4511857695225003d-003)/I1 ,'pt_src_vars(): I1')
call test_value((I2-1.4062132997373915d-003)/I2 ,'pt_src_vars(): I2')
call test_value((I3-1.1400361746955842d-002)/I3 ,'pt_src_vars(): I3')
call test_value((I4+6.4425920722921085d-003)/I4 ,'pt_src_vars(): I4')
call test_value((I5-1.5490988410148277d-002)/I5 ,'pt_src_vars(): I5')
call test_value((J1+2.3037701706801377d-003)/J1 ,'pt_src_vars(): J1')
call test_value((J2-8.3831231772032186d-004)/J2 ,'pt_src_vars(): J2')
call test_value((J3-2.9153389890690740d-003)/J3 ,'pt_src_vars(): J3')
call test_value((J4+1.6706841541556938d-003)/J4 ,'pt_src_vars(): J4')
call test_value((K1+1.7024420377114214d-003)/K1 ,'pt_src_vars(): K1')
call test_value((K2-1.3075030688791331d-004)/K2 ,'pt_src_vars(): K2')
call test_value((K3+5.4300227400023550d-003)/K3 ,'pt_src_vars(): K3')
call test_value((U2-0.18917678161950324d0)/U2 ,'pt_src_vars(): U2')
call test_value((V2-1.4014492135054364d0)/V2 ,'pt_src_vars(): V2')
call test_value((W2-1.1288694024054116d0)/W2 ,'pt_src_vars(): W2')
call test_value((U3-1.3427079288808756d0)/U3 ,'pt_src_vars(): U3')
call test_value((V3-5.0611389600960157d0)/V3 ,'pt_src_vars(): V3')
call test_value((W3-1.6847280722065445d0)/W3 ,'pt_src_vars(): W3')
write(0,*) 'subroutine pt_src_vars() passed unit test'
write(0,*)

moment = 0.0d0
moment(1) = 40.0d0
call o92_pt_disp(disp,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value((disp(1)+9.4474186064630631d-4)/disp(1) ,'o92_pt_disp(): disp(1)')
call test_value((disp(2)+1.0229885427654646d-3)/disp(2) ,'o92_pt_disp(): disp(2)')
call test_value((disp(3)+7.4201246532591128d-4)/disp(3) ,'o92_pt_disp(): disp(3)')
write(0,*) 'subroutine o92_pt_disp() passed strike-slip unit test'
call o92_pt_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz,&
                     sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value((uxx+2.2862094355883986d-004)/uxx ,'o92_pt_partials(): uxx')
call test_value((uxy+1.4246398245262621d-004)/uxy ,'o92_pt_partials(): uxy')
call test_value((uxz-6.2590905699873601d-005)/uxz ,'o92_pt_partials(): uxz')
call test_value((uyx+2.0511017314976633d-004)/uyx ,'o92_pt_partials(): uyx')
call test_value((uyy+3.0067825210386353d-004)/uyy ,'o92_pt_partials(): uyy')
call test_value((uyz-1.6933001353308135d-004)/uyz ,'o92_pt_partials(): uyz')
call test_value((uzx+6.2590905699873859d-005)/uzx ,'o92_pt_partials(): uzx')
call test_value((uzy+1.6933001353308129d-004)/uzy ,'o92_pt_partials(): uzy')
call test_value((uzz-1.7643306522090115d-004)/uzz ,'o92_pt_partials(): uzz')
write(0,*) 'subroutine o92_pt_partials() passed strike-slip unit test'
call o92_pt_strain(strain,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value((strain(1,1)+2.2862094355883986d-004)/strain(1,1) ,'o92_pt_strain(): strain(1,1)')
call test_value((strain(1,2)+1.7378707780119627d-004)/strain(1,2) ,'o92_pt_strain(): strain(1,2)')
call test_value((strain(1,3)+1.2874900798265365d-019)/strain(1,3) ,'o92_pt_strain(): strain(1,3)')
call test_value((strain(2,1)+1.7378707780119627d-004)/strain(2,1) ,'o92_pt_strain(): strain(2,1)')
call test_value((strain(2,2)+3.0067825210386353d-004)/strain(2,2) ,'o92_pt_strain(): strain(2,2)')
call test_value((strain(2,3)-2.7105054312137611d-020)/strain(2,3) ,'o92_pt_strain(): strain(2,3)')
call test_value((strain(3,1)+1.2874900798265365d-019)/strain(3,1) ,'o92_pt_strain(): strain(3,1)')
call test_value((strain(3,2)-2.7105054312137611d-020)/strain(3,2) ,'o92_pt_strain(): strain(3,2)')
call test_value((strain(3,3)-1.7643306522090115d-004)/strain(3,3) ,'o92_pt_strain(): strain(3,3)')
write(0,*) 'subroutine o92_pt_strain() passed strike-slip unit test'
moment = 0.0d0
moment(2) = 40.0d0
call o92_pt_disp(disp,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value((disp(1)+1.1723215364075840d-3)/disp(1) ,'o92_pt_disp(): disp(1)')
call test_value((disp(2)+2.0819985410726201d-3)/disp(2) ,'o92_pt_disp(): disp(2)')
call test_value((disp(3)+2.5315947108865925d-3)/disp(3) ,'o92_pt_disp(): disp(3)')
write(0,*) 'subroutine o92_pt_disp() passed dip-slip unit test'
call o92_pt_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz,&
                     sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value((uxx+1.5259408689661085d-004)/uxx ,'o92_pt_partials(): uxx')
call test_value((uxy+3.5441800891339358d-004)/uxy ,'o92_pt_partials(): uxy')
call test_value((uxz+8.7070529207727333d-004)/uxz ,'o92_pt_partials(): uxz')
call test_value((uyx-6.9826591999645748d-004)/uyx ,'o92_pt_partials(): uyx')
call test_value((uyy+1.1537526800074350d-003)/uyy ,'o92_pt_partials(): uyy')
call test_value((uyz-6.3453607666584869d-004)/uyz ,'o92_pt_partials(): uyz')
call test_value((uzx-8.7070529207727322d-004)/uzx ,'o92_pt_partials(): uzx')
call test_value((uzy+6.3453607666584847d-004)/uzy ,'o92_pt_partials(): uzy')
call test_value((uzz-4.3544892230134886d-004)/uzz ,'o92_pt_partials(): uzz')
write(0,*) 'subroutine o92_pt_partials() passed dip-slip unit test'
call o92_pt_strain(strain,sta_coord,evdp,dip,moment,lambda,shear_modulus)
call test_value((strain(1,1)+1.5259408689661085d-004)/strain(1,1) ,'o92_pt_strain(): strain(1,1)')
call test_value((strain(1,2)-1.7192395554153195d-004)/strain(1,2) ,'o92_pt_strain(): strain(1,2)')
call test_value((strain(1,3)+5.4210108624275222d-020)/strain(1,3) ,'o92_pt_strain(): strain(1,3)')
call test_value((strain(2,1)-1.7192395554153195d-004)/strain(2,1) ,'o92_pt_strain(): strain(2,1)')
call test_value((strain(2,2)+1.1537526800074350d-003)/strain(2,2) ,'o92_pt_strain(): strain(2,2)')
call test_value((strain(2,3)-1.0842021724855044d-019)/strain(2,3) ,'o92_pt_strain(): strain(2,3)')
call test_value((strain(3,1)+5.4210108624275222d-020)/strain(3,1) ,'o92_pt_strain(): strain(3,1)')
call test_value((strain(3,2)-1.0842021724855044d-019)/strain(3,2) ,'o92_pt_strain(): strain(3,2)')
call test_value((strain(3,3)-4.3544892230134886d-004)/strain(3,3) ,'o92_pt_strain(): strain(3,3)')
write(0,*) 'subroutine o92_pt_strain() passed dip-slip unit test'
! moment = 0.0d0
! moment(3) = 40.0d0
! call o92_pt_disp(disp,sta_coord,evdp,dip,moment,lambda,shear_modulus)
! call test_value((disp(1)+3.5716333691232548d-4)/disp(1) ,'o92_pt_disp(): disp(1)')
! call test_value((disp(2)-3.5310854913572190d-4)/disp(2) ,'o92_pt_disp(): disp(2)')
! call test_value((disp(3)+2.0068126969335952d-4)/disp(3) ,'o92_pt_disp(): disp(3)')
! write(0,*) 'subroutine o92_pt_disp() passed tensile-slip unit test'
write(0,*) 'example from Okada (1985) passed unit test'
write(0,*)

! sta_coord(1) = -12.0d0
! sta_coord(2) = 6.2d0
! sta_coord(3) = 13.0d0
! evdp = 7.0d0
! moment(1) = -27.0d0
! moment(2) = 13.0d0
! moment(3) = 0.0d0
! call o92_pt_disp(disp,sta_coord,evdp,dip,moment,lambda,shear_modulus)
! call test_value((disp(1)-1.9245175992850267d-004)/disp(1) ,'o92_pt_disp(): disp(1)')
! call test_value((disp(2)+1.9098044959056914d-004)/disp(2) ,'o92_pt_disp(): disp(2)')
! call test_value((disp(3)-8.1576627770328566d-005)/disp(3) ,'o92_pt_disp(): disp(3)')

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
call o92pt(ux,uy,uz,xin,yin,stdp,evdp,dipin,rakin,area,slip,vp,vs,dens)
write(0,*) ux,uy,uz
call o92ptstn(strain,xin,yin,stdp,evdp,dipin,rakin,area,slip,vp,vs,dens)
write(0,*) strain(1,:)
write(0,*) strain(2,:)
write(0,*) strain(3,:)

write(0,*) 'okada92_module unit test passed'
end



subroutine test_value(value,string)
implicit none
double precision :: value
character(len=*) :: string
if (dabs(value).gt.1.0d-10) then
    write(0,*) 'tri_disloc unit test failed for '//trim(string)
    stop
else
    write(0,*) 'tri_disloc unit test passed for '//trim(string)
endif
return
end
