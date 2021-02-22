program main
use io, only: stdout
use test, only: test_value
use tri_disloc
implicit none
double precision :: disp(3), strain(3,3), tri_coord(3,4), sta_coord(3), slip(3)
double precision :: origin(3), normal(3), pointOfInterest(3), depthOfPlaneAtPoint
double precision :: p1(3), p2(3), p3(3), p4(3), p5(3), p6(3),center(3)
double precision :: strike(3), updip(3)

y1=1.0d0
y2=2.0d0
y3=3.0d0
a=4.0d0
b=0.32d0
nu=0.25d0
B1=0.1d0
B2=0.2d0
B3=0.3d0

call angular_disloc_vars()
call test_value(sinb,0.31456656061611776d0,'angular_disloc_vars(): sinb')
call test_value(cosb,0.94923541808244083d0,'angular_disloc_vars(): cosb')
call test_value(cotb,3.0175979806093984d0,'angular_disloc_vars(): cotb')
call test_value(z1,5.5357362340875405d-003,'angular_disloc_vars(): z1')
call test_value(z3,3.1622728148634400d0,'angular_disloc_vars(): z3')
call test_value(y3bar,11.000000000000000d0,'angular_disloc_vars(): y3bar')
call test_value(z1bar,4.4094675848597360d0,'angular_disloc_vars(): z1bar')
call test_value(z3bar,10.127023038290732d0,'angular_disloc_vars(): z3bar')
call test_value(R,3.7416573867739413d0,'angular_disloc_vars(): R')
call test_value(R2,14.000000000000000d0,'angular_disloc_vars(): R2')
call test_value(R3,52.383203414835180d0,'angular_disloc_vars(): R3')
call test_value(R5,733.36484780769251d0,'angular_disloc_vars(): R5')
call test_value(Rbar,11.224972160321824d0,'angular_disloc_vars(): Rbar')
call test_value(R2bar,126.00000000000000d0,'angular_disloc_vars(): R2bar')
call test_value(R3bar,1414.3464922005498d0,'angular_disloc_vars(): R3bar')
call test_value(R4bar,15876.000000000000d0,'angular_disloc_vars(): R4bar')
call test_value(R5bar,178207.65801726928d0,'angular_disloc_vars(): R5bar')
call test_value(F,1.0152035956202865d0,'angular_disloc_vars(): F')
call test_value(Fbar,2.9255063385751567d-002,'angular_disloc_vars(): Fbar')
write(stdout,*) 'subroutine angular_disloc_vars() passed unit test'
write(stdout,*)


pointOfInterest(1) = 7.5d0
pointOfInterest(2) = 3.6d0
pointOfInterest(3) = 3.2d0
origin(1) =  1.5d0
origin(2) = -2.7d0
origin(3) =  3.2d0
normal(1) =  0.6d0
normal(2) = -0.69282032d0
normal(3) = -0.4d0
call verticalLinePlaneIntersect(depthOfPlaneAtPoint,pointOfInterest,origin,normal)
call test_value(depthOfPlaneAtPoint,1.28807996d0,'verticalLinePlaneIntersect()')
write(stdout,*) 'subroutine verticalLinePlaneIntersect() passed unit test'
write(stdout,*)


call test_value(v1B1(),0.035526830150259d0,'v1B1()')
call test_value(v2B1(),0.007506783445002d0,'v2B1()')
call test_value(v3B1(),0.003462691554128d0,'v3B1()')
call test_value(v1B2(),0.010574942146227d0,'v1B2()')
call test_value(v2B2(),0.120415543165994d0,'v2B2()')
call test_value(v3B2(),-0.001863937837612d0,'v3B2()')
call test_value(v1B3(),0.003083428958074d0,'v1B3()')
call test_value(v2B3(),-0.074712144051146d0,'v2B3()')
call test_value(v3B3(),0.099928400921349d0,'v3B3()')
write(stdout,*) 'functions viBj() for displacements passed unit tests'
write(stdout,*)


call ang_disloc_disp(disp,y1,y2,y3,a,b,nu,B1,B2,B3)
call test_value(disp(1),0.006592700131694d0,'ang_disloc_disp(): disp(1)')
call test_value(disp(2),0.002420143762355d0,'ang_disloc_disp(): disp(2)')
call test_value(disp(3),0.029952001864295d0,'ang_disloc_disp(): disp(3)')
write(stdout,*) 'subroutine ang_disloc_disp() passed unit test'
write(stdout,*)

call test_value(e11_termB1a(),0.008400368790923d0,'e11_termB1a()')
call test_value(e11_termB1b(),5.858996497052338d-06,'e11_termB1b()')
call test_value(e11_termB2a(),0.021993745473790d0,'e11_termB2a()')
call test_value(e11_termB2b(),6.367379316451323d-05,'e11_termB2b()')
call test_value(e11_termB3a(),(-0.014208301891610d0),'e11_termB3a()')
call test_value(e11_termB3b(),1.690248925356104d-04,'e11_termB3b()')
call test_value(e22_termB1a(),(-0.020631844330829d0),'e22_termB1a()')
call test_value(e22_termB1b(),6.527964013639348d-06,'e22_termB1b()')
call test_value(e22_termB2a(),(-0.041802214117666d0),'e22_termB2a()')
call test_value(e22_termB2b(),3.614748844589103d-04,'e22_termB2b()')
call test_value(e22_termB3a(),0.009510531315538d0,'e22_termB3a()')
call test_value(e22_termB3b(),5.062474014402907d-04,'e22_termB3b()')
call test_value(e33_termB1a(),0.004012383803250d0,'e33_termB1a()')
call test_value(e33_termB1b(),(-5.173267965506375d-05),'e33_termB1b()')
call test_value(e33_termB2a(),(0.001589791985872d0),'e33_termB2a()')
call test_value(e33_termB2b(),(-0.001669794447960d0),'e33_termB2b()')
call test_value(e33_termB3a(),0.020233119241681d0,'e33_termB3a()')
call test_value(e33_termB3b(),(-0.001419832184834d0),'e33_termB3b()')
call test_value(e12_termB1a(),(-0.008579559459718d0),'e12_termB1a()')
call test_value(e12_termB1b(),0.023064656529654d0,'e12_termB1b()')
call test_value(e12_termB2a(),(-0.006203119918170d0),'e12_termB2a()')
call test_value(e12_termB2b(),(-0.018531757523155d0),'e12_termB2b()')
call test_value(e12_termB3a(),6.131170880527790d-04,'e12_termB3a()')
call test_value(e12_termB3b(),(-5.409727323827170d-04),'e12_termB3b()')
call test_value(e13_termB1a(),0.004475322868068d0,'e13_termB1a()')
call test_value(e13_termB1b(),(-0.007189035380886d0),'e13_termB1b()')
call test_value(e13_termB2a(),1.032714122351998d-04,'e13_termB2a()')
call test_value(e13_termB2b(),2.438873637085293d-04,'e13_termB2b()')
call test_value(e13_termB3a(),0.001970489774467d0,'e13_termB3a()')
call test_value(e13_termB3b(),(-0.004134810449308d0),'e13_termB3b()')
call test_value(e23_termB1a(),(-5.876101709989005d-04),'e23_termB1a()')
call test_value(e23_termB1b(),7.070956009099477d-04,'e23_termB1b()')
call test_value(e23_termB2a(),0.020210757289498d0,'e23_termB2a()')
call test_value(e23_termB2b(),(-0.003151829718226d0),'e23_termB2b()')
call test_value(e23_termB3a(),(-0.002693279090971d0),'e23_termB3a()')
call test_value(e23_termB3b(),(-0.010174116926604d0),'e23_termB3b()')
write(stdout,*) 'functions eij_termBk() for strains passed unit tests'
write(stdout,*)

call ang_disloc_strain(strain,y1,y2,y3,a,b,nu,B1,B2,B3)
call test_value(strain(1,1),1.0403235324105481d-3,'ang_disloc_strain(): strain(1,1)')
call test_value(strain(1,2),-3.4768224745703653d-3,'ang_disloc_strain(): strain(1,2)')
call test_value(strain(1,3),-8.5123569854524128d-4,'ang_disloc_strain(): strain(1,3)')
call test_value(strain(2,1),-3.4768224745703653d-3,'ang_disloc_strain(): strain(2,1)')
call test_value(strain(2,2),-7.3456458682292997d-3,'ang_disloc_strain(): strain(2,2)')
call test_value(strain(2,3),-4.3648474802720008d-4,'ang_disloc_strain(): strain(2,3)')
call test_value(strain(3,1),-8.5123569854524128d-4,'ang_disloc_strain(): strain(3,1)')
call test_value(strain(3,2),-4.3648474802720008d-4,'ang_disloc_strain(): strain(3,2)')
call test_value(strain(3,3),6.0240507369960488d-3,'ang_disloc_strain(): strain(3,3)')
write(stdout,*) 'subroutine ang_disloc_strain() passed unit test'
write(stdout,*)

! Example in Meade (2007), Figure 7
tri_coord(1,1) = 20.0d0
tri_coord(2,1) = 50.0d0
tri_coord(3,1) = 5.0d0
tri_coord(1,2) = 80.0d0
tri_coord(2,2) = 80.0d0
tri_coord(3,2) = 0.0d0
tri_coord(1,3) = 50.0d0
tri_coord(2,3) = 30.0d0
tri_coord(3,3) = 20.0d0
tri_coord(:,4) = 0.0d0
sta_coord(1) = 1.0d0
sta_coord(2) = 2.0d0
sta_coord(3) = 3.0d0
slip(1) = -0.8d0
slip(2) = 0.2d0
slip(3) = 0.6d0
call tri_disloc_disp(disp,sta_coord,tri_coord,nu, slip)
call test_value(disp(1),0.018484710048163d0,'tri_disloc_disp(): disp(1)')
call test_value(disp(2),0.024890597153371d0,'tri_disloc_disp(): disp(2)')
call test_value(disp(3),-1.195473671420089d-4,'tri_disloc_disp(): disp(3)')
write(stdout,*) 'subroutine tri_disloc_strain() passed unit test'
write(stdout,*)

call tri_disloc_strain(strain,sta_coord,tri_coord,nu, slip)
call test_value(strain(1,1),4.466731272802531d-4,'tri_disloc_strain(): strain(1,1)')
call test_value(strain(2,2),3.342701850596506d-4,'tri_disloc_strain(): strain(2,2)')
call test_value(strain(3,3),-2.531425586976375d-4,'tri_disloc_strain(): strain(3,3)')
call test_value(strain(1,2),7.153107700399887d-4,'tri_disloc_strain(): strain(1,2)')
call test_value(strain(1,3),-9.831385046414070d-5,'tri_disloc_strain(): strain(1,3)')
call test_value(strain(2,3),-9.961453222958287d-5,'tri_disloc_strain(): strain(2,3)')
write(stdout,*) 'subroutine tri_disloc_strain() passed unit test'
write(stdout,*)

sta_coord(1) = 51.0d0
sta_coord(2) = 52.0d0
sta_coord(3) = 45.0d0
call tri_disloc_disp(disp,sta_coord,tri_coord,nu, slip)
call test_value(disp(1),-0.005500148465194d0,'tri_disloc_disp(): disp(1)')
call test_value(disp(2),-0.009000676865995d0,'tri_disloc_disp(): disp(2)')
call test_value(disp(3),-0.022649202364603d0,'tri_disloc_disp(): disp(3)')
write(stdout,*) 'subroutine tri_disloc_strain() passed unit test'
write(stdout,*)

call tri_disloc_strain(strain,sta_coord,tri_coord,nu, slip)
call test_value(strain(1,1),-3.930417208817573d-4,'tri_disloc_strain(): strain(1,1)')
call test_value(strain(2,2),-6.756624218076405d-4,'tri_disloc_strain(): strain(2,2)')
call test_value(strain(3,3),0.001610270039414d0,'tri_disloc_strain(): strain(3,3)')
call test_value(strain(1,2),-1.833690817949515d-4,'tri_disloc_strain(): strain(1,2)')
call test_value(strain(1,3),-2.461216982562142d-4,'tri_disloc_strain(): strain(1,3)')
call test_value(strain(2,3),1.335040694307682d-4,'tri_disloc_strain(): strain(2,3)')
write(stdout,*) 'subroutine tri_disloc_strain() passed unit test'
write(stdout,*)


p1(1) =  -4.4600515603497746d0
p1(2) =  -8.6569884424234189d-2
p1(3) =   5.0199612757737722d0
p2(1) =  -9.5108247453189811d0
p2(2) =  -8.4314941422125749d0
p2(3) =  -8.1220467753405305d0
p3(1) =  -7.2401514824260911d0
p3(2) =  -5.2259626871753104d0
p3(3) =   7.2451208789322159d0
call tri_center(center,p1,p2,p3)
call test_value(center(1),-7.0703425960316153d0,'tri_center(): center(1)')
call test_value(center(2),-4.5813422379373732d0,'tri_center(): center(2)')
call test_value(center(3), 1.3810117931218191d0,'tri_center(): center(3)')
write(stdout,*) 'subroutine tri_center() passed unit test'
write(stdout,*)


call tri_geometry(normal,strike,updip,p1,p2,p3)
call test_value(normal(1), 0.87409206677578921d0,'tri_geometry(): normal(1)')
call test_value(normal(2),-0.48495276544177424d0,'tri_geometry(): normal(2)')
call test_value(normal(3), 2.7997751516951008d-2,'tri_geometry(): normal(3)')
call test_value(strike(1), 0.48514294821244353d0,'tri_geometry(): strike(1)')
call test_value(strike(2), 0.87443485737917515d0,'tri_geometry(): strike(2)')
call test_value(strike(3),                   0d0,'tri_geometry(): strike(3)')
call test_value(updip(1),-2.4482209854662641d-2,'tri_geometry(): updip(1)')
call test_value(updip(2), 1.3582911714253024d-2,'tri_geometry(): updip(2)')
call test_value(updip(3), 0.99960798611755564d0,'tri_geometry(): updip(3)')


call tri_geo2cart(p4,p5,p6,p1,p2,p3,'km')
call test_value(p4(1), 290.54853475303207d0,'tri_geo2cart(): point1(1)')
call test_value(p4(2), 499.43977645236936d0,'tri_geo2cart(): point1(2)')
call test_value(p4(3), 5.0199612757737722d0,'tri_geo2cart(): point1(3)')
call test_value(p5(1),-268.63732246707610d0,'tri_geo2cart(): point2(1)')
call test_value(p5(2),-428.70221439679267d0,'tri_geo2cart(): point2(2)')
call test_value(p5(3),-8.1220467753405305d0,'tri_geo2cart(): point2(3)')
call test_value(p6(1),-18.803795453256310d0,'tri_geo2cart(): point3(1)')
call test_value(p6(2),-71.680853403417458d0,'tri_geo2cart(): point3(2)')
call test_value(p6(3), 7.2451208789322159d0,'tri_geo2cart(): point3(3)')


write(stdout,*) 'tri_disloc unit test passed'
end




! subroutine test_value(value,string)
! implicit none
! double precision :: value
! character(len=*) :: string
! if (dabs(value).gt.1.0d-10) then
!     write(0,*) 'tri_disloc unit test failed for '//trim(string)
!     stop
! else
!     write(0,*) 'tri_disloc unit test passed for '//trim(string)
! endif
! return
! end
