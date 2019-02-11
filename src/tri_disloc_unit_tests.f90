program main
use tri_disloc_mod
implicit none
double precision :: disp(3), strain(3,3), tri_coord(3,4), sta_coord(3), slip(3)

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
call test_value((sinb-0.31456656061611776d0)/sinb,'angular_disloc_vars(): sinb')
call test_value((cosb-0.94923541808244083d0)/cosb,'angular_disloc_vars(): cosb')
call test_value((cotb-3.0175979806093984d0)/cotb,'angular_disloc_vars(): cotb')
call test_value((z1-5.5357362340875405d-003)/z1,'angular_disloc_vars(): z1')
call test_value((z3-3.1622728148634400d0)/z3,'angular_disloc_vars(): z3')
call test_value((y3bar-11.000000000000000d0)/y3bar,'angular_disloc_vars(): y3bar')
call test_value((z1bar-4.4094675848597360d0)/z1bar,'angular_disloc_vars(): z1bar')
call test_value((z3bar-10.127023038290732d0)/z3bar,'angular_disloc_vars(): z3bar')
call test_value((R-3.7416573867739413d0)/R,'angular_disloc_vars(): R')
call test_value((R2-14.000000000000000d0)/R2,'angular_disloc_vars(): R2')
call test_value((R3-52.383203414835180d0)/R3,'angular_disloc_vars(): R3')
call test_value((R5-733.36484780769251d0)/R5,'angular_disloc_vars(): R5')
call test_value((Rbar-11.224972160321824d0)/Rbar,'angular_disloc_vars(): Rbar')
call test_value((R2bar-126.00000000000000d0)/R2bar,'angular_disloc_vars(): R2bar')
call test_value((R3bar-1414.3464922005498d0)/R3bar,'angular_disloc_vars(): R3bar')
call test_value((R4bar-15876.000000000000d0)/R4bar,'angular_disloc_vars(): R4bar')
call test_value((R5bar-178207.65801726928d0)/R5bar,'angular_disloc_vars(): R5bar')
call test_value((F-1.0152035956202865d0)/F,'angular_disloc_vars(): F')
call test_value((Fbar-2.9255063385751567d-002)/Fbar,'angular_disloc_vars(): Fbar')
write(0,*) 'subroutine angular_disloc_vars() passed unit test'
write(0,*)

call test_value((v1B1()-0.035526830150259d0)/v1B1(),'v1B1()')
call test_value((v2B1()-0.007506783445002d0)/v2B1(),'v2B1()')
call test_value((v3B1()-0.003462691554128d0)/v3B1(),'v3B1()')
call test_value((v1B2()-0.010574942146227d0)/v1B2(),'v1B2()')
call test_value((v2B2()-0.120415543165994d0)/v2B2(),'v2B2()')
call test_value((v3B2()+0.001863937837612d0)/v3B2(),'v3B2()')
call test_value((v1B3()-0.003083428958074d0)/v1B3(),'v1B3()')
call test_value((v2B3()+0.074712144051146d0)/v2B3(),'v2B3()')
call test_value((v3B3()-0.099928400921349d0)/v3B3(),'v3B3()')
write(0,*) 'functions viBj() for displacements passed unit tests'
write(0,*)

call ang_disloc_disp(disp,y1,y2,y3,a,b,nu,B1,B2,B3)
call test_value((disp(1)-0.006592700131694d0)/disp(1),'ang_disloc_disp(): disp(1)')
call test_value((disp(2)-0.002420143762355d0)/disp(2),'ang_disloc_disp(): disp(2)')
call test_value((disp(3)-0.029952001864295d0)/disp(3),'ang_disloc_disp(): disp(3)')
write(0,*) 'subroutine ang_disloc_disp() passed unit test'
write(0,*)

call test_value((e11_termB1a()-0.008400368790923d0)/e11_termB1a()     ,'e11_termB1a()')
call test_value((e11_termB1b()-5.858996497052338d-06)/e11_termB1b()   ,'e11_termB1b()')
call test_value((e11_termB2a()-0.021993745473790d0)/e11_termB2a()     ,'e11_termB2a()')
call test_value((e11_termB2b()-6.367379316451323d-05)/e11_termB2b()   ,'e11_termB2b()')
call test_value((e11_termB3a()-(-0.014208301891610d0))/e11_termB3a()  ,'e11_termB3a()')
call test_value((e11_termB3b()-1.690248925356104d-04)/e11_termB3b()   ,'e11_termB3b()')
call test_value((e22_termB1a()-(-0.020631844330829d0))/e22_termB1a()  ,'e22_termB1a()')
call test_value((e22_termB1b()-6.527964013639348d-06)/e22_termB1b()   ,'e22_termB1b()')
call test_value((e22_termB2a()-(-0.041802214117666d0))/e22_termB2a()  ,'e22_termB2a()')
call test_value((e22_termB2b()-3.614748844589103d-04)/e22_termB2b()   ,'e22_termB2b()')
call test_value((e22_termB3a()-0.009510531315538d0)/e22_termB3a()     ,'e22_termB3a()')
call test_value((e22_termB3b()-5.062474014402907d-04)/e22_termB3b()   ,'e22_termB3b()')
call test_value((e33_termB1a()-0.004012383803250d0)/e33_termB1a()     ,'e33_termB1a()')
call test_value((e33_termB1b()-(-5.173267965506375d-05))/e33_termB1b(),'e33_termB1b()')
call test_value((e33_termB2a()-(0.001589791985872d0))/e33_termB2a()   ,'e33_termB2a()')
call test_value((e33_termB2b()-(-0.001669794447960d0))/e33_termB2b()  ,'e33_termB2b()')
call test_value((e33_termB3a()-0.020233119241681d0)/e33_termB3a()     ,'e33_termB3a()')
call test_value((e33_termB3b()-(-0.001419832184834d0))/e33_termB3b()  ,'e33_termB3b()')
call test_value((e12_termB1a()-(-0.008579559459718d0))/e12_termB1a()  ,'e12_termB1a()')
call test_value((e12_termB1b()-0.023064656529654d0)/e12_termB1b()     ,'e12_termB1b()')
call test_value((e12_termB2a()-(-0.006203119918170d0))/e12_termB2a()  ,'e12_termB2a()')
call test_value((e12_termB2b()-(-0.018531757523155d0))/e12_termB2b()  ,'e12_termB2b()')
call test_value((e12_termB3a()-6.131170880527790d-04)/e12_termB3a()   ,'e12_termB3a()')
call test_value((e12_termB3b()-(-5.409727323827170d-04))/e12_termB3b(),'e12_termB3b()')
call test_value((e13_termB1a()-0.004475322868068d0)/e13_termB1a()     ,'e13_termB1a()')
call test_value((e13_termB1b()-(-0.007189035380886d0))/e13_termB1b()  ,'e13_termB1b()')
call test_value((e13_termB2a()-1.032714122351998d-04)/e13_termB2a()   ,'e13_termB2a()')
call test_value((e13_termB2b()-2.438873637085293d-04)/e13_termB2b()   ,'e13_termB2b()')
call test_value((e13_termB3a()-0.001970489774467d0)/e13_termB3a()     ,'e13_termB3a()')
call test_value((e13_termB3b()-(-0.004134810449308d0))/e13_termB3b()  ,'e13_termB3b()')
call test_value((e23_termB1a()-(-5.876101709989005d-04))/e23_termB1a(),'e23_termB1a()')
call test_value((e23_termB1b()-7.070956009099477d-04)/e23_termB1b()   ,'e23_termB1b()')
call test_value((e23_termB2a()-0.020210757289498d0)/e23_termB2a()     ,'e23_termB2a()')
call test_value((e23_termB2b()-(-0.003151829718226d0))/e23_termB2b()  ,'e23_termB2b()')
call test_value((e23_termB3a()-(-0.002693279090971d0))/e23_termB3a()  ,'e23_termB3a()')
call test_value((e23_termB3b()-(-0.010174116926604d0))/e23_termB3b()  ,'e23_termB3b()')
write(0,*) 'functions eij_termBk() for strains passed unit tests'
write(0,*)

call ang_disloc_strain(strain,y1,y2,y3,a,b,nu,B1,B2,B3)
call test_value((strain(1,1)-1.0403235324105481d-3)/strain(1,1) ,'ang_disloc_strain(): strain(1,1)')
call test_value((strain(1,2)+3.4768224745703653d-3)/strain(1,2) ,'ang_disloc_strain(): strain(1,2)')
call test_value((strain(1,3)+8.5123569854524128d-4)/strain(1,3) ,'ang_disloc_strain(): strain(1,3)')
call test_value((strain(2,1)+3.4768224745703653d-3)/strain(2,1) ,'ang_disloc_strain(): strain(2,1)')
call test_value((strain(2,2)+7.3456458682292997d-3)/strain(2,2) ,'ang_disloc_strain(): strain(2,2)')
call test_value((strain(2,3)+4.3648474802720008d-4)/strain(2,3) ,'ang_disloc_strain(): strain(2,3)')
call test_value((strain(3,1)+8.5123569854524128d-4)/strain(3,1) ,'ang_disloc_strain(): strain(3,1)')
call test_value((strain(3,2)+4.3648474802720008d-4)/strain(3,2) ,'ang_disloc_strain(): strain(3,2)')
call test_value((strain(3,3)-6.0240507369960488d-3)/strain(3,3) ,'ang_disloc_strain(): strain(3,3)')
write(0,*) 'subroutine ang_disloc_strain() passed unit test'
write(0,*)

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
call test_value((disp(1)-0.018484710048163d0)/disp(1) ,'tri_disloc_disp(): disp(1)')
call test_value((disp(2)-0.024890597153371d0)/disp(2) ,'tri_disloc_disp(): disp(2)')
call test_value((disp(3)+1.195473671420089d-4)/disp(3) ,'tri_disloc_disp(): disp(3)')
write(0,*) 'subroutine tri_disloc_strain() passed unit test'
write(0,*)

call tri_disloc_strain(strain,sta_coord,tri_coord,nu, slip)
call test_value((strain(1,1)-4.466731272802531d-4)/strain(1,1) ,'tri_disloc_strain(): strain(1,1)')
call test_value((strain(2,2)-3.342701850596506d-4)/strain(2,2) ,'tri_disloc_strain(): strain(2,2)')
call test_value((strain(3,3)+2.531425586976375d-4)/strain(3,3) ,'tri_disloc_strain(): strain(3,3)')
call test_value((strain(1,2)-7.153107700399887d-4)/strain(1,2) ,'tri_disloc_strain(): strain(1,2)')
call test_value((strain(1,3)+9.831385046414070d-5)/strain(1,3) ,'tri_disloc_strain(): strain(1,3)')
call test_value((strain(2,3)+9.961453222958287d-5)/strain(2,3) ,'tri_disloc_strain(): strain(2,3)')
write(0,*) 'subroutine tri_disloc_strain() passed unit test'
write(0,*)

sta_coord(1) = 51.0d0
sta_coord(2) = 52.0d0
sta_coord(3) = 45.0d0
call tri_disloc_disp(disp,sta_coord,tri_coord,nu, slip)
call test_value((disp(1)+0.005500148465194d0)/disp(1) ,'tri_disloc_disp(): disp(1)')
call test_value((disp(2)+0.009000676865995d0)/disp(2) ,'tri_disloc_disp(): disp(2)')
call test_value((disp(3)+0.022649202364603d0)/disp(3) ,'tri_disloc_disp(): disp(3)')
write(0,*) 'subroutine tri_disloc_strain() passed unit test'
write(0,*)

call tri_disloc_strain(strain,sta_coord,tri_coord,nu, slip)
call test_value((strain(1,1)+3.930417208817573d-4)/strain(1,1) ,'tri_disloc_strain(): strain(1,1)')
call test_value((strain(2,2)+6.756624218076405d-4)/strain(2,2) ,'tri_disloc_strain(): strain(2,2)')
call test_value((strain(3,3)-0.001610270039414d0)/strain(3,3) ,'tri_disloc_strain(): strain(3,3)')
call test_value((strain(1,2)+1.833690817949515d-4)/strain(1,2) ,'tri_disloc_strain(): strain(1,2)')
call test_value((strain(1,3)+2.461216982562142d-4)/strain(1,3) ,'tri_disloc_strain(): strain(1,3)')
call test_value((strain(2,3)-1.335040694307682d-4)/strain(2,3) ,'tri_disloc_strain(): strain(2,3)')
write(0,*) 'subroutine tri_disloc_strain() passed unit test'
write(0,*)

write(0,*) 'tri_disloc unit test passed'
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
