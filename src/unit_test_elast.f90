program main

use test, only: test_value
use io, only: stdout

use elast

implicit none

double precision :: strain(3,3), stress(3,3), normal(3), traction(3), tn, ts, td, max_shear
double precision :: poisson, lame, shearmod
integer :: ierr

! subroutine strain2stress(strain,lame_param,shear_modulus,stress)
strain(1,1) =  3.6339344890890711d-002
strain(2,2) =  5.5369613435628490d-002
strain(3,3) =  9.7093061723786969d-002
strain(1,2) =  4.3088477821184440d-002
strain(1,3) =  8.8046778867983427d-002
strain(2,3) =  2.2125123077825250d-003
strain(2,1) = strain(1,2)
strain(3,1) = strain(1,3)
strain(3,2) = strain(2,3)
call strain2stress(strain,42.0d9,37.0d9,stress)
call test_value(stress(1,1),10618796364.038771d0,'strain2stress: stress(1,1)')
call test_value(stress(2,1),3188547358.7676487d0,'strain2stress: stress(2,1)')
call test_value(stress(3,1),6515461636.2307739d0,'strain2stress: stress(3,1)')
call test_value(stress(1,2),3188547358.7676487d0,'strain2stress: stress(1,2)')
call test_value(stress(2,2),12027036236.349367d0,'strain2stress: stress(2,2)')
call test_value(stress(3,2),163725910.77590686d0,'strain2stress: stress(3,2)')
call test_value(stress(1,3),6515461636.2307739d0,'strain2stress: stress(1,3)')
call test_value(stress(2,3),163725910.77590686d0,'strain2stress: stress(2,3)')
call test_value(stress(3,3),15114571409.673096d0,'strain2stress: stress(3,3)')

! subroutine stress2traction(stress_tensor,normal_vector,traction_vector)
normal(1) =  0.517088d0
normal(2) = -0.775632d0
normal(3) =  0.361961d0
call stress2traction(stress,normal,traction)
call test_value(traction(1), 5376051707.6165333d0,'stress2traction: traction(1)')
call test_value(traction(2),-7620526371.6986370d0,'stress2traction: traction(2)')
call test_value(traction(3), 8712954690.2444344d0,'stress2traction: traction(3)')

! subroutine traction_components(trac_vec,nor_vec,trac_nor,trac_str,trac_dip)
call traction_components(traction,normal,tn,ts,td)
call test_value(tn,11844356671.507679d0,'traction_components: normal traction')
call test_value(ts,246037931.97703600d0,'traction_components: strike traction')
call test_value(td,4747682533.5481758d0,'traction_components: dip traction')

! subroutine max_shear_stress(stress,max_shear)
call max_shear_stress(stress,max_shear)
call test_value(max_shear,7611365305.0594463d0,'max_shear_stress')

! subroutine read_halfspace_file(halfspace_file,poisson,shearmod,lame,ierr)
open(unit=11,file='unit_test_elast_read_halfspace_file',status='unknown')
write(11,*) 'lame 40e9 shear 30e9'
close(11)
call read_halfspace_file('unit_test_elast_read_halfspace_file',poisson,shearmod,lame,ierr)
call test_value(poisson,0.28571428571428570d0,'read_halfspace_file: lame/shear->poisson')
call test_value(shearmod,30d9,'read_halfspace_file: lame/shear->shear modulus')
call test_value(lame,40d9,'read_halfspace_file: lame/shear->lame')

open(unit=11,file='unit_test_elast_read_halfspace_file',status='unknown')
write(11,*) 'shear 40e9 lame 30e9'
close(11)
call read_halfspace_file('unit_test_elast_read_halfspace_file',poisson,shearmod,lame,ierr)
call test_value(poisson,0.21428571428571427d0,'read_halfspace_file: shear/lame->poisson')
call test_value(shearmod,40d9,'read_halfspace_file: shear/lame->shear modulus')
call test_value(lame,30d9,'read_halfspace_file: shear/lame->lame')

open(unit=11,file='unit_test_elast_read_halfspace_file',status='unknown')
write(11,*) 'poisson 0.25 lame 20e9'
close(11)
call read_halfspace_file('unit_test_elast_read_halfspace_file',poisson,shearmod,lame,ierr)
call test_value(poisson,0.25d0,'read_halfspace_file: poisson/lame->poisson')
call test_value(shearmod,20d9,'read_halfspace_file: poisson/lame->shear modulus')
call test_value(lame,20d9,'read_halfspace_file: poisson/lame->lame')

open(unit=11,file='unit_test_elast_read_halfspace_file',status='unknown')
write(11,*) 'shearmod 60e9 poisson 0.25'
close(11)
call read_halfspace_file('unit_test_elast_read_halfspace_file',poisson,shearmod,lame,ierr)
call test_value(poisson,0.25d0,'read_halfspace_file: shearmod/poisson->poisson')
call test_value(shearmod,60d9,'read_halfspace_file: shearmod/poisson->shear modulus')
call test_value(lame,60d9,'read_halfspace_file: shearmod/poisson->lame')

open(unit=11,file='unit_test_elast_read_halfspace_file',status='unknown')
write(11,*) 'young 100e9 poisson 0.25'
close(11)
call read_halfspace_file('unit_test_elast_read_halfspace_file',poisson,shearmod,lame,ierr)
call test_value(poisson,0.25d0,'read_halfspace_file: young/poisson->poisson')
call test_value(shearmod,40d9,'read_halfspace_file: young/poisson->shear modulus')
call test_value(lame,40d9,'read_halfspace_file: young/poisson->lame')

open(unit=11,file='unit_test_elast_read_halfspace_file',status='unknown')
close(11,status='delete')

write(stdout,*) 'elast_module unit test passed'
end
