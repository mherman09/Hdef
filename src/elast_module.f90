module elast

public :: strain2stress
public :: stress2traction
public :: traction_components
public :: max_shear_stress

public :: read_halfspace_file

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

subroutine strain2stress(strain,lame_param,shear_modulus,stress)
!----
! Calculate the stress tensor from the strain tensor and elastic moduli, assuming a linear,
! isotropic medium.
!----

implicit none

! Arguments
double precision :: stress(3,3), strain(3,3), lame_param, shear_modulus

! Local variables
double precision :: diag

diag = strain(1,1) + strain(2,2) + strain(3,3)

! Diagonal terms
stress(1,1) = lame_param*diag + 2.0d0*shear_modulus*strain(1,1)
stress(2,2) = lame_param*diag + 2.0d0*shear_modulus*strain(2,2)
stress(3,3) = lame_param*diag + 2.0d0*shear_modulus*strain(3,3)

! Off-diagonal terms
stress(1,2) = 2.0d0*shear_modulus*strain(1,2)
stress(1,3) = 2.0d0*shear_modulus*strain(1,3)
stress(2,3) = 2.0d0*shear_modulus*strain(2,3)
stress(2,1) = stress(1,2)
stress(3,1) = stress(1,3)
stress(3,2) = stress(2,3)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine stress2traction(stress_tensor,normal_vector,traction_vector)
!----
! Calculate the traction vector resolved onto a plane from the stress tensor and the normal vector
! to the plane.
!----

use io, only: stderr
implicit none

! Arguments
double precision :: stress_tensor(3,3), normal_vector(3), traction_vector(3)

! Local variables
integer :: i, j
double precision :: magnitude

! Initialize traction vector
traction_vector = 0.0d0

! Normalize the normal vector
magnitude = normal_vector(1)*normal_vector(1) + &
            normal_vector(2)*normal_vector(2) + &
            normal_vector(3)*normal_vector(3)

if (magnitude.lt.1.0d-8) then
    write(stderr,*) 'stress2traction: normal vector magnitude is zero, setting traction to zero'
    return
endif

normal_vector = normal_vector/magnitude

! traction_vector = stress_tensor*normal_vector
do i = 1,3
    do j = 1,3
        traction_vector(i) = traction_vector(i) + stress_tensor(i,j)*normal_vector(j)
    enddo
enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine traction_components(trac_vec,nor_vec,trac_nor,trac_str,trac_dip)
!----
! Calculate the normal, strike-parallel, and up-dip components of the traction vector on a plane
! defined by its normal vector.
!----

use geom, only: normal2strike, normal2updip
use algebra, only: dot_product

implicit none

! Arguments
double precision :: trac_vec(3), nor_vec(3), trac_nor, trac_str, trac_dip

! Local variables
double precision :: str_vec(3), dip_vec(3)

! Calculate unit vectors parallel to strike and up-dip
call normal2strike(nor_vec,str_vec)
call normal2updip(nor_vec,dip_vec)

! Project the traction vector onto normal, along-strike, and horizontal up-dip directions
call dot_product(trac_vec,nor_vec,trac_nor)
call dot_product(trac_vec,str_vec,trac_str)
call dot_product(trac_vec,dip_vec,trac_dip)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine max_shear_stress(stress,max_shear)
!----
! Calculate maximum shear stress from a stress tensor
!----

implicit none

! Arguments
double precision :: stress(3,3), max_shear

! Local variables
double precision :: s11_s22, s11_s33, s22_s33, s12s12, s13s13, s23s23

s11_s22 = stress(1,1)-stress(2,2)
s11_s33 = stress(1,1)-stress(3,3)
s22_s33 = stress(2,2)-stress(3,3)

s12s12 = stress(1,2)*stress(1,2)
s13s13 = stress(1,3)*stress(1,3)
s23s23 = stress(2,3)*stress(2,3)

max_shear = (1.0d0/6.0d0)*(s11_s22*s11_s22+s11_s33*s11_s33+s22_s33*s22_s33) + s12s12+s13s13+s23s23

! Make sure it has units of Pa
max_shear = dsqrt(max_shear)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine read_halfspace_file(halfspace_file,poisson,shearmod,lame,ierr)
!----
! Read half-space elastic moduli from file (if provided), in format:
!     lbl1 value1 lbl2 value2 [lbl3 value3]
!
! For example, if you have the shear modulus and Young's modulus, the file would be:
!     shearmod 40e9 young 8e10
!----

use io, only: stderr, line_count, fileExists

implicit none

! Arguments
character(len=*) :: halfspace_file
double precision :: poisson, shearmod, lame

! Local variables
integer :: i, j, ios, ierr
character(len=512) :: input_line
character(len=32) :: lbl(4)
character(len=8) :: halfspace_file_mode
double precision :: val(3)


ierr = 0

! If the file is not defined, use defaults set in gcmdln
if (halfspace_file.eq.'') then
    return
endif


! Check whether the half-space file exists
if (.not.fileExists(halfspace_file)) then
    write(stderr,*) 'read_halfspace: no half-space file named '//trim(halfspace_file)//' found'
    write(stderr,*) 'using default elastic moduli values'
    return
endif


! Initialize modulus variables
lbl = ''
val = 0.0d0


! Open the file and read the half-space elastic moduli
halfspace_file_mode = 'modern'
open(unit=15,file=halfspace_file,status='old')
read(15,'(A)') input_line

! For back-compatibility, can read parameters in these formats (print warning):
!     vp vs dens
!     "lame" lame shear_mod
read(input_line,*,end=1501,err=1500,iostat=ios) lbl(1:4)
1501 if (ios.ne.0) then
    halfspace_file_mode = 'legacy'
    ! Found fewer than 4 inputs
    if (lbl(1).eq.'lame'.or.lbl(1).eq.'Lame') then
        read(lbl(2),*,err=1500,iostat=ios) val(1)
        read(lbl(3),*,err=1500,iostat=ios) val(2)
        lbl(1) = 'lame'
        lbl(2) = 'shearmod'
        lbl(3:4) = ''
        write(stderr,*) 'read_halfspace: WARNING: you are using a legacy half-space format ("lame" lame shear)'
        write(stderr,*) 'In the future, you should use labeled inputs, e.g.: lame 40e9 shearmod 40e9'
    else
        read(lbl(1),*,err=1500,iostat=ios) val(1)
        read(lbl(2),*,err=1500,iostat=ios) val(2)
        read(lbl(3),*,err=1500,iostat=ios) val(3)
        lbl(1) = 'vp'
        lbl(2) = 'vs'
        lbl(3) = 'dens'
        write(stderr,*) 'read_halfspace: WARNING: you are using a legacy half-space format (vp vs dens)'
        write(stderr,*) 'In the future, you should use labeled inputs, e.g.: vp 6800 vs 3926 dens 3000'
    endif
endif

! Read in updated format, which can take three (labeled) inputs for a minimum of 4 items on line
! Examples: shear 40e9 lame 35e9
!           vp 6500 vs 3700 dens 3000
if (halfspace_file_mode.eq.'modern') then
    read(input_line,*,iostat=ios) lbl(1), val(1), lbl(2), val(2), lbl(3), val(3)

    ! Then check for two inputs
    if (ios.ne.0) then
        read(input_line,*,iostat=ios) lbl(1), val(1), lbl(2), val(2)
    endif
endif

1500 if (ios.ne.0) then
    write(stderr,*) 'read_halfspace: error reading modulus labels and values'
    write(stderr,*) 'offending line: '//trim(input_line)
    ierr = 1
    return
endif

close(15)


! Check labels make sense
do i = 1,2
    if (lbl(i).eq.'mu' .or. lbl(i).eq.'shear_modulus' .or. lbl(i).eq.'shear' .or. &
                                                                     lbl(i).eq.'shearmod') then
        lbl(i) = 'shearmod'
    elseif (lbl(i).eq.'lame' .or. lbl(i).eq.'lambda') then
        lbl(i) = 'lame'
    elseif (lbl(i).eq.'poisson' .or. lbl(i).eq.'nu') then
        lbl(i) = 'poisson'
    elseif (lbl(i).eq.'young') then
        lbl(i) = 'young'
    elseif (lbl(i).eq.'vp'.or.lbl(i).eq.'vs'.or.lbl(i).eq.'dens') then
        ! Okay
    else
        write(stderr,*) 'read_halfspace: elastic modulus "',trim(lbl(i)),'" not implemented'
        write(stderr,*) 'specify two of: shear_modulus, lame, poisson, or young OR vp vs dens'
        ierr = 2
        return
    endif
enddo


! Calculate Poisson's ratio, Lame parameter, and shear modulus
do i = 1,2

    j = mod(i,2)+1

    if (lbl(i).eq.'shearmod'.and.lbl(j).eq.'lame') then
        ! Calculate Poisson's ratio from Lame's parameter and shear modulus
        shearmod = val(i)
        lame = val(j)
        poisson = lame/(2.0d0*(lame+shearmod))

    elseif (lbl(i).eq.'shearmod'.and.lbl(j).eq.'poisson') then
        ! Calculate Lame's parameter from Poisson's ratio and shear modulus
        shearmod = val(i)
        poisson = val(j)
        lame = 2.0d0*shearmod*poisson/(1.0d0-2.0d0*poisson)

    elseif (lbl(i).eq.'shearmod'.and.lbl(j).eq.'young') then
        ! Calculate Lame's parameter and Poisson's ratio from Young's modulus and shear modulus
        shearmod = val(i)
        lame = shearmod*(val(j)-2.0d0*shearmod)/(3.0d0*shearmod-val(j))
        poisson = lame/(2.0d0*(lame+shearmod))

    elseif (lbl(i).eq.'lame'.and.lbl(j).eq.'poisson') then
        ! Calculate shear modulus from Lame's parameter and Poisson's ratio
        lame = val(i)
        poisson = val(j)
        shearmod = lame*(1.0d0-2.0d0*poisson)/(2.0d0*poisson)

    elseif (lbl(i).eq.'lame'.and.lbl(j).eq.'young') then
        ! Calculate shear modulus and Poisson's ratio from Lame's parameter and Young's modulus
        lame = val(i)
        poisson = 2.0d0*lame/(val(j)+lame+sqrt(val(j)*val(j)+9.0d0*lame*lame+2.0d0*lame*val(j)))
        shearmod = lame*(1.0d0-2.0d0*poisson)/(2.0d0*poisson)

    elseif (lbl(i).eq.'poisson'.and.lbl(j).eq.'young') then
        ! Calculate shear modulus and Lame's parameter from Poisson's ratio and Young's modulus
        poisson = val(i)
        shearmod = val(j)/(2.0d0*(1.0d0+poisson))
        lame = 2.0d0*shearmod*poisson/(1.0d0-2.0d0*poisson)

    endif
enddo

if (lbl(1).eq.'vp'.and.lbl(2).eq.'vs'.and.lbl(3).eq.'dens') then
    ! Calculate shear modulus, Lame's parameter, and Poisson's ratio from seismic velocities and density
    shearmod = val(2)*val(2)*val(3)
    lame = val(1)*val(1)*val(3) - 2.0d0*shearmod
    poisson = lame/(2.0d0*(lame+shearmod))
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!















! !--------------------------------------------------------------------------------------------------!
!
!     subroutine calc_traction_components(traction,normal,strike,updip,traction_components)
!     implicit none
!     ! I/O variables
!     double precision :: traction(3), normal(3), strike(3), updip(3), traction_components(3)
!     ! Local variables
!     integer :: i
!
!     ! Project the traction vector on normal, along-strike, and horizontal up-dip directions
!     traction_components = 0.0d0
!     do i = 1,3
!         traction_components(1) = traction_components(1) + traction(i)*normal(i)
!         traction_components(2) = traction_components(2) + traction(i)*strike(i)
!         traction_components(3) = traction_components(3) + traction(i)*updip(i)
!     enddo
!     return
!     end subroutine calc_traction_components
!
!
! !--------------------------------------------------------------------------------------------------!
!
!     subroutine calc_plane_unit_vectors(strike,dip,unit_normal,unit_strike,unit_updip)
!     implicit none
!     ! I/O variables
!     double precision :: strike, dip, unit_normal(3), unit_strike(3), unit_updip(3)
!     ! Local variables
!     double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r=pi/180.0d0, pi2=pi/2.0d0
!     double precision :: str_rad, dip_rad
!
!     str_rad = strike*d2r
!     dip_rad = dip*d2r
!
!     unit_normal(1) = dsin(dip_rad)*dsin(str_rad+pi2)
!     unit_normal(2) = dsin(dip_rad)*dcos(str_rad+pi2)
!     unit_normal(3) = dcos(dip_rad)
!
!     unit_strike(1) = dsin(str_rad)
!     unit_strike(2) = dcos(str_rad)
!     unit_strike(3) = 0.0d0
!
!     unit_updip(1) = dcos(dip_rad)*dsin(str_rad-pi2)
!     unit_updip(2) = dcos(dip_rad)*dcos(str_rad-pi2)
!     unit_updip(3) = dsin(dip_rad)
!
!     return
!     end
!
! !--------------------------------------------------------------------------------------------------!
!
!     subroutine rotate_strain(strain,strike)
!     !----
!     ! Rotate strain matrix from (x=str, y=updip horizontal, z=up) to (x=E, y=N, z=up)
!     !----
!     implicit none
!     ! I/O variables
!     double precision, dimension(3,3) :: strain
!     double precision :: strike
!     ! Local variables
!     double precision, dimension(3,3) :: rot, rot_transpose, tmp
!     double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r=pi/180.0d0
!     rot(1,1) = dcos(d2r*strike-pi/2.0d0)
!     rot(2,2) = rot(1,1)
!     rot(1,2) = dsin(d2r*strike-pi/2.0d0)
!     rot(2,1) = -rot(1,2)
!     rot(1,3) = 0.0d0
!     rot(2,3) = 0.0d0
!     rot(3,1) = 0.0d0
!     rot(3,2) = 0.0d0
!     rot(3,3) = 1.0d0
!     call mattr(rot_transpose,rot)
!     call matmult(tmp,rot,strain)
!     call matmult(strain,tmp,rot_transpose)
!     return
!     end
!
! !--------------------------------------------------------------------------------------------------!
!
!     subroutine mattr(matout,matin)
!     !----
!     ! Transpose a 3x3 matrix
!     !----
!     implicit none
!     ! I/O variables
!     double precision, dimension(3,3) :: matout, matin
!     ! Local variables
!     integer :: i, j
!     do i = 1,3
!         do j = 1,3
!             matout(i,j) = matin(j,i)
!         enddo
!     enddo
!     return
!     end
!
! !--------------------------------------------------------------------------------------------------!
!
!     subroutine matmult(matout,mat1,mat2)
!     !----
!     ! Multiply two 3x3 matrices
!     !----
!     implicit none
!     ! I/O variables
!     double precision, dimension(3,3) :: mat1, mat2, matout
!     ! Local variables
!     integer :: i, j, k
!     do i = 1,3
!         do j = 1,3
!             matout(i,j) = 0.0d0
!             do k = 1,3
!                 matout(i,j) = matout(i,j) + mat1(i,k)*mat2(k,j)
!             enddo
!         enddo
!     enddo
!     return
!     end
! !--------------------------------------------------------------------------------------------------!
!
!
!
! subroutine calc_tractions(trac_vector,shear,shearmax,normal,coulomb,stress,trg)
! !----
! ! Calculate various tractions resolved onto a fault geometry from a stress tensor
! !----
!
! use trig, only: pi, d2r
! implicit none
!
! ! Arguments
! double precision :: trac_vector(3), shear, shearmax, normal, coulomb, stress(3,3), trg(4)
!
! ! Local variables
! integer :: i
! double precision :: str, dip, rak, fric
! double precision :: n(3), s(3), r(3)
!
! str = trg(1)*d2r
! dip = trg(2)*d2r
! rak = trg(3)*d2r
! fric = trg(4)
!
! ! Unit normal vector to target plane
! n(1) = dsin(dip)*dsin(str+pi/2.0d0)
! n(2) = dsin(dip)*dcos(str+pi/2.0d0)
! n(3) = dcos(dip)
!
! ! Traction vector is stress matrix times normal vector: t = S*n
! do i = 1,3
!     trac_vector(i) = stress(i,1)*n(1) + stress(i,2)*n(2) + stress(i,3)*n(3)
! enddo
!
! ! Normal component of traction is parallel to unit normal vector; take dot product
! normal = 0.0d0
! do i = 1,3
!     normal = normal + trac_vector(i)*n(i)
! enddo
!
! ! Shear component of traction is difference between total traction vector and normal traction
! s(1) = trac_vector(1) - normal*n(1)
! s(2) = trac_vector(2) - normal*n(2)
! s(3) = trac_vector(3) - normal*n(3)
! shearmax = dsqrt(s(1)*s(1) + s(2)*s(2) + s(3)*s(3))
!
! ! Compute unit slip vector (parallel to rake)
! r(1) =  dcos(str-pi/2.0d0)*dcos(rak) + dsin(str-pi/2.0d0)*dsin(rak)*dcos(dip)
! r(2) = -dsin(str-pi/2.0d0)*dcos(rak) + dcos(str-pi/2.0d0)*dsin(rak)*dcos(dip)
! r(3) =                                                    dsin(rak)*dsin(dip)
!
! ! Project shear component on slip vector
! shear = s(1)*r(1) + s(2)*r(2) + s(3)*r(3)
!
! ! Coulomb stress (recall sign convention: pos = dilation)
! coulomb = shear + fric*normal
!
! return
! end subroutine calc_tractions

end module elast
