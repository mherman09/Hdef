module elast

! TODO
! Remove redundant subroutines
!     - trac_vector & calc_tractions
!     - calc_strain_to_stress & strain2stress

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

function trac_vector(stress_tensor,normal_vec)
!----
! Compute the traction vector resolved onto a plane from the stress tensor and normal vector to
! the plane.
!----

implicit none

! I/O variables
double precision :: stress_tensor(3,3), normal_vec(3), trac_vector(3)

! Local variables
integer :: i, j
double precision :: magnitude

! Initialize traction
trac_vector = 0.0d0

! Normalize the normal vector
magnitude = normal_vec(1)*normal_vec(1) + normal_vec(2)*normal_vec(2) + normal_vec(3)*normal_vec(3)
if (magnitude.lt.1.0d-8) then
    write(0,*) 'traction: normal vector magnitude is zero, setting traction to zero'
    return
endif
normal_vec = normal_vec/magnitude

! traction_vector = stress_tensor*normal_vector
do i = 1,3
    do j = 1,3
        trac_vector(i) = trac_vector(i) + stress_tensor(i,j)*normal_vec(j)
    enddo
enddo

return
end function trac_vector

!--------------------------------------------------------------------------------------------------!

    subroutine calc_strain_to_stress(strain,vp,vs,dens,stress)
    implicit none
    ! I/O variables
    double precision :: strain(3,3), stress(3,3), vp, vs, dens
    ! Local variables
    double precision :: diag, shear_modulus, lambda

    diag = strain(1,1) + strain(2,2) + strain(3,3)
    shear_modulus = vs*vs*dens
    lambda = vp*vp*dens - 2.0d0*shear_modulus

    stress(1,1) = lambda*diag + 2.0d0*shear_modulus*strain(1,1)
    stress(2,2) = lambda*diag + 2.0d0*shear_modulus*strain(2,2)
    stress(3,3) = lambda*diag + 2.0d0*shear_modulus*strain(3,3)
    stress(1,2) = 2.0d0*shear_modulus*strain(1,2)
    stress(1,3) = 2.0d0*shear_modulus*strain(1,3)
    stress(2,3) = 2.0d0*shear_modulus*strain(2,3)
    stress(2,1) = stress(1,2)
    stress(3,1) = stress(1,3)
    stress(3,2) = stress(2,3)

    return
    end subroutine calc_strain_to_stress

!--------------------------------------------------------------------------------------------------!

    subroutine calc_traction_components(traction,normal,strike,updip,traction_components)
    implicit none
    ! I/O variables
    double precision :: traction(3), normal(3), strike(3), updip(3), traction_components(3)
    ! Local variables
    integer :: i

    ! Project the traction vector on normal, along-strike, and horizontal up-dip directions
    traction_components = 0.0d0
    do i = 1,3
        traction_components(1) = traction_components(1) + traction(i)*normal(i)
        traction_components(2) = traction_components(2) + traction(i)*strike(i)
        traction_components(3) = traction_components(3) + traction(i)*updip(i)
    enddo
    return
    end subroutine calc_traction_components


!--------------------------------------------------------------------------------------------------!

    subroutine calc_plane_unit_vectors(strike,dip,unit_normal,unit_strike,unit_updip)
    implicit none
    ! I/O variables
    double precision :: strike, dip, unit_normal(3), unit_strike(3), unit_updip(3)
    ! Local variables
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r=pi/180.0d0, pi2=pi/2.0d0
    double precision :: str_rad, dip_rad

    str_rad = strike*d2r
    dip_rad = dip*d2r

    unit_normal(1) = dsin(dip_rad)*dsin(str_rad+pi2)
    unit_normal(2) = dsin(dip_rad)*dcos(str_rad+pi2)
    unit_normal(3) = dcos(dip_rad)

    unit_strike(1) = dsin(str_rad)
    unit_strike(2) = dcos(str_rad)
    unit_strike(3) = 0.0d0

    unit_updip(1) = dcos(dip_rad)*dsin(str_rad-pi2)
    unit_updip(2) = dcos(dip_rad)*dcos(str_rad-pi2)
    unit_updip(3) = dsin(dip_rad)

    return
    end

!--------------------------------------------------------------------------------------------------!

    subroutine rotate_strain(strain,strike)
    !----
    ! Rotate strain matrix from (x=str, y=updip horizontal, z=up) to (x=E, y=N, z=up)
    !----
    implicit none
    ! I/O variables
    double precision, dimension(3,3) :: strain
    double precision :: strike
    ! Local variables
    double precision, dimension(3,3) :: rot, rot_transpose, tmp
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r=pi/180.0d0
    rot(1,1) = dcos(d2r*strike-pi/2.0d0)
    rot(2,2) = rot(1,1)
    rot(1,2) = dsin(d2r*strike-pi/2.0d0)
    rot(2,1) = -rot(1,2)
    rot(1,3) = 0.0d0
    rot(2,3) = 0.0d0
    rot(3,1) = 0.0d0
    rot(3,2) = 0.0d0
    rot(3,3) = 1.0d0
    call mattr(rot_transpose,rot)
    call matmult(tmp,rot,strain)
    call matmult(strain,tmp,rot_transpose)
    return
    end

!--------------------------------------------------------------------------------------------------!

    subroutine mattr(matout,matin)
    !----
    ! Transpose a 3x3 matrix
    !----
    implicit none
    ! I/O variables
    double precision, dimension(3,3) :: matout, matin
    ! Local variables
    integer :: i, j
    do i = 1,3
        do j = 1,3
            matout(i,j) = matin(j,i)
        enddo
    enddo
    return
    end

!--------------------------------------------------------------------------------------------------!

    subroutine matmult(matout,mat1,mat2)
    !----
    ! Multiply two 3x3 matrices
    !----
    implicit none
    ! I/O variables
    double precision, dimension(3,3) :: mat1, mat2, matout
    ! Local variables
    integer :: i, j, k
    do i = 1,3
        do j = 1,3
            matout(i,j) = 0.0d0
            do k = 1,3
                matout(i,j) = matout(i,j) + mat1(i,k)*mat2(k,j)
            enddo
        enddo
    enddo
    return
    end
!--------------------------------------------------------------------------------------------------!

subroutine strain2stress(stress,strain,lame_param,shear_modulus)
!----
! Calculate stress tensor from strain tensor, assuming linear isotropic medium
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
end subroutine strain2stress

!--------------------------------------------------------------------------------------------------!

subroutine max_shear_stress(estress,stress)
!----
! Calculate maximum shear stress
!----

implicit none

! Arguments
double precision :: estress, stress(3,3)

! Local variables
double precision :: s11_s22, s11_s33, s22_s33, s12s12, s13s13, s23s23

s11_s22 = stress(1,1)-stress(2,2)
s11_s33 = stress(1,1)-stress(3,3)
s22_s33 = stress(2,2)-stress(3,3)

s12s12 = stress(1,2)*stress(1,2)
s13s13 = stress(1,3)*stress(1,3)
s23s23 = stress(2,3)*stress(2,3)

estress = (1.0d0/6.0d0)*(s11_s22*s11_s22+s11_s33*s11_s33+s22_s33*s22_s33) + s12s12+s13s13+s23s23

! Make sure it has units of Pa
estress = dsqrt(estress)

return
end subroutine max_shear_stress

!--------------------------------------------------------------------------------------------------!

subroutine calc_tractions(trac_vector,shear,shearmax,normal,coulomb,stress,trg)
!----
! Calculate various tractions resolved onto a fault geometry from a stress tensor
!----

use trig, only: pi, d2r
implicit none

! Arguments
double precision :: trac_vector(3), shear, shearmax, normal, coulomb, stress(3,3), trg(4)

! Local variables
integer :: i
double precision :: str, dip, rak, fric
double precision :: n(3), s(3), r(3)

str = trg(1)*d2r
dip = trg(2)*d2r
rak = trg(3)*d2r
fric = trg(4)

! Unit normal vector to target plane
n(1) = dsin(dip)*dsin(str+pi/2.0d0)
n(2) = dsin(dip)*dcos(str+pi/2.0d0)
n(3) = dcos(dip)

! Traction vector is stress matrix times normal vector: t = S*n
do i = 1,3
    trac_vector(i) = stress(i,1)*n(1) + stress(i,2)*n(2) + stress(i,3)*n(3)
enddo

! Normal component of traction is parallel to unit normal vector; take dot product
normal = 0.0d0
do i = 1,3
    normal = normal + trac_vector(i)*n(i)
enddo

! Shear component of traction is difference between total traction vector and normal traction
s(1) = trac_vector(1) - normal*n(1)
s(2) = trac_vector(2) - normal*n(2)
s(3) = trac_vector(3) - normal*n(3)
shearmax = dsqrt(s(1)*s(1) + s(2)*s(2) + s(3)*s(3))

! Compute unit slip vector (parallel to rake)
r(1) =  dcos(str-pi/2.0d0)*dcos(rak) + dsin(str-pi/2.0d0)*dsin(rak)*dcos(dip)
r(2) = -dsin(str-pi/2.0d0)*dcos(rak) + dcos(str-pi/2.0d0)*dsin(rak)*dcos(dip)
r(3) =                                                    dsin(rak)*dsin(dip)

! Project shear component on slip vector
shear = s(1)*r(1) + s(2)*r(2) + s(3)*r(3)

! Coulomb stress (recall sign convention: pos = dilation)
coulomb = shear + fric*normal

return
end subroutine calc_tractions

end module elast
