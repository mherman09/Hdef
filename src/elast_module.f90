module elast
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

end module elast
