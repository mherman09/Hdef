module algebra

public :: normalize
public :: dot_product
public :: cross_product
public :: rotate_vector_angle_axis
public :: rotate_matrix_angle_axis

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

subroutine normalize(vector)
implicit none
double precision :: vector(3)
double precision :: magnitude
magnitude = sqrt(vector(1)*vector(1)+vector(2)*vector(2)+vector(3)*vector(3))
vector = vector/magnitude
return
end

!--------------------------------------------------------------------------------------------------!

subroutine dot_product(vec1,vec2,dot)
implicit none
double precision :: vec1(3), vec2(3), dot
dot = vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)
return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine cross_product(vec1,vec2,cross)
implicit none
double precision :: vec1(3), vec2(3), cross(3)
cross(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
cross(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
cross(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine rotate_vector_angle_axis(vec_in,angle,axis,vec_out,ierr)
!----
! Rotate a vector counter-clockwise through a coordinate system by an angle (in degrees) about an
! axis of the coordinate system.
!----

use io, only: stderr
use trig, only: d2r

implicit none

! Arguments
double precision :: vec_in(3), angle, vec_out(3)
character(len=*) :: axis
integer :: ierr

! Local variables
double precision :: vec_tmp(3), rot_matrix(3,3), cosa, sina


! Initialize error
ierr = 0

! Cosine and sine of angle
cosa = cos(angle*d2r)
sina = sin(angle*d2r)

vec_tmp = vec_in

! Build the rotation matrix
rot_matrix = 0.0d0
if (axis.eq.'x') then
    rot_matrix(1,1) = 1.0d0
    rot_matrix(2,2) =  cosa
    rot_matrix(3,3) =  cosa
    rot_matrix(2,3) = -sina
    rot_matrix(3,2) =  sina
elseif (axis.eq.'y') then
    rot_matrix(2,2) = 1.0d0
    rot_matrix(3,3) =  cosa
    rot_matrix(1,1) =  cosa
    rot_matrix(3,1) = -sina
    rot_matrix(1,3) =  sina
elseif (axis.eq.'z') then
    rot_matrix(3,3) = 1.0d0
    rot_matrix(1,1) =  cosa
    rot_matrix(2,2) =  cosa
    rot_matrix(1,2) = -sina
    rot_matrix(2,1) =  sina
else
    write(stderr,*) 'rotate_vector_angle_axis: no axis named "',trim(axis),'"'
    ierr = 1
    return
endif

! Multiply matrix and input vector to get rotated vector
vec_out(1) = rot_matrix(1,1)*vec_tmp(1) + rot_matrix(1,2)*vec_tmp(2) + rot_matrix(1,3)*vec_tmp(3)
vec_out(2) = rot_matrix(2,1)*vec_tmp(1) + rot_matrix(2,2)*vec_tmp(2) + rot_matrix(2,3)*vec_tmp(3)
vec_out(3) = rot_matrix(3,1)*vec_tmp(1) + rot_matrix(3,2)*vec_tmp(2) + rot_matrix(3,3)*vec_tmp(3)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine rotate_matrix_angle_axis(matrix_in,angle,axis,matrix_out,ierr)
!----
! Rotate a matrix counter-clockwise through a coordinate system by an angle (in degrees) about an
! axis of the coordinate system.
!----

use io, only: stderr
use trig, only: d2r

implicit none

! Arguments
double precision :: matrix_in(3,3), angle, matrix_out(3,3)
character(len=*) :: axis
integer :: ierr

! Local variables
double precision :: rot_matrix(3,3), rot_matrix_trans(3,3), matrix_tmp(3,3), cosa, sina
integer :: i, j


! Initialize error
ierr = 0

! Cosine and sine of angle
cosa = cos(angle*d2r)
sina = sin(angle*d2r)

! Build the rotation matrix
rot_matrix = 0.0d0
if (axis.eq.'x') then
    rot_matrix(1,1) = 1.0d0
    rot_matrix(2,2) =  cosa
    rot_matrix(3,3) =  cosa
    rot_matrix(2,3) = -sina
    rot_matrix(3,2) =  sina
elseif (axis.eq.'y') then
    rot_matrix(2,2) = 1.0d0
    rot_matrix(3,3) =  cosa
    rot_matrix(1,1) =  cosa
    rot_matrix(3,1) = -sina
    rot_matrix(1,3) =  sina
elseif (axis.eq.'z') then
    rot_matrix(3,3) = 1.0d0
    rot_matrix(1,1) =  cosa
    rot_matrix(2,2) =  cosa
    rot_matrix(1,2) = -sina
    rot_matrix(2,1) =  sina
else
    write(stderr,*) 'rotate_matrix_angle_axis: no axis named "',trim(axis),'"'
    ierr = 1
    return
endif

! And its transpose (inverse)
rot_matrix_trans = rot_matrix
rot_matrix_trans(1,2) = rot_matrix(2,1)
rot_matrix_trans(1,3) = rot_matrix(3,1)
rot_matrix_trans(2,1) = rot_matrix(1,2)
rot_matrix_trans(2,3) = rot_matrix(3,2)
rot_matrix_trans(3,1) = rot_matrix(1,3)
rot_matrix_trans(3,2) = rot_matrix(2,3)

! Multiply input matrix by rotation matrix and transpose to get rotated matrix: R*M*RT
do i = 1,3
    do j = 1,3
        matrix_tmp(i,j) = rot_matrix(i,1)*matrix_in(1,j) + &
                          rot_matrix(i,2)*matrix_in(2,j) + &
                          rot_matrix(i,3)*matrix_in(3,j)
    enddo
enddo
do i = 1,3
    do j = 1,3
        matrix_out(i,j) = matrix_tmp(i,1)*rot_matrix_trans(1,j) + &
                          matrix_tmp(i,2)*rot_matrix_trans(2,j) + &
                          matrix_tmp(i,3)*rot_matrix_trans(3,j)
    enddo
enddo

return
end subroutine

end module
