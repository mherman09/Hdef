module algebra

public :: normalize
public :: dot_product
public :: cross_product
public :: rotate_vector_angle_axis
public :: rotate_matrix_angle_axis
public :: jacobi
public :: eig_sort
public :: mean
public :: stdev
public :: coefvar

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

subroutine normalize(vector)
implicit none
double precision, intent(inout) :: vector(3)
double precision :: magnitude
magnitude = sqrt(vector(1)*vector(1)+vector(2)*vector(2)+vector(3)*vector(3))
vector = vector/magnitude
return
end

!--------------------------------------------------------------------------------------------------!

subroutine dot_product(vec1,vec2,dot)
implicit none
double precision, intent(in)  :: vec1(3), vec2(3)
double precision, intent(out) :: dot
dot = vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)
return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine cross_product(vec1,vec2,cross)
implicit none
double precision, intent(in)  :: vec1(3), vec2(3)
double precision, intent(out) :: cross(3)
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
double precision, intent(in)  :: vec_in(3), angle
double precision, intent(out) :: vec_out(3)
character(len=*), intent(in)  :: axis
integer, intent(out)          :: ierr

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
double precision, intent(in)  :: matrix_in(3,3), angle
double precision, intent(out) :: matrix_out(3,3)
character(len=*), intent(in)  :: axis
integer, intent(out)          :: ierr

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

!--------------------------------------------------------------------------------------------------!

subroutine jacobi(a,n,np,d,v,nrot,ierr)
!----
! Computes all eigenvalues and eigenvectors of a real symmetric matrix a, which is of size n by n,
! stored in a physical np by np array. On output, elements of a above the diagonal are destroyed. d
! returns the eigenvalues of a in its first n elements. v is a matrix with the same logical and
! physical dimensions as a, whose columns contain, on output, the normalized eigenvectors of a.
! nrot returns the number of Jacobi rotations that were required.
!
! Numerical Recipes (Press et al.)
!----

use io, only: stderr

implicit none

! Arguments
integer, intent(in)  :: n, np
integer, intent(out) :: nrot, ierr
double precision :: a(np,np), d(np), v(np,np)

! Local variables
integer, parameter :: NMAX = 500
integer :: i, ip, iq, j
double precision :: c, g, h, s, sm, t, tau, theta, tresh, b(NMAX), z(NMAX)
double precision, parameter :: eps = 1.0d-12


ierr = 0

! Initialize v to the identity matrix
v = 0.0d0
do ip = 1,n
    v(ip,ip) = 1.0d0
enddo

! Initialize b and d to the diagonal of a
! Vector z will accumulate terms of the form tapq as in equation (11.1.14)
do ip = 1,n
    b(ip) = a(ip,ip)
    d(ip) = b(ip)
    z(ip) = 0.0d0
enddo

nrot = 0
do i = 1,50

    ! Sum off-diagonal elements
    sm = 0.0d0
    do ip = 1,n-1
        do iq = ip+1,n
            sm = sm + abs(a(ip,iq))
        enddo
    enddo

    ! The normal return, which relies on quadratic convergence
    if (abs(sm).lt.eps) then
        return
    endif

    if (i.lt.4) then
        tresh = 0.2d0*sm/n**2
    else
        tresh = 0.0d0
    endif

    do ip = 1,n-1
        do iq = ip+1,n

            g = 100.0d0*abs(a(ip,iq))

            ! After four sweeps, skip the rotation if the off-diagonal element is small
            if ((i.gt.4) .and. abs(g).lt.eps) then
                a(ip,iq) = 0.0d0
            elseif (abs(a(ip,iq)).gt.tresh) then
                h = d(iq)-d(ip)
                if (abs(g).lt.eps) then
                    t = a(ip,iq)/h  ! t = 1/(2theta)
                else
                    theta = 0.5d0*h/a(ip,iq) !Equation (11.1.10)
                    t = 1.0d0/(abs(theta)+sqrt(1.0d0+theta**2))
                    if (theta.lt.0.0d0) then
                        t = -t
                    endif
                endif
                c = 1.0d0/sqrt(1.0d0+t**2)
                s = t*c
                tau = s/(1.0d0+c)
                h = t*a(ip,iq)
                z(ip) = z(ip)-h
                z(iq) = z(iq)+h
                d(ip) = d(ip)-h
                d(iq) = d(iq)+h
                a(ip,iq) = 0.0d0
                do j = 1,ip-1 ! Case of rotations 1 < j < p
                    g = a(j,ip)
                    h = a(j,iq)
                    a(j,ip) = g-s*(h+g*tau)
                    a(j,iq) = h+s*(g-h*tau)
                enddo
                do j = ip+1,iq-1 !Case of rotations p < j < q
                    g = a(ip,j)
                    h = a(j,iq)
                    a(ip,j) = g-s*(h+g*tau)
                    a(j,iq) = h+s*(g-h*tau)
                enddo
                do j = iq+1,n !Case of rotations q < j ≤ n
                    g = a(ip,j)
                    h = a(iq,j)
                    a(ip,j) = g-s*(h+g*tau)
                    a(iq,j) = h+s*(g-h*tau)
                enddo
                do j = 1,n
                    g = v(j,ip)
                    h = v(j,iq)
                    v(j,ip) = g-s*(h+g*tau)
                    v(j,iq) = h+s*(g-h*tau)
                enddo
                nrot = nrot + 1
            endif
        enddo
    enddo

    ! Update d with the sum of tapq, and reinitialize z
    do ip = 1,n
        b(ip) = b(ip)+z(ip)
        d(ip) = b(ip)
        z(ip) = 0.0d0
    enddo
enddo

write(stderr,*) 'Warning: jacobi: too many iterations'
ierr = 1

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine eig_sort(d,v,n,np)
!----
! Given the eigenvalues d and eigenvectors v as output from jacobi, this routine sorts the
! eigenvalues into ascending order, and rearranges the columns of v correspondingly. The method
! is straight insertion.
!
! Numerical Recipes (Press et al.)
!----

implicit none

! Arguments
integer, intent(in)             :: n, np
double precision, intent(inout) :: d(np), v(np,np)

! local variables
integer :: i, j, k
double precision :: p

do i = 1,n-1
    k = i
    p = d(i)
    do j = i+1,n
        if (d(j).le.p) then
            k = j
            p = d(j)
        endif
    enddo
    if (k.ne.i) then
        d(k) = d(i)
        d(i) = p
        do j = 1,n
            p = v(j,i)
            v(j,i) = v(j,k)
            v(j,k) = p
        enddo
    endif
enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

function mean(values,nvalues,ierr)
!----
! Calculate mean of a set of values
!----

implicit none

! Arguments
integer, intent(in)          :: nvalues
integer, intent(out)         :: ierr
double precision, intent(in) :: values(nvalues)
double precision             :: mean

! Local variables
integer :: i

! Initialize variables
ierr = 0
mean = 0.0d0

! Calculate sum of values
do i = 1,nvalues
    mean = mean + values(i)
enddo

! Divide by number of values
mean = mean/dble(nvalues)

return
end function

!--------------------------------------------------------------------------------------------------!

function stdev(values,nvalues,ierr)
!----
! Calculate standard deviation of a set of values
!----

implicit none

! Arguments
integer, intent(in)          :: nvalues
integer, intent(out)         :: ierr
double precision, intent(in) :: values(nvalues)
double precision             :: stdev

! Local variables
integer :: i
double precision :: m

! Initialize variables
ierr = 0
stdev = 0.0d0

! Calculate mean
m = mean(values,nvalues,ierr)

! Calculate sum of differences squared
do i = 1,nvalues
    stdev = stdev + (values(i)-m)**2
enddo

! Divide by number of values, take square root
stdev = sqrt(stdev/dble(nvalues-1))

return
end function

!--------------------------------------------------------------------------------------------------!

function coefvar(values,nvalues,ierr)
!----
! Calculate coefficient of variation of a set of values
!----

implicit none

! Arguments
integer, intent(in)          :: nvalues
integer, intent(out)         :: ierr
double precision, intent(in) :: values(nvalues)
double precision             :: coefvar

! Local variables
double precision :: s, m

! Initialize variables
ierr = 0
coefvar = 0.0d0

! Calculate mean
m = mean(values,nvalues,ierr)

! Calculate standard deviation
s = stdev(values,nvalues,ierr)

! Divide standard vevation by mean
coefvar = s/m

return
end function

end module
