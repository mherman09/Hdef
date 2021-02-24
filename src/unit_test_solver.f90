program main

use test, only: test_value
use io, only: stdout

use solver

implicit none

integer :: i, j, ierr, nrows, ncols, imat5x4(5,4), imat2x3(2,3)
double precision :: vec2(2), vec5(5), vec4_constraint(4)
double precision :: mat2x3(2,3), mat5x2(5,2), mat5x4(5,4)
logical :: isConst(4)

#ifdef USE_LAPACK
double precision :: mat3x3(3,3), vec3(3), mat3x3_b(3,3)
#elif USE_SUPERLU
double precision :: mat3x3(3,3), vec3(3)
#endif


!---- Test: loading 1-D double precision array
vec5 = 0.0d0
vec2(1) =  8.0d0
vec2(2) = -3.0d0
call load_dp_array1(vec5,5,vec2,2,3,'vec2',ierr)
call test_value(vec5(1), 0.0d0,'load_dp_array1: vec5(1)')
call test_value(vec5(2), 0.0d0,'load_dp_array1: vec5(2)')
call test_value(vec5(3), 8.0d0,'load_dp_array1: vec5(3)')
call test_value(vec5(4),-3.0d0,'load_dp_array1: vec5(4)')
call test_value(vec5(5), 0.0d0,'load_dp_array1: vec5(5)')
call load_dp_array1(vec5,5,vec2,2,5,'vec2',ierr)
call test_value(ierr,1,'load_dp_array1: overflow')
write(stdout,*) 'subroutine load_dp_array1() passed unit test'


!---- Test: loading 2-D double precision array
mat5x4 = 0.0d0
do i = 1,2
    do j = 1,3
        mat2x3(i,j) = dble(i)*10.0d0 + dble(j)
    enddo
enddo
call load_dp_array2(mat5x4,5,4,mat2x3,2,3,3,2,'mat2x3',ierr)
call test_value(mat5x4(1,1), 0.0d0,'load_dp_array2: mat5x4(1,1)')
call test_value(mat5x4(1,2), 0.0d0,'load_dp_array2: mat5x4(1,2)')
call test_value(mat5x4(1,3), 0.0d0,'load_dp_array2: mat5x4(1,3)')
call test_value(mat5x4(1,4), 0.0d0,'load_dp_array2: mat5x4(1,4)')
call test_value(mat5x4(2,1), 0.0d0,'load_dp_array2: mat5x4(2,1)')
call test_value(mat5x4(2,2), 0.0d0,'load_dp_array2: mat5x4(2,2)')
call test_value(mat5x4(2,3), 0.0d0,'load_dp_array2: mat5x4(2,3)')
call test_value(mat5x4(2,4), 0.0d0,'load_dp_array2: mat5x4(2,4)')
call test_value(mat5x4(3,1), 0.0d0,'load_dp_array2: mat5x4(3,1)')
call test_value(mat5x4(3,2),11.0d0,'load_dp_array2: mat5x4(3,2)')
call test_value(mat5x4(3,3),12.0d0,'load_dp_array2: mat5x4(3,3)')
call test_value(mat5x4(3,4),13.0d0,'load_dp_array2: mat5x4(3,4)')
call test_value(mat5x4(4,1), 0.0d0,'load_dp_array2: mat5x4(4,1)')
call test_value(mat5x4(4,2),21.0d0,'load_dp_array2: mat5x4(4,2)')
call test_value(mat5x4(4,3),22.0d0,'load_dp_array2: mat5x4(4,3)')
call test_value(mat5x4(4,4),23.0d0,'load_dp_array2: mat5x4(4,4)')
call test_value(mat5x4(5,1), 0.0d0,'load_dp_array2: mat5x4(5,1)')
call test_value(mat5x4(5,2), 0.0d0,'load_dp_array2: mat5x4(5,2)')
call test_value(mat5x4(5,3), 0.0d0,'load_dp_array2: mat5x4(5,3)')
call test_value(mat5x4(5,4), 0.0d0,'load_dp_array2: mat5x4(5,4)')
call load_dp_array2(mat5x4,5,4,mat2x3,2,3,3,3,'mat2x3',ierr)
call test_value(ierr,1,'load_dp_array2: overflow')
write(stdout,*) 'subroutine load_dp_array2() passed unit test'


!---- Test: loading 2-D integer array
imat5x4 = 0
do i = 1,2
    do j = 1,3
        imat2x3(i,j) = i*10 + j
    enddo
enddo
call load_int_array2(imat5x4,5,4,imat2x3,2,3,3,2,'imat2x3',ierr)
call test_value(imat5x4(1,1), 0,'load_int_array2: imat5x4(1,1)')
call test_value(imat5x4(1,2), 0,'load_int_array2: imat5x4(1,2)')
call test_value(imat5x4(1,3), 0,'load_int_array2: imat5x4(1,3)')
call test_value(imat5x4(1,4), 0,'load_int_array2: imat5x4(1,4)')
call test_value(imat5x4(2,1), 0,'load_int_array2: imat5x4(2,1)')
call test_value(imat5x4(2,2), 0,'load_int_array2: imat5x4(2,2)')
call test_value(imat5x4(2,3), 0,'load_int_array2: imat5x4(2,3)')
call test_value(imat5x4(2,4), 0,'load_int_array2: imat5x4(2,4)')
call test_value(imat5x4(3,1), 0,'load_int_array2: imat5x4(3,1)')
call test_value(imat5x4(3,2),11,'load_int_array2: imat5x4(3,2)')
call test_value(imat5x4(3,3),12,'load_int_array2: imat5x4(3,3)')
call test_value(imat5x4(3,4),13,'load_int_array2: imat5x4(3,4)')
call test_value(imat5x4(4,1), 0,'load_int_array2: imat5x4(4,1)')
call test_value(imat5x4(4,2),21,'load_int_array2: imat5x4(4,2)')
call test_value(imat5x4(4,3),22,'load_int_array2: imat5x4(4,3)')
call test_value(imat5x4(4,4),23,'load_int_array2: imat5x4(4,4)')
call test_value(imat5x4(5,1), 0,'load_int_array2: imat5x4(5,1)')
call test_value(imat5x4(5,2), 0,'load_int_array2: imat5x4(5,2)')
call test_value(imat5x4(5,3), 0,'load_int_array2: imat5x4(5,3)')
call test_value(imat5x4(5,4), 0,'load_int_array2: imat5x4(5,4)')
call load_int_array2(imat5x4,5,4,imat2x3,2,3,3,3,'imat2x3',ierr)
call test_value(ierr,1,'load_int_array2: overflow')
write(stdout,*) 'subroutine load_int_array2() passed unit test'


!---- Test: load equation constraints
nrows = 5
ncols = 4
vec5 = 0.0d0
vec4_constraint(1) = 0.0d0
vec4_constraint(2) = 1.0d0
vec4_constraint(3) = 0.0d0
vec4_constraint(4) = 2.0d0
isConst(1) = .false.
isConst(2) = .true.
isConst(3) = .false.
isConst(4) = .true.
! do i = 1,5
!     write(0,*) mat5x4(i,:)
! enddo
call load_dp_array_constraints(mat5x4,vec5,nrows,ncols,vec4_constraint,isConst,ierr)
! do i = 1,5
!     write(0,*) mat5x4(i,:)
! enddo
! do i = 1,5
!     write(0,*) vec5(i)
! enddo
call test_value(mat5x4(1,1), 0.0d0,'load_dp_array_constraints: mat5x4(1,1)')
call test_value(mat5x4(1,2), 0.0d0,'load_dp_array_constraints: mat5x4(1,2)')
call test_value(mat5x4(2,1), 0.0d0,'load_dp_array_constraints: mat5x4(2,1)')
call test_value(mat5x4(2,2), 0.0d0,'load_dp_array_constraints: mat5x4(2,2)')
call test_value(mat5x4(3,1), 0.0d0,'load_dp_array_constraints: mat5x4(3,1)')
call test_value(mat5x4(3,2),12.0d0,'load_dp_array_constraints: mat5x4(3,2)')
call test_value(mat5x4(4,1), 0.0d0,'load_dp_array_constraints: mat5x4(4,1)')
call test_value(mat5x4(4,2),22.0d0,'load_dp_array_constraints: mat5x4(4,2)')
call test_value(mat5x4(5,1), 0.0d0,'load_dp_array_constraints: mat5x4(5,1)')
call test_value(mat5x4(5,2), 0.0d0,'load_dp_array_constraints: mat5x4(5,2)')
call test_value(vec5(1),  0.0d0,'load_dp_array_constraints: vec5(1)')
call test_value(vec5(2),  0.0d0,'load_dp_array_constraints: vec5(2)')
call test_value(vec5(3),-37.0d0,'load_dp_array_constraints: vec5(3)')
call test_value(vec5(4),-67.0d0,'load_dp_array_constraints: vec5(4)')
call test_value(vec5(5),  0.0d0,'load_dp_array_constraints: vec5(5)')
call test_value(ncols,2,'load_dp_array_constraints: ncols output')
write(stdout,*) 'subroutine load_dp_array_constraints() passed unit test'


!---- Test: linear equation solver
! dgesv
#ifdef USE_LAPACK
mat3x3(1,1) = 1.0d0
mat3x3(1,2) = 2.0d0
mat3x3(1,3) = 3.0d0
mat3x3(2,1) = 4.0d0
mat3x3(2,2) = 5.0d0
mat3x3(2,3) = 4.0d0
mat3x3(3,1) = 3.0d0
mat3x3(3,2) = 2.0d0
mat3x3(3,3) = 1.0d0
vec3(1) = 11.0d0
vec3(2) = 12.0d0
vec3(3) = 13.0d0
call solve_dgesv(mat3x3,vec3,vec3,3,ierr)
! do i = 1,3
!     write(0,*) vec3(i)
! enddo
call test_value(vec3(1),  9.5d0,'solve_dgesv: vec3(1)')
call test_value(vec3(2),-12.0d0,'solve_dgesv: vec3(2)')
call test_value(vec3(3),  8.5d0,'solve_dgesv: vec3(3)')
write(stdout,*) 'subroutine solve_dgesv() (linear equation solver) passed unit test'
mat3x3_b(1,1) = -7.0d0
mat3x3_b(1,2) =  4.0d0
mat3x3_b(1,3) =  1.0d0
mat3x3_b(2,1) = -3.0d0
mat3x3_b(2,2) = -2.0d0
mat3x3_b(2,3) =  8.0d0
mat3x3_b(3,1) =  4.0d0
mat3x3_b(3,2) = -4.0d0
mat3x3_b(3,3) =  0.0d0
call solve_dgesv_nrhs(mat3x3,mat3x3_b,mat3x3_b,3,3,ierr)
call test_value(mat3x3_b(1,1), 2.375d0,'solve_dgesv_nrhs: mat3x3_b(1,1)')
call test_value(mat3x3_b(1,2),-1.000d0,'solve_dgesv_nrhs: mat3x3_b(1,2)')
call test_value(mat3x3_b(1,3),-3.625d0,'solve_dgesv_nrhs: mat3x3_b(1,3)')
call test_value(mat3x3_b(2,1), 0.000d0,'solve_dgesv_nrhs: mat3x3_b(2,1)')
call test_value(mat3x3_b(2,2),-2.000d0,'solve_dgesv_nrhs: mat3x3_b(2,1)')
call test_value(mat3x3_b(2,3), 7.000d0,'solve_dgesv_nrhs: mat3x3_b(2,1)')
call test_value(mat3x3_b(3,1),-3.125d0,'solve_dgesv_nrhs: mat3x3_b(3,1)')
call test_value(mat3x3_b(3,2), 3.000d0,'solve_dgesv_nrhs: mat3x3_b(3,2)')
call test_value(mat3x3_b(3,3),-3.125d0,'solve_dgesv_nrhs: mat3x3_b(3,3)')

! dsysv
mat3x3(1,1) = -2.0d0
mat3x3(1,2) =  1.0d0
mat3x3(1,3) =  3.0d0
mat3x3(2,1) =  1.0d0
mat3x3(2,2) =  1.0d0
mat3x3(2,3) = -1.0d0
mat3x3(3,1) =  3.0d0
mat3x3(3,2) = -1.0d0
mat3x3(3,3) =  2.0d0
vec3(1) = 5.0d0
vec3(2) = 9.0d0
vec3(3) = 2.0d0
call solve_dsysv(mat3x3,vec3,vec3,3,ierr)
call test_value(vec3(1),2.526315789473684d0,'solve_dsysv: vec3(1)')
call test_value(vec3(2),7.368421052631579d0,'solve_dsysv: vec3(2)')
call test_value(vec3(3),0.894736842105263d0,'solve_dsysv: vec3(3)')
mat3x3_b(1,1) = -7.0d0
mat3x3_b(1,2) =  4.0d0
mat3x3_b(1,3) =  1.0d0
mat3x3_b(2,1) = -3.0d0
mat3x3_b(2,2) = -2.0d0
mat3x3_b(2,3) =  8.0d0
mat3x3_b(3,1) =  4.0d0
mat3x3_b(3,2) = -4.0d0
mat3x3_b(3,3) =  0.0d0
call solve_dsysv_nrhs(mat3x3,mat3x3_b,mat3x3_b,3,3,ierr)
call test_value(mat3x3_b(1,1),0.42105263157894735d0,'solve_dsysv_nrhs: mat3x3_b(1,1)')
call test_value(mat3x3_b(1,2),-1.5789473684210527d0,'solve_dsysv_nrhs: mat3x3_b(1,2)')
call test_value(mat3x3_b(1,3),2.0526315789473681d0,'solve_dsysv_nrhs: mat3x3_b(1,3)')
call test_value(mat3x3_b(2,1),-4.1052631578947363d0,'solve_dsysv_nrhs: mat3x3_b(2,1)')
call test_value(mat3x3_b(2,2),-0.10526315789473685d0,'solve_dsysv_nrhs: mat3x3_b(2,1)')
call test_value(mat3x3_b(2,3),5.7368421052631575d0,'solve_dsysv_nrhs: mat3x3_b(2,1)')
call test_value(mat3x3_b(3,1),-0.68421052631578960d0,'solve_dsysv_nrhs: mat3x3_b(3,1)')
call test_value(mat3x3_b(3,2),0.31578947368421056d0,'solve_dsysv_nrhs: mat3x3_b(3,2)')
call test_value(mat3x3_b(3,3),-0.21052631578947367d0,'solve_dsysv_nrhs: mat3x3_b(3,3)')
#endif

!---- Test: sparse linear equation solver
#ifdef USE_SUPERLU
mat3x3(1,1) = 1.0d0
mat3x3(1,2) = 0.0d0
mat3x3(1,3) = 3.0d0
mat3x3(2,1) = 0.0d0
mat3x3(2,2) = 5.0d0
mat3x3(2,3) = 0.0d0
mat3x3(3,1) = 3.0d0
mat3x3(3,2) = 0.0d0
mat3x3(3,3) = 1.0d0
vec3(1) = 11.0d0
vec3(2) = 12.0d0
vec3(3) = 13.0d0
call solve_dgssv(mat3x3,vec3,vec3,3,3,ierr)
call test_value(vec3(1),3.5d0,'solve_dgssv: vec3(1)')
call test_value(vec3(2),2.4d0,'solve_dgssv: vec3(2)')
call test_value(vec3(3),2.5d0,'solve_dgssv: vec3(3)')
! do i = 1,3
!     write(0,*) vec3(i)
! enddo
write(stdout,*) 'subroutine solve_dgssv() (sparse linear equation solver) passed unit test'
#endif

!---- Test: least squares solver
#ifdef USE_LAPACK
mat5x2(:,1) = 1.0d0
mat5x2(1,2) = 3.0d0
mat5x2(2,2) = 4.0d0
mat5x2(3,2) = 5.0d0
mat5x2(4,2) = 6.0d0
mat5x2(5,2) = 7.0d0
vec5(1) = 2.1d0
vec5(2) = 4.1d0
vec5(3) = 5.8d0
vec5(4) = 7.9d0
vec5(5) = 10.1d0
call solve_dgels(mat5x2,vec5,vec2,5,2,ierr)
! do i = 1,2
!     write(0,*) vec2(i)
! enddo
call test_value(vec2(1),-3.90d0,'solve_dgels: vec2(1)')
call test_value(vec2(2), 1.98d0,'solve_dgels: vec2(2)')
write(stdout,*) 'subroutine solve_dgels() (least-squares solver) passed unit test'
#endif

!---- Test: non-negative least squares solver
mat5x2(:,1) = 1.0d0
mat5x2(1,2) = 3.0d0
mat5x2(2,2) = 4.0d0
mat5x2(3,2) = 5.0d0
mat5x2(4,2) = 6.0d0
mat5x2(5,2) = 7.0d0
vec5(1) = 10.1d0
vec5(2) = 12.1d0
vec5(3) = 13.8d0
vec5(4) = 15.9d0
vec5(5) = 18.1d0
call solve_nnls(mat5x2,vec5,vec2,5,2,ierr)
! do i = 1,2
!     write(0,*) vec2(i)
! enddo
call test_value(vec2(1), 4.10d0,'solve_nnls: vec2(1)')
call test_value(vec2(2), 1.98d0,'solve_nnls: vec2(2)')
write(stdout,*) 'subroutine solve_nnls() (non-negative least-squares solver) passed unit test'



write(stdout,*) 'solver_module unit test passed'
end
