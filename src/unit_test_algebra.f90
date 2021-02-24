program main

use test, only: test_value
use io, only: stdout

use algebra

implicit none

integer :: ierr, nrot
double precision :: vec1(3), vec2(3), vector_6(6), mat1(3,3), mat2(3,3)
double precision :: dot, cross(3)
double precision :: mean_value, stdev_value, coefvar_value

! subroutine normalize()
vec1(1) =   8.6d0
vec1(2) = -14.9d0
vec1(3) =   3.2d0
call normalize(vec1)
call test_value(vec1(1), 0.49146066786720721d0,'normalize(): component 1')
call test_value(vec1(2),-0.85148418037457996d0,'normalize(): component 2')
call test_value(vec1(3), 0.18286908571803062d0,'normalize(): component 3')
write(stdout,*) 'subroutine normalize() passed unit test'
write(stdout,*)

! subroutine dot_product()
vec1(1) =  18.6d0
vec1(2) = -34.7d0
vec1(3) =  -6.6d0
vec2(1) =  15.4d0
vec2(2) =   0.7d0
vec2(3) =   9.0d0
call dot_product(vec1,vec2,dot)
call test_value(dot,202.75000000000003d0,'dot_product()')
write(stdout,*) 'subroutine dot_product() passed unit test'
write(stdout,*)

! subroutine cross_product()
call cross_product(vec1,vec2,cross)
call test_value(cross(1),-307.68000000000001d0,'cross_product(): component 1')
call test_value(cross(2),-269.04000000000002d0,'cross_product(): component 2')
call test_value(cross(3), 547.40000000000009d0,'cross_product(): component 3')
write(stdout,*) 'subroutine cross_product() passed unit test'
write(stdout,*)

! subroutine rotate_vector_angle_axis(vec_in,angle,axis,vec_out,ierr)
vec1(1) = -0.34d0
vec1(2) =  0.75d0
vec1(3) = -0.06d0
call rotate_vector_angle_axis(vec1,  32.0d0,'x',vec2,ierr)
call test_value(vec2(1),               -0.34d0,'rotate_vector_angle_axis(): x-axis: component 1')
call test_value(vec2(2), 0.66783122797131178d0,'rotate_vector_angle_axis(): x-axis: component 2')
call test_value(vec2(3), 0.34655656240551813d0,'rotate_vector_angle_axis(): x-axis: component 3')
call rotate_vector_angle_axis(vec1,-249.2d0,'y',vec2,ierr)
call test_value(vec2(1), 6.4646826635005955d-2,'rotate_vector_angle_axis(): y-axis: component 1')
call test_value(vec2(2),                0.75d0,'rotate_vector_angle_axis(): y-axis: component 2')
call test_value(vec2(3), 0.33914714771913312d0,'rotate_vector_angle_axis(): y-axis: component 3')
call rotate_vector_angle_axis(vec1, 106.1d0,'z',vec2,ierr)
call test_value(vec2(1),-0.62629738349538722d0,'rotate_vector_angle_axis(): z-axis: component 1')
call test_value(vec2(2),-0.53465090239036517d0,'rotate_vector_angle_axis(): z-axis: component 2')
call test_value(vec2(3),               -0.06d0,'rotate_vector_angle_axis(): z-axis: component 3')
write(stdout,*) 'subroutine rotate_vector_angle_axis() passed unit test'
write(stdout,*)

! subroutine rotate_matrix_angle_axis(matrix_in,angle,axis,matrix_out,ierr)
mat1(1,1) =  -4.1065187541649522d0
mat1(1,2) =   1.7393039774189756d0
mat1(1,3) =  -7.5180411049161950d0
mat1(2,1) =   4.2831518746279666d0
mat1(2,2) =   6.9335695426045376d0
mat1(2,3) =  -7.4966824230646907d0
mat1(3,1) =   3.2585177718455043d0
mat1(3,2) =   5.9082031688210321d0
mat1(3,3) = -0.82932751247140502d0
call rotate_matrix_angle_axis(mat1,158.20628169005411d0,'x',mat2,ierr)
call test_value(mat2(1,1),-4.1065187541649522d0,'rotate_matrix_angle_axis(): x-axis: (1,1)')
call test_value(mat2(1,2), 1.1762034338405825d0,'rotate_matrix_angle_axis(): x-axis: (1,2)')
call test_value(mat2(1,3), 7.6264451655602041d0,'rotate_matrix_angle_axis(): x-axis: (1,3)')
call test_value(mat2(2,1),-5.1867971579695755d0,'rotate_matrix_angle_axis(): x-axis: (2,1)')
call test_value(mat2(2,2), 5.3159466459258233d0,'rotate_matrix_angle_axis(): x-axis: (2,2)')
call test_value(mat2(2,3),-9.9538352315436462d0,'rotate_matrix_angle_axis(): x-axis: (2,3)')
call test_value(mat2(3,1),-1.4354313960063656d0,'rotate_matrix_angle_axis(): x-axis: (3,1)')
call test_value(mat2(3,2), 3.4510503603420775d0,'rotate_matrix_angle_axis(): x-axis: (3,2)')
call test_value(mat2(3,3),0.78829538420731038d0,'rotate_matrix_angle_axis(): x-axis: (3,3)')
call rotate_matrix_angle_axis(mat1,12.976505091183620d0,'y',mat2,ierr)
call test_value(mat2(1,1),-4.8733275833542962d0,'rotate_matrix_angle_axis(): y-axis: (1,1)')
call test_value(mat2(1,2), 3.0215818039024152d0,'rotate_matrix_angle_axis(): y-axis: (1,2)')
call test_value(mat2(1,3),-6.5861566427086036d0,'rotate_matrix_angle_axis(): y-axis: (1,3)')
call test_value(mat2(2,1), 2.4903785677930772d0,'rotate_matrix_angle_axis(): y-axis: (2,1)')
call test_value(mat2(2,2), 6.9335695426045376d0,'rotate_matrix_angle_axis(): y-axis: (2,2)')
call test_value(mat2(2,3),-8.2670219500430164d0,'rotate_matrix_angle_axis(): y-axis: (2,3)')
call test_value(mat2(3,1), 4.1904022340530958d0,'rotate_matrix_angle_axis(): y-axis: (3,1)')
call test_value(mat2(3,2), 5.3667575324639163d0,'rotate_matrix_angle_axis(): y-axis: (3,2)')
call test_value(mat2(3,3),-6.2518683282060783d-2,'rotate_matrix_angle_axis(): y-axis: (3,3)')
call rotate_matrix_angle_axis(mat1,296.12107901952908d0,'z',mat2,ierr)
call test_value(mat2(1,1), 7.1742833398756733d0,'rotate_matrix_angle_axis(): z-axis: (1,1)')
call test_value(mat2(1,2), 1.2484031956652555d0,'rotate_matrix_angle_axis(): z-axis: (1,2)')
call test_value(mat2(1,3),-10.040978097812472d0,'rotate_matrix_angle_axis(): z-axis: (1,3)')
call test_value(mat2(2,1), 3.7922510928742463d0,'rotate_matrix_angle_axis(): z-axis: (2,1)')
call test_value(mat2(2,2),-4.3472325514360870d0,'rotate_matrix_angle_axis(): z-axis: (2,2)')
call test_value(mat2(2,3), 3.4496301608646154d0,'rotate_matrix_angle_axis(): z-axis: (2,3)')
call test_value(mat2(3,1), 6.7393988148516133d0,'rotate_matrix_angle_axis(): z-axis: (3,1)')
call test_value(mat2(3,2),-0.32450942648333436d0,'rotate_matrix_angle_axis(): z-axis: (3,2)')
call test_value(mat2(3,3),-0.82932751247140502d0,'rotate_matrix_angle_axis(): z-axis: (3,3)')
write(stdout,*) 'subroutine rotate_matrix_angle_axis() passed unit test'
write(stdout,*)

! subroutine jacobi(a,n,np,d,v,nrot,ierr)
! subroutine eig_sort(d,v,n,np)
mat1(1,1) = -2.2833110157453653d0
mat1(2,2) = -11.608231238033319d0
mat1(3,3) =  8.4575844946838608d0
mat1(1,2) =  6.6226220737607093d0
mat1(1,3) =  6.4092119904680871d0
mat1(2,3) = -4.3740580973518988d0
mat1(2,1) = mat1(1,2)
mat1(3,1) = mat1(1,3)
mat1(3,2) = mat1(2,3)
call jacobi(mat1,3,3,vec1,mat2,nrot,ierr)
call eig_sort(vec1,mat2,3,3)
if (mat2(1,1).gt.0.0d0) then
    mat2(1,1) = -mat2(1,1)
    mat2(2,1) = -mat2(2,1)
    mat2(3,1) = -mat2(3,1)
endif
if (mat2(1,2).lt.0.0d0) then
    mat2(1,2) = -mat2(1,2)
    mat2(2,2) = -mat2(2,2)
    mat2(3,2) = -mat2(3,2)
endif
if (mat2(1,3).lt.0.0d0) then
    mat2(1,3) = -mat2(1,3)
    mat2(2,3) = -mat2(2,3)
    mat2(3,3) = -mat2(3,3)
endif
call test_value(vec1(1),-16.939897642693150d0,'jacobi(): eigenvalue 1')
call test_value(vec1(2),-1.2480176224910310d-2,'jacobi(): eigenvalue 2')
call test_value(vec1(3),11.518420059823237d0,'jacobi(): eigenvalue 3')
call test_value(mat2(1,1),-0.49130341155129770d0,'jacobi(): eigenvector 1, component 1')
call test_value(mat2(2,1), 0.82912602200394458d0,'jacobi(): eigenvector 1, component 2')
call test_value(mat2(3,1), 0.26677892989134439d0,'jacobi(): eigenvector 1, component 3')
call test_value(mat2(1,2), 0.77537716392007883d0,'jacobi(): eigenvector 2, component 1')
call test_value(mat2(2,2), 0.55587251633836765d0,'jacobi(): eigenvector 2, component 2')
call test_value(mat2(3,2),-0.29965980586476104d0,'jacobi(): eigenvector 2, component 3')
call test_value(mat2(1,3), 0.39675081785588218d0,'jacobi(): eigenvector 3, component 1')
call test_value(mat2(2,3),-5.9630405126627532d-2,'jacobi(): eigenvector 3, component 2')
call test_value(mat2(3,3), 0.91598744713839986d0,'jacobi(): eigenvector 3, component 3')
write(stdout,*) 'subroutines jacobi() and eig_sort() passed unit test'
write(stdout,*)

! subroutine mean()
! subroutine stdev()
! subroutine coefvar()
vector_6(1) =  727.7d0
vector_6(2) = 1086.5d0
vector_6(3) = 1091.0d0
vector_6(4) = 1361.3d0
vector_6(5) = 1490.5d0
vector_6(6) = 1956.1d0
mean_value = mean(vector_6,6,ierr)
stdev_value = stdev(vector_6,6,ierr)
coefvar_value = coefvar(vector_6,6,ierr)
call test_value(mean_value,1285.5166666666667d0,'mean')
call test_value(stdev_value,420.96248961952256d0,'stdev')
call test_value(coefvar_value,0.32746560237999445d0,'coefvar')
write(stdout,*) 'functions mean(), stdev(), and coefvar() passed unit test'
write(stdout,*)


write(stdout,*) 'algebra_module unit test passed'
end
