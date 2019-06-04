module fitutil

character(len=512) :: input_file
character(len=512) :: output_file

integer :: poly_order
double precision :: period(4)
integer :: nperiod
double precision :: exp_const(4)
integer :: nexponent
logical :: doPolynomial
logical :: doSinusoid
logical :: doExponential
logical :: findExpConstant

logical :: labelCoefficients

character(len=512) :: pre_file
logical :: printPredicted

end module

!==================================================================================================!

program main

use io, only: stdin, stdout, stderr, fileExists, line_count
use trig, only: pi
use solver, only: solve_dgels

use fitutil, only: input_file, &
                   output_file, &
                   poly_order, &
                   period, &
                   nperiod, &
                   exp_const, &
                   nexponent, &
                   doPolynomial, &
                   doSinusoid, &
                   doExponential, &
                   findExpConstant, &
                   labelCoefficients, &
                   pre_file, &
                   printPredicted

implicit none

! Local variables
integer :: i, j, ios, ierr, nrows, ncols, ptr_poly, ptr_sin, ptr_exp, npoly, nsin, nexp, luin, luout
integer :: iperiod, iexponent
integer, parameter :: NROWS_STDIN = 100000
character(len=512) :: line
double precision, allocatable :: obs(:,:), A(:,:), b(:), x(:)
double precision :: pre, amp, pha


#ifndef USE_LAPACK
    call usage('fitutil: does not work if compiled without LAPACK libraries')
#endif


! Parse command line
call gcmdln()

! Make sure something needs to be fit...
if (.not.doPolynomial .and. .not.doSinusoid .and. .not.doExponential) then
    call usage('fitutil: no fitting function specified')
endif


! Determine input/output streams
if (input_file.eq.''.or.input_file.eq.'none'.or.input_file.eq.'stdin') then
    luin = stdin
else
    if (.not.fileExists(input_file)) then
        call usage('fitutil: no file found named "'//trim(input_file)//'"')
    endif
    luin = 11
endif

if (output_file.eq.''.or.output_file.eq.'none'.or.output_file.eq.'stdout') then
    luout = stdout
else
    luout = 12
    open(unit=luout,file=output_file,status='unknown')
endif


! Calculate number of parameters to fit and set up pointers to array
ncols = 0
ptr_poly = 0
ptr_sin = 0
ptr_exp = 0
npoly = 0
nsin = 0
nexp = 0

if (doPolynomial) then
    ptr_poly = ncols + 1
    npoly = poly_order + 1
    ncols = ncols + npoly
endif

if (doSinusoid) then
    ptr_sin = ncols + 1
    nsin = 2*nperiod
    ncols = ncols + nsin
endif

if (doExponential) then
    ptr_exp = ncols + 1
    if (findExpConstant) then
        if (ptr_exp.gt.1) then
            call usage('fitutil: cannot linearize model matrix with exponential and polynomial '//&
                       'or sinusoid')
        endif
        nexp = 2
    else
        nexp = nexponent
    endif
    ncols = ncols + nexp
endif


! Count number of input values and read into obs array
nrows = 0
if (luin.eq.stdin) then
    allocate(obs(NROWS_STDIN,2))
    do
        read(luin,'(A)',iostat=ios) line
        if (ios.ne.0) then
            exit
        else
            nrows = nrows + 1
            if (nrows.gt.NROWS_STDIN) then
                write(stderr,*) 'fitutil: max number of lines using stdin is ',NROWS_STDIN
                call usage(     'Read input from file to increase number of inputs')
            endif
            read(line,*,iostat=ios) obs(nrows,1),obs(nrows,2)
            if (ios.ne.0) then
                write(stderr,*) 'fitutil: error parsing x-y from line ',nrows
                call usage(     'Offending line: '//trim(line))
            endif
        endif
    enddo

else
    nrows = line_count(input_file)
    open(unit=luin,file=input_file,status='old')
    allocate(obs(nrows,2))
    do i = 1,nrows
        read(luin,'(A)',iostat=ios,err=9001,end=9001) line
        read(line,*,iostat=ios) obs(i,1),obs(i,2)
        9001 if (ios.ne.0) then
            write(stderr,*) 'fitutil: error parsing x-y from line ',i
            call usage(     'Offending line: '//trim(line))
        endif
    enddo
endif


! Load model matrix and data vector

! Allocate memory to model matrix, data vector, and solution vector
allocate(A(nrows,ncols))
allocate(b(nrows))
allocate(x(ncols))

! Load arrays
do i = 1,nrows
    if (doPolynomial) then
        do j = 1,poly_order+1
            A(i,j) = obs(i,1)**(j-1)
        enddo
    endif

    if (doSinusoid) then
        do j = 1,nperiod
            A(i,ptr_sin+2*(j-1)) = sin(2.0d0*pi/period(j)*obs(i,1))
            A(i,ptr_sin+2*(j-1)+1) = cos(2.0d0*pi/period(j)*obs(i,1))
        enddo
    endif

    if (doExponential) then
        if (findExpConstant) then
            A(i,ptr_exp) = obs(i,1)
            A(i,ptr_exp+1) = 1.0d0
        else
            do j = 1,nexponent
                A(i,ptr_exp+j-1) = exp(exp_const(j)*obs(i,1))
            enddo
        endif
    endif

    b(i) = obs(i,2)
enddo

if (doExponential.and.findExpConstant) then
    b = log(b)
endif


! Solve system of equations Ax = b for x
call solve_dgels(A,b,x,nrows,ncols,ierr)


! Write results to file or standard output
i = 1
do while (i.le.ncols)
    ! Polynomial
    if (i.le.ptr_poly+npoly-1) then
        if (labelCoefficients) then
            write(luout,'(" Order_",I0.2,X,1PE16.8)') i-1,x(i)
        else
            write(luout,*) x(i)
        endif

    ! Sinusoid
    elseif (i.le.ptr_sin+nsin-1) then
        amp = sqrt(x(i)**2+x(i+1)**2)
        pha = atan2(x(i+1),x(i))
        if (labelCoefficients) then
            write(luout,'(" SinAmplitude_",I0.2,X,1PE16.8)') i-ptr_sin+1,amp
            write(luout,'(" SinPhase_",I0.2,X,1PE16.8)') i-ptr_sin+1,pha
        else
            write(luout,*) amp
            write(luout,*) pha
        endif
        i = i + 1

    ! Exponential
    elseif (i.le.ptr_exp+nexp-1) then
        if (findExpConstant) then
            if (labelCoefficients) then
                write(luout,'(" ExpAmplitude_",I0.2,X,1PE16.8)') i,exp(x(i+1))
                write(luout,'(" ExpConst_",I0.2,X,1PE16.8)') i,x(i)
            else
                write(luout,*) exp(x(i+1))
                write(luout,*) x(i)
            endif
            i = i + 1
        else
            if (labelCoefficients) then
                write(luout,'( " ExpAmplitude_",I0.2,X,1PE16.8)') i-ptr_exp+1,x(i)
            else
                write(luout,*) x(i)
            endif
        endif
    endif
    i = i + 1
enddo


! Print predicted results
if (printPredicted) then
    open(unit=31,file=pre_file,status='unknown')
    do i = 1,nrows
        pre = 0.0d0
        j = 1
        do while (j.le.ncols)
            ! Polynomial
            if (j.le.ptr_poly+npoly-1) then
                pre = pre + x(j)*obs(i,1)**(j-1)

            ! Sinusoid
            elseif (j.le.ptr_sin+nsin-1) then
                iperiod = (j-ptr_sin)/nsin + 1
                amp = sqrt(x(j)**2+x(j+1)**2)
                pha = atan2(x(j+1),x(j))
                pre = pre + amp*sin(2.0d0*pi/period(iperiod)*obs(i,1)+pha)
                j = j + 1

            ! Exponential
            elseif (j.le.ptr_exp+nexp-1) then
                if (findExpConstant) then
                    pre = pre + exp(x(i+1))*exp(x(i)*obs(i,1))
                    j = j + 1
                else
                    iexponent = (j-ptr_exp)/nexponent+1
                    pre = pre + x(j)*exp(exp_const(iexponent)*obs(i,1))
                endif
            endif
            j = j + 1
        enddo
        write(31,*) pre
    enddo
    close(31)
endif


if (luin.ne.stdin) then
    close(11)
endif
if (luout.ne.stdout) then
    close(12)
endif

end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use fitutil, only: input_file, &
                   output_file, &
                   poly_order, &
                   period, &
                   nperiod, &
                   exp_const, &
                   nexponent, &
                   doPolynomial, &
                   doSinusoid, &
                   doExponential, &
                   findExpConstant, &
                   labelCoefficients, &
                   pre_file, &
                   printPredicted


implicit none

! Local variables
integer :: i, ios, narg
character(len=512) :: tag
double precision :: dp


input_file = ''
output_file = ''
poly_order = 0
period = 0.0d0
nperiod = 0
exp_const = 0.0d0
nexponent = 0
doPolynomial = .false.
doSinusoid = .false.
doExponential = .false.
findExpConstant = .false.
labelCoefficients = .true.
pre_file = ''
printPredicted = .false.

! Number of arguments
narg = command_argument_count()
if (narg.eq.0) then
  call usage('')
endif

! Read arguments
i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    if (trim(tag).eq.'-f') then
        i = i + 1
        call get_command_argument(i,input_file)

    elseif (tag.eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file)

    elseif (tag.eq.'-poly') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) poly_order
        doPolynomial = .true.

    elseif (tag.eq.'-sin') then
        i = i + 1
        call get_command_argument(i,tag)
        nperiod = nperiod + 1
        if (nperiod.gt.4) then
            call usage('fitutil: maximum number of sinusoidal periods is 4')
        endif
        read(tag,*) period(nperiod)
        doSinusoid = .true.

    elseif (tag.eq.'-exp') then
        ! Default: fit single exponential amplitude and exponential constant to data
        findExpConstant = .true.

        ! Check if the exponential constant is fixed
        i = i + 1
        if (i.le.narg) then
            call get_command_argument(i,tag)
            read(tag,*,iostat=ios) dp
            if (ios.eq.0) then
                ! If this can be read as a real number, then the exponential constant is fixed
                nexponent = nexponent + 1
                if (nexponent.gt.4) then
                    call usage('fitutil: maximum number of exponential constants is 4')
                endif
                read(tag,*,iostat=ios) exp_const(nexponent)
                findExpConstant = .false.
            else
                i = i - 1
            endif
        endif

        doExponential = .true.

    elseif (tag.eq.'-label_coeff') then
        i = i + 1
        call get_command_argument(i,tag)
        if (tag.eq.'Y'.or.tag.eq.'y'.or.tag.eq.'1'.or.tag.eq.'on') then
            labelCoefficients = .true.
        elseif (tag.eq.'N'.or.tag.eq.'n'.or.tag.eq.'0'.or.tag.eq.'off') then
            labelCoefficients = .false.
        endif

    elseif (tag.eq.'-pre') then
        i = i + 1
        call get_command_argument(i,pre_file)
        printPredicted = .true.

    else
        call usage('fitutil: no option '//trim(tag))
    endif

    i = i + 1

enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine usage(str)

use io, only: stderr

implicit none

character(len=*) :: str

if (str.ne.'') then
    write(stderr,*) trim(str)
    write(stderr,*)
endif

write(stderr,*) 'Usage: fitutil ...options...'
write(stderr,*)
write(stderr,*) '-f FILE           X-Y data to fit (default: stdin)'
write(stderr,*) '-o FILE           Inverted coefficients (default: stdout)'
write(stderr,*) '-poly N           Polynomial fit of order N'
write(stderr,*) '-sin T            Sinusoid fit with period T (repeat up to 4 times)'
write(stderr,*) '-exp [TAU]        Exponential fit (optional: up to 4 fixed exponential constants)'
write(stderr,*) '-label_coeff Y|N  Label output coefficients (default: Y)'
write(stderr,*) '-pre PFILE        Predicted values at input x-coordinates'
write(stderr,*)

stop
end subroutine
