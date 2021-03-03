module grid

character(len=512) :: output_file                 ! Name of output file (default: stdout)

character(len=32) :: grid_mode                    ! uniform, dipping_plane, cross_section
character(len=32) :: spacing_mode                 ! linear, exponential
character(len=32) :: coord_type                   ! cartesian, geographic

double precision :: x_beg
double precision :: x_end
double precision :: x_inc
integer :: nx
logical :: isXDefined

double precision :: y_beg
double precision :: y_end
double precision :: y_inc
integer :: ny
logical :: isYDefined

double precision :: z_beg
double precision :: z_end
double precision :: z_inc
integer :: nz
logical :: isZDefined

double precision :: x0                            ! x-coordinate of plane origin
double precision :: y0                            ! y-coordinate of plane origin
double precision :: z0                            ! z-coordinate of plane origin
double precision :: strike                        ! strike of plane
double precision :: dip                           ! dip of plane
double precision :: azimuth                       ! azimuth of plane

end module

!==================================================================================================!

program main
!----
! Generate Cartesian or geographic grids
!----

use io, only: stdout, stderr, verbosity
use trig, only: d2r
use geom, only: distaz2lola, lola2distaz
use earth, only: radius_earth_km

use grid, only: output_file, &
                grid_mode, &
                spacing_mode, &
                coord_type, &
                x_beg, &
                x_end, &
                x_inc, &
                nx, &
                isXDefined, &
                y_beg, &
                y_end, &
                y_inc, &
                ny, &
                isYDefined, &
                z_beg, &
                z_end, &
                z_inc, &
                nz, &
                isZDefined, &
                x0, &
                y0, &
                z0, &
                strike, &
                dip, &
                azimuth

implicit none

! Local variables
integer :: i, j, k, lu_out, ierr
double precision :: x, y, z, dx, dy, ddip, dist, angle, tand

! Parse command line
call gcmdln()
if (verbosity.eq.1) then
    write(stdout,'(A,A)')       ' output_file:  ',trim(output_file)
    write(stdout,'(A,A)')       ' grid_mode:    ',trim(grid_mode)
    write(stdout,'(A,A)')       ' spacing_mode: ',trim(spacing_mode)
    write(stdout,'(A,A)')       ' coord_type:   ',trim(coord_type)
    write(stdout,'(A,1PE14.6)') ' x_beg:        ',x_beg
    write(stdout,'(A,1PE14.6)') ' x_end:        ',x_end
    write(stdout,'(A,1PE14.6)') ' x_inc:        ',x_inc
    write(stdout,'(A,I14)')     ' nx:           ',nx
    write(stdout,'(A,L14)')     ' isXDefined:   ',isXDefined
    write(stdout,'(A,1PE14.6)') ' y_beg:        ',y_beg
    write(stdout,'(A,1PE14.6)') ' y_end:        ',y_end
    write(stdout,'(A,1PE14.6)') ' y_inc:        ',y_inc
    write(stdout,'(A,I14)')     ' ny:           ',ny
    write(stdout,'(A,L14)')     ' isYDefined:   ',isYDefined
    write(stdout,'(A,1PE14.6)') ' z_beg:        ',z_beg
    write(stdout,'(A,1PE14.6)') ' z_end:        ',z_end
    write(stdout,'(A,1PE14.6)') ' z_inc:        ',z_inc
    write(stdout,'(A,I14)')     ' nz:           ',nz
    write(stdout,'(A,L14)')     ' isZDefined:   ',isZDefined
    write(stdout,'(A,1PE14.6)') ' x0:           ',x0
    write(stdout,'(A,1PE14.6)') ' y0:           ',y0
    write(stdout,'(A,1PE14.6)') ' z0:           ',z0
    write(stdout,'(A,1PE14.6)') ' strike:       ',strike
    write(stdout,'(A,1PE14.6)') ' dip:          ',dip
    write(stdout,'(A,1PE14.6)') ' azimuth:      ',azimuth
    write(stdout,*)
endif


! Everything is hunky dory to start out
ierr = 0


! Determine the number of points and spacing between points in each dimension

! Calculate nx and dx for first dimension (x)
if (isXDefined) then
    call get_spacing(nx,x_inc,x_beg,x_end,spacing_mode,'x')
else
    call usage('grid: x values are required')
endif

! Calculate ny and dy for second dimension (y), if necessary
if (isYDefined) then
    call get_spacing(ny,y_inc,y_beg,y_end,spacing_mode,'y')
endif

! Calculate nz and dz for third dimension (z), if necessary
if (isZDefined.and.(grid_mode.eq.'uniform'.or.grid_mode.eq.'cross_section')) then
    call get_spacing(nz,z_inc,z_beg,z_end,spacing_mode,'z')
endif


! Where to send output?
if (output_file.eq.'stdout') then
    lu_out = stdout ! print to terminal
else
    lu_out = 33 ! write to file
    open(unit=lu_out,file=output_file,status='unknown')
endif


! Output depends on which dimensions have been specified and the grid mode

if (isXDefined.and..not.isYDefined.and..not.isZDefined) then

    ! Only x coordinate is defined

    if (grid_mode.eq.'uniform') then
        do i = 1,nx+1
            call get_value(x,i,x_inc,x_beg,spacing_mode,'x')
            write(lu_out,*) x
        enddo
    else
        ierr = 1
    endif

elseif (isXDefined.and.isYDefined.and..not.isZDefined) then

    ! X and y coordinates are defined

    if (grid_mode.eq.'uniform') then
        ! Print only x and y
        do i = 1,nx+1
            call get_value(x,i,x_inc,x_beg,spacing_mode,'x')
            do j = 1,ny+1
                call get_value(y,j,y_inc,y_beg,spacing_mode,'y')
                write(lu_out,*) x, y
            enddo
        enddo

    elseif (grid_mode.eq.'dipping_plane') then
        tand = dtan(dip*d2r)

        do i = 1,nx+1
            call get_value(x,i,x_inc,x_beg,spacing_mode,'x')
            do j = 1,ny+1
                call get_value(y,j,y_inc,y_beg,spacing_mode,'y')

                ! Calculate z on dipping plane
                if (coord_type.eq.'cartesian') then
                    dx = x-x0
                    dy = y-y0
                    angle = atan2(dx,dy)                 ! angle clockwise from +y axis
                    dist = sqrt(dx*dx+dy*dy)             ! distance from origin
                elseif (coord_type.eq.'geographic') then
                    call lola2distaz(x0,y0,x,y,dist,angle,'radians','radians',ierr)
                    if (ierr.ne.0) then
                        call usage('grid: error calculating dist and az')
                    endif
                    dist = dist*radius_earth_km
                    angle = angle
                endif

                ddip = dist*dsin(angle-strike*d2r)   ! horizontal, down-dip distance
                if (z0.lt.0.0d0) then
                    ! If input z0 is negative, interpret it as a depth with positive up
                    z = z0-ddip*tand
                else
                    ! If input z0 is positive or zero, interpret it as a depth with positive down
                    z = z0+ddip*tand
                endif

                write(lu_out,*) x, y, z
            enddo
        enddo
    else
        ierr = 1
    endif

elseif (isXDefined.and.isYDefined.and.isZDefined) then

    ! X, y, and z coordinates are defined

    if (grid_mode.eq.'uniform') then

        ! Print x, y, and z
        do i = 1,nx+1
            call get_value(x,i,x_inc,x_beg,spacing_mode,'x')
            do j = 1,ny+1
                call get_value(y,j,y_inc,y_beg,spacing_mode,'y')
                do k = 1,nz+1
                    call get_value(z,k,z_inc,z_beg,spacing_mode,'z')
                    write(lu_out,*) x, y, z
                enddo
            enddo
        enddo


    elseif (grid_mode.eq.'cross_section') then
        do i = 1,nx+1
            call get_value(dx,i,x_inc,x_beg,spacing_mode,'x') ! distance parallel to cross-section plane
            if (coord_type.eq.'cartesian') then
                x = x0 + dx*dsin(azimuth*d2r) ! move global x along plane
                y = y0 + dx*dcos(azimuth*d2r) ! move global y along plane
            elseif (coord_type.eq.'geographic') then
                call distaz2lola(x0,y0,dx/radius_earth_km,azimuth,x,y,'radians','degrees',ierr)
                if (ierr.ne.0) then
                    call usage('grid: error computing longitude and latitude')
                endif
            endif

            do j = 1,ny+1
                call get_value(dy,j,y_inc,y_beg,spacing_mode,'y') ! distance perpendicular to cross-section plane
                if (coord_type.eq.'cartesian') then
                    x = x + dy*dsin((azimuth-90.0d0)*d2r) ! add to global x
                    y = y + dy*dcos((azimuth-90.0d0)*d2r) ! add to global y
                elseif (coord_type.eq.'geographic') then
                    call distaz2lola(x,y,dy/radius_earth_km,azimuth-90.0d0,x,y,'radians','degrees',&
                                     ierr)
                    if (ierr.ne.0) then
                        call usage('grid: error computing longitude and latitude')
                    endif
                endif

                do k = 1,nz+1
                    call get_value(z,k,z_inc,z_beg,spacing_mode,'z')
                    write(lu_out,*) x, y, z, dx, dy
                enddo
            enddo
        enddo

    else
        ierr = 1
    endif

elseif (isXDefined.and..not.isYDefined.and.isZDefined) then

    ! X and z coordinates are defined

    if (grid_mode.eq.'cross_section') then

        dy = 0.0d0

        do i = 1,nx+1
            call get_value(dx,i,x_inc,x_beg,spacing_mode,'x')  ! distance parallel to cross-section plane
            if (coord_type.eq.'cartesian') then
                x = x0 + dx*dsin(azimuth*d2r) ! move global x along plane
                y = y0 + dx*dcos(azimuth*d2r) ! move global y along plane
            elseif (coord_type.eq.'geographic') then
                call distaz2lola(x0,y0,dx/radius_earth_km,azimuth,x,y,'radians','degrees',ierr)
                if (ierr.ne.0) then
                    call usage('grid: error computing longitude and latitude')
                endif
            endif

            do k = 1,nz+1
                call get_value(z,k,z_inc,z_beg,spacing_mode,'z')
                write(lu_out,*) x, y, z, dx, dy
            enddo
        enddo

    else
        ierr = 1
    endif

else
    write(stderr,'(A,L14)')     ' isXDefined:  ',isXDefined
    write(stderr,'(A,L14)')     ' isYDefined:  ',isYDefined
    write(stderr,'(A,L14)')     ' isZDefined:  ',isZDefined
    call usage('grid: WHOOPS! No output defined for this combination of inputs')
endif

if (ierr.eq.1) then
    write(stderr,*) 'grid: coordinate definitions and grid_mode are incompatible'
    write(stderr,*) 'grid_mode = uniform requires: x, x-y, or x-y-z'
    write(stderr,*) 'grid_mode = dipping_plane requires: x-y'
    write(stderr,*) 'grid_mode = cross_section requires: x-z or x-y-z'
    write(stderr,*) 'isXDefined:  ',isXDefined
    write(stderr,*) 'isYDefined:  ',isYDefined
    write(stderr,*) 'isZDefined:  ',isZDefined
    write(stderr,*) 'grid_mode:    "',trim(grid_mode),'"'
    call usage('exiting')
endif


if (output_file.ne.'stdout') then
    close(lu_out)
endif

end

!--------------------------------------------------------------------------------------------------!

subroutine get_spacing(n,d,beg,end,mode,dim)
!----
! Compute the number of grid points or spacing/ratio between grid points
!----

use io, only: stderr
implicit none

! Arguments
integer :: n
double precision :: d, beg, end
character(len=*) :: mode, dim

! Local variables
integer, parameter :: nmax = 1000000


! Verify that mode is available
if (mode.eq.'linear') then
    ! all good
elseif (mode.eq.'exponential') then
    ! all good
else
    call usage('get_spacing: no spacing mode named "'//trim(mode)//'"')
endif


! Verify that one of number and increment is defined
if (n.gt.0.and.abs(d).lt.1.0d19) then
    write(stderr,*) 'get_spacing: both number and increment are specified for ',trim(dim)
    call usage('no_details')
elseif (n.le.0.and.abs(d).gt.1.0d19) then
    write(stderr,*) 'get_spacing: neither number nor increment are specified for ',trim(dim)
    call usage('use -n'//trim(dim)//' or -d'//trim(dim))
endif


! If only first value is set, then make sure we only print it
if (dabs(end).gt.1.0d98) then
    if (n.gt.1.or.dabs(d).lt.1.0d20) then
        write(stderr,*) 'get_spacing: only first ',trim(dim),' value defined, setting n',trim(dim),' to 1'
    endif
    n = 1
    d = -1.0d20
endif


! Compute spacing from number of points or vice versa
if (n.gt.0) then

    ! input n is number of grid points, including starting and ending points
    ! subtract one to compute even spacing
    n = n - 1

    if (n.eq.0) then
        if (mode.eq.'linear') then
            d = 0.0d0
        elseif (mode.eq.'exponential') then
            d = 1.0d0
        endif
    else
        if (mode.eq.'linear') then
            ! difference, d, between points is the same
            d = (end-beg)/dble(n)
        elseif (mode.eq.'exponential') then
            ! ratio, d, between points is the same
            d = exp(log(end/beg)/dble(n))
        endif
    endif

elseif (dabs(d).lt.1.0d19) then

    ! Input d is the spacing between points (difference or ratio)

    if (mode.eq.'linear') then

        ! Check that end value will be reached
        if (abs(d).lt.1.0d-20) then
            ! Spacing between points is zero
            write(stderr,*) 'get_spacing: ',trim(dim),' increment cannot be ==0 in linear mode'
            call usage('no_details')

        elseif ((end-beg)*d.lt.0) then
            ! Increment has wrong sign to get from beg to end
            write(stderr,1001) trim(dim),d,beg,end
            call usage('no_details')

        endif

        n = int((end-beg)/d)

    elseif (mode.eq.'exponential') then

        ! Check that end value will be reached
        if (d-1.0d0.lt.1.0d-20) then
            ! Ratio between points is one
            write(stderr,*) 'get_spacing: ',trim(dim),' increment cannot be ==1 in exponential mode'
            call usage('no_details')

        elseif (d.le.0.0d0) then
            ! Ratio between points is zero or negative
            write(stderr,*) 'get_spacing: ',trim(dim),' increment cannot be <=0 in exponential mode'
            call usage('no_details')

        elseif (dabs(beg).lt.1.0d-20) then
            ! Starting value is zero
            write(stderr,*) 'get_spacing: starting ',trim(dim),' value cannot be ==0 in ',&
                            'exponential mode'
            call usage('no_details')

        elseif (log(end/beg)*log(d).lt.0) then
            ! Ratio has wrong magnitude to get from beg to end
            write(stderr,1002) trim(dim),d,beg,end
            call usage('no_details')

        endif

        n = int(log(end/beg)/log(d))

    endif
endif


! Be reasonable about number of grid points
if (n.gt.nmax) then
    write(stderr,1011) nmax,trim(dim)
    call usage('no_details')
endif


1001 format(' get_spacing: ',A1,' inc=',1PE10.2,' has wrong sign for bounds [',E10.2,',',E10.2,']')
1002 format(' get_spacing: ',A1,' inc=',1PE10.2,' has wrong magnitude for bounds [',&
            E10.2,',',E10.2,']')
1011 format(' get_spacing: you requested more than ',I10,X,A4,' points')

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine get_value(v,i,d,beg,mode,dim)
!----
! Compute the value of the ith grid point
!----

use io, only: stderr
implicit none

! Arguments
integer :: i
double precision :: v, d, beg
character(len=*) :: mode, dim

if (mode.eq.'linear') then
    v = beg + d*dble(i-1)
elseif (mode.eq.'exponential') then
    v = beg * d**dble(i-1)
else
    call usage('get_value: no spacing mode named "'//trim(mode)//'"')
endif

if (abs(v).gt.1.0d99) then
    write(stderr,*) 'get_value: computed a giant value in dimension ',dim
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use io, only: stderr, verbosity, isNumeric

use grid, only: output_file, &
                grid_mode, &
                spacing_mode, &
                coord_type, &
                x_beg, &
                x_end, &
                x_inc, &
                nx, &
                isXDefined, &
                y_beg, &
                y_end, &
                y_inc, &
                ny, &
                isYDefined, &
                z_beg, &
                z_end, &
                z_inc, &
                nz, &
                isZDefined, &
                x0, &
                y0, &
                z0, &
                strike, &
                dip, &
                azimuth

implicit none

! Local variables
character(len=512) :: tag
integer :: i, j, narg, ios

! Initialize command line variables

! Output
output_file = 'stdout'

! Grid mode: uniform, dipping_plane, cross_section
grid_mode = 'uniform'

! Spacing mode: linear, exponential
spacing_mode = 'linear'

! Coordinate type: cartesian, geographic
coord_type = 'cartesian'

! Checks for output
isXDefined = .false.
isYDefined = .false.
isZDefined = .false.

! x parameters
x_beg = -1.0d99
x_end = -1.0d99
x_inc = -1.0d20
nx = -1

! y parameters
y_beg = -1.0d99
y_end = -1.0d99
y_inc = -1.0d20
ny = -1

! z parameters (uniform)
z_beg = -1.0d99
z_end = -1.0d99
z_inc = -1.0d20
nz = -1

! z parameters (dipping_plane, cross_section)
x0 = 0.0d0
y0 = 0.0d0
z0 = 0.0d0
strike = 0.0d0
dip = 0.0d0
azimuth = 0.0d0

narg = command_argument_count()
if (narg.eq.0) then
    call usage('')
endif

i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    !---- X dimension
    if (trim(tag).eq.'-x') then
        ! First argument (starting x value) is required
        i = i + 1
        call get_command_argument(i,tag)
        if (.not.isNumeric(tag)) call usage('grid: x_beg is not numeric')
        read (tag,*) x_beg
        isXDefined = .true.
        ! Second argument (ending x value) is optional
        i = i + 1
        if (i.gt.narg) then
            ! No more arguments
            nx = 1
            cycle
        else
            call get_command_argument(i,tag)
            ! if (.not.isNumeric(tag)) call usage('grid: x_end is not numeric')
            read (tag,*,iostat=ios) x_end
            if (ios.ne.0) then
                ! No second value provided
                nx = 1
                i = i - 1
            endif
        endif

    elseif (trim(tag).eq.'-nx') then
        x_inc = -1.0d20 ! Overrides x-increment definition
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) nx
    elseif (trim(tag).eq.'-dx') then
        nx = -1 ! Overrides nx definition
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) x_inc


    !---- Y dimension
    elseif (trim(tag).eq.'-y') then
        ! First argument (starting y value) is required
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) y_beg
        isYDefined = .true.
        ! Second argument (ending y value) is optional
        i = i + 1
        if (i.gt.narg) then
            ! No more arguments
            ny = 1
            cycle
        else
            call get_command_argument(i,tag)
            read (tag,*,iostat=ios) y_end
            if (ios.ne.0) then
                ! No second value provided
                ny = 1
                i = i - 1
            endif
        endif

    elseif (trim(tag).eq.'-ny') then
        y_inc = -1.0d20 ! Overrides y-increment definition
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) ny
    elseif (trim(tag).eq.'-dy') then
        ny = -1 ! Overrides ny definition
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) y_inc


    !---- Z dimension
    elseif (trim(tag).eq.'-z') then
        ! First argument (starting z value) is required
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) z_beg
        isZDefined = .true.
        if (grid_mode.ne.'cross_section') then
            grid_mode = 'uniform'
        endif
        ! Second argument (ending z value) is optional
        i = i + 1
        if (i.gt.narg) then
            ! No more arguments
            nz = 1
            cycle
        else
            call get_command_argument(i,tag)
            read (tag,*,iostat=ios) z_end
            if (ios.ne.0) then
                ! No second value provided
                nz = 1
                i = i - 1
            endif
        endif

    elseif (trim(tag).eq.'-nz') then
        z_inc = -1.0d20 ! Overrides z-increment definition
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) nz
    elseif (trim(tag).eq.'-dz') then
        nz = -1 ! Overrides nz definition
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) z_inc

    elseif (trim(tag).eq.'-xsec') then
        ! Vertical plane
        ! x -> distance along azimuth
        ! y -> distance perpendicular to plane
        ! z -> depth
        grid_mode = 'cross_section'
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) x0
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) y0
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) azimuth
    elseif (trim(tag).eq.'-xz') then
        write(stderr,*) 'grid: the -xz option is no longer needed, and is therefore deprecated'
        write(stderr,*) '-xsec prints all of the information'

    elseif (trim(tag).eq.'-dip') then
        ! Dipping plane
        grid_mode = 'dipping_plane'
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) x0
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) y0
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) z0
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) strike
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) dip

    elseif (trim(tag).eq.'-geo') then
        coord_type = 'geographic'

    elseif (trim(tag).eq.'-exp') then
        spacing_mode = 'exponential'

    elseif (tag(1:2).eq.'-R') then
        ! Input x, y ranges in GMT limits format
        ! Remove "-R" from string
        j = index(tag,'-R')
        tag(j:j+1) = '  '
        ! Remove three "/" from string
        j = index(tag,'/')
        tag(j:j) = ' '
        j = index(tag,'/')
        tag(j:j) = ' '
        j = index(tag,'/')
        tag(j:j) = ' '
        read(tag,*) x_beg, x_end, y_beg, y_end
        isXDefined = .true.
        isYDefined = .true.

    elseif (trim(tag).eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file)

    elseif (trim(tag).eq.'-v'.or.trim(tag).eq.'verbosity') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) verbosity

    else
        call usage('grid: no option '//tag)
    endif
    i = i + 1
enddo

return
end subroutine


!--------------------------------------------------------------------------------------------------!

subroutine usage(str)

use io, only: stderr

implicit none

! Arguments
character(len=*) :: str

! Local variables
integer :: i

if (str.ne.'') then
    i = index(str,'no_details')
    if (i.ne.0) then
        if (str.ne.'') then
            write(stderr,*) trim(str(1:i-1))
        endif
        call error_exit(1)
        write(stderr,*)
    endif
    write(stderr,*) trim(str)
    write(stderr,*)
endif

write(stderr,*) 'Usage: grid ...options...'
write(stderr,*)
write(stderr,*) '-o OFILE               Output to file (default prints to standard output)'
write(stderr,*) '-x XBEG [XEND]         Starting and (optional) ending x values'
write(stderr,*) '-nx NX                 Number of x points, including starting point'
write(stderr,*) '-dx DX                 Increment between x points'
write(stderr,*) '-y YBEG [YEND]         Starting and (optional) ending y values'
write(stderr,*) '-ny NY                 Number of y points, including starting point'
write(stderr,*) '-dy DY                 Increment between y points'
write(stderr,*) '-Rx1/x2/y1/y2          Define x, y limits with GMT limits format'
write(stderr,*) '-z ZBEG [ZEND]         Starting and (optional) ending z values'
write(stderr,*) '-nz NZ                 Number of z points, including starting point'
write(stderr,*) '-dz DZ                 Increment between z points'
write(stderr,*) '-dip X0 Y0 Z0 STR DIP  Z-values lie on plane containing (X0,Y0,Z0) with STR/DIP'
write(stderr,*) '-xsec X0 Y0 AZ         Vertical cross section through (X0,Y0) with strike AZ ',&
                                        '(x=az; y=az+90; z=up)'
write(stderr,*) '-geo                   Treat -dip and -xsec coordinates as lon/lat'
write(stderr,*) '-exp                   Make grid increment exponential'
write(stderr,*)
write(stderr,*) 'See man page for details'
write(stderr,*)

call error_exit(1)
end subroutine
