!--------------------------------------------------------------------------------------------------!
! Module: ffm
!
! Routines to parse finite fault models in different standard formats and variables to contain the
! fault data.
!
! Standard finite fault model formats implemented in this module are:
! - U.S. Geological Survey .param format (previously known as static_out format)
! - SRCMOD database FSP format
! - GMT psmeca format
! - GMT psmeca format, replacing magnitude with fault slip, width, and length
!
!--------------------------------------------------------------------------------------------------!

module ffm

type ffm_data
    integer :: nflt
    integer :: nseg
    double precision, allocatable :: subflt(:,:)
    integer, allocatable :: seg(:)
    double precision :: hylo
    double precision :: hyla
    double precision :: hydp
    double precision :: mag
    double precision :: mom
end type

public :: ffm_data

public :: read_usgs_param
public :: read_srcmod_fsp
public :: get_ffm_outline

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

subroutine read_usgs_param(usgs_param_file,usgs_param_data,ierr)
!----
! Read a finite fault model in the standard .param format published by the U.S. Geological Survey.
!----

use io, only: stderr

implicit none

! Arguments
character(len=*) :: usgs_param_file
type(ffm_data) :: usgs_param_data
integer :: ierr

! Local variables
logical :: ex
integer :: i, j, n, nx, ny, ios, iseg
character(len=1) :: dum
character(len=16) :: dx_char, dy_char
double precision :: lon, lat, dep, strike, dip, rake, slip, dx, dy, trup

! As we begin our journey, everything seems okay...
ierr = 0

! Open the FFM file
inquire(exist=ex,file=usgs_param_file)
if (.not.ex) then
    write(stderr,*) 'read_usgs_param: no USGS FFM file found named ',trim(usgs_param_file)
    ierr = 1
    return
endif

open(unit=34,file=usgs_param_file,status='old',err=1001,iostat=ios)
1001 if (ios.ne.0) then
    write(stderr,*) 'read_usgs_param: error opening USGS FFM file ',trim(usgs_param_file)
    ierr = 1
    return
endif

! Header indicates the number of fault segments
read (34,*) dum,dum,dum,dum,usgs_param_data%nseg

! Count the total number of sub-faults
usgs_param_data%nflt = 0
do iseg = 1,usgs_param_data%nseg
    read(34,*) dum,dum,dum,dum,nx,dum,dum,dum,ny
    usgs_param_data%nflt = usgs_param_data%nflt + nx*ny
    do j = 1,8+nx*ny
        read(34,*)
    enddo
enddo
rewind(34)

! Allocate memory to the sub-fault array
if (allocated(usgs_param_data%subflt)) then
    deallocate(usgs_param_data%subflt)
endif
allocate(usgs_param_data%subflt(usgs_param_data%nflt,10))

! Allocate memory to the segment array
if (allocated(usgs_param_data%seg)) then
    deallocate(usgs_param_data%seg)
endif
allocate(usgs_param_data%seg(usgs_param_data%nflt))

! Now fill the sub-fault array
n = 0

! First line with number of segments
read (34,*) dum,dum,dum,dum,usgs_param_data%nseg

! Read details of each segment
do iseg = 1,usgs_param_data%nseg

    ! Sub-fault dimensions
    read (34,*) dum,dum,dum,dum,nx,dum,dx_char,dum,ny,dum,dy_char
    j = index(dx_char,'km')
    dx_char(j:j+1) = ''
    read(dx_char,*) dx
    j = index(dy_char,'km')
    dy_char(j:j+1) = ''
    read(dy_char,*) dy
    usgs_param_data%subflt(n+1:n+nx*ny,8) = dy*1.0d3
    usgs_param_data%subflt(n+1:n+nx*ny,9) = dx*1.0d3

    ! Hypocenter
    read (34,*) (dum,j=1,10),usgs_param_data%hylo,dum,usgs_param_data%hyla

    ! Not sure what this block is actually...
    do j = 1,7
        read(34,*) dum
    enddo

    ! Sub-fault locations and slip
    do i = n+1,n+nx*ny
        read (34,*) lat,lon,dep,slip,rake,strike,dip,trup
        usgs_param_data%subflt(i,1) = lon
        usgs_param_data%subflt(i,2) = lat
        usgs_param_data%subflt(i,3) = dep*1.0d3
        usgs_param_data%subflt(i,4) = strike
        usgs_param_data%subflt(i,5) = dip
        usgs_param_data%subflt(i,6) = rake
        usgs_param_data%subflt(i,7) = slip*1.0d-2
        usgs_param_data%subflt(i,10) = trup
        usgs_param_data%seg(i) = iseg
     enddo

     n = n + nx*ny
enddo

close(34)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine read_srcmod_fsp(srcmod_fsp_file,srcmod_fsp_data,ierr)
!----
! Read a finite fault model in the standard FSP format published in the SRCMOD database.
!----

use io, only: stderr

implicit none

! Arguments
character(len=*) :: srcmod_fsp_file
type(ffm_data) :: srcmod_fsp_data
integer :: ierr

! Local variables
logical :: ex
integer :: ios, iline, i, j, n, ncols
character(len=256) :: line
integer, parameter :: NFMAX = 25
character(len=32) :: linef(NFMAX)
character(len=32) :: linef_sf_cols(NFMAX)
double precision :: strike, dip, rake, wid, len

! Let us begin
ierr = 0
ncols = 0

! Open the FSP file
inquire(exist=ex,file=srcmod_fsp_file)
if (.not.ex) then
    write(stderr,*) 'read_srcmod_fsp: no FSP file found named ',trim(srcmod_fsp_file)
    ierr = 1
    return
endif

open(unit=35,file=srcmod_fsp_file,status='old',err=1001,iostat=ios)
1001 if (ios.ne.0) then
    write(stderr,*) 'read_srcmod_fsp: error opening FSP file ',trim(srcmod_fsp_file)
    ierr = 1
    return
endif


! Count the number of sub-faults in the FSP file
srcmod_fsp_data%nflt = 0
do
    linef = ''
    read(35,*,iostat=ios) linef(1)
    if (ios.ne.0) then
        exit
    endif
    if (index(linef(1),'%').eq.1.or.linef(1).eq.'') then
        ! Header or empty line: do not count
    else
        srcmod_fsp_data%nflt = srcmod_fsp_data%nflt + 1
    endif
enddo
rewind(35)


! Allocate memory to the sub-fault array
if (allocated(srcmod_fsp_data%subflt)) then
    deallocate(srcmod_fsp_data%subflt)
endif
allocate(srcmod_fsp_data%subflt(srcmod_fsp_data%nflt,10))
srcmod_fsp_data%subflt = 0.0d0

! Allocate memory to the segment array
if (allocated(srcmod_fsp_data%seg)) then
    deallocate(srcmod_fsp_data%seg)
endif
allocate(srcmod_fsp_data%seg(srcmod_fsp_data%nflt))
srcmod_fsp_data%nseg = 1
srcmod_fsp_data%seg = 1


! Read the FSP file
iline = 1
n = 0
do
    ! Read entire line, checking for end of file and errors
    read(35,'(A)',iostat=ios,end=1002,err=1003) line
    ! write(stderr,*) 'line: ',trim(line)

    ! Ignore empty lines
    if (line.eq.'') then
        cycle
    endif

    ! Read elements of line
    linef = ''
    read(line,*,end=101) (linef(i),i=1,NFMAX)
    101 continue

    ! Is this line a header?
    if (linef(1)(1:1).eq.'%') then

        ! Found a header line: parse it
        ! write(stderr,*) 'read_fsp: line ',iline,' is a header'

        i = 2
        do while (i.le.NFMAX-1)
            if (linef(2).eq.'Loc'.and.linef(i).eq.'LAT'.and.linef(i+1).eq.'=') then
                i = i + 2
                ! write(stderr,*) 'hyla ',linef(i)
                read(linef(i),*) srcmod_fsp_data%hyla

            elseif (linef(2).eq.'Loc'.and.linef(i).eq.'LON'.and.linef(i+1).eq.'=') then
                i = i + 2
                ! write(stderr,*) 'hylo ',linef(i)
                read(linef(i),*) srcmod_fsp_data%hylo

            elseif (linef(2).eq.'Loc'.and.linef(i).eq.'DEP'.and.linef(i+1).eq.'=') then
                i = i + 2
                ! write(stderr,*) 'hydp ',linef(i)
                read(linef(i),*) srcmod_fsp_data%hydp

            elseif (linef(i).eq.'Mw'.and.linef(i+1).eq.'=') then
                i = i + 2
                ! write(stderr,*) 'mag ',linef(i)
                read(linef(i),*) srcmod_fsp_data%mag

            elseif (linef(i).eq.'Mo'.and.linef(i+1).eq.'=') then
                i = i + 2
                ! write(stderr,*) 'mom ',linef(i)
                read(linef(i),*) srcmod_fsp_data%mom

            elseif ((linef(i).eq.'STRK'.or.linef(i).eq.'STRIKE').and.linef(i+1).eq.'=') then
                i = i + 2
                ! write(stderr,*) 'str ',linef(i)
                read(linef(i),*) strike
                srcmod_fsp_data%subflt(n+1:srcmod_fsp_data%nflt,4) = strike

            elseif (linef(i).eq.'DIP'.and.linef(i+1).eq.'=') then
                i = i + 2
                ! write(stderr,*) 'dip ',linef(i)
                read(linef(i),*) dip
                srcmod_fsp_data%subflt(n+1:srcmod_fsp_data%nflt,5) = dip

            elseif (linef(i).eq.'RAKE'.and.linef(i+1).eq.'=') then
                i = i + 2
                ! write(stderr,*) 'rak ',linef(i)
                read(linef(i),*) rake
                srcmod_fsp_data%subflt(n+1:srcmod_fsp_data%nflt,6) = rake

            elseif (linef(i).eq.'Dx'.and.linef(i+1).eq.'=') then
                i = i + 2
                ! write(stderr,*) 'len ',linef(i)
                read(linef(i),*) len
                srcmod_fsp_data%subflt(n+1:srcmod_fsp_data%nflt,9) = len*1.0d3

            elseif (linef(i).eq.'Dz'.and.linef(i+1).eq.'=') then
                i = i + 2
                ! write(stderr,*) 'wid ',linef(i)
                read(linef(i),*) wid
                srcmod_fsp_data%subflt(n+1:srcmod_fsp_data%nflt,8) = wid*1.0d3

            elseif (linef(i).eq.'SEGMENT'.and.linef(i+1).eq.'#') then
                i = i + 2
                ! write(stderr,*) 'wid ',linef(i)
                j = index(linef(i),':')
                linef(i)(j:j) = ' '
                read(linef(i),*) srcmod_fsp_data%nseg
                srcmod_fsp_data%seg(n+1:srcmod_fsp_data%nflt) = srcmod_fsp_data%nseg

            elseif (linef(i).eq.'LAT'.and.linef(i+1).eq.'LON') then
                ! This line describes the available sub-fault data: save the list
                ! write(0,*) iline,trim(line)
                linef_sf_cols = linef
                do j = i+2,NFMAX
                    if (linef(j).eq.'') then
                        ncols = j-1
                        exit
                    endif
                    ncols = j
                enddo
            endif
            i = i + 1
        enddo

    else

       ! This is a sub-fault line; parse it with previously read order
       ! write(0,*) 'iline',iline,' ncols',ncols

        n = n + 1
        do i = 2,ncols
            if (linef_sf_cols(i).eq.'LON') then
                read(linef(i-1),*) srcmod_fsp_data%subflt(n,1)

            elseif (linef_sf_cols(i).eq.'LAT') then
                read(linef(i-1),*) srcmod_fsp_data%subflt(n,2)

            elseif (linef_sf_cols(i).eq.'Z') then
                read(linef(i-1),*) srcmod_fsp_data%subflt(n,3)
                srcmod_fsp_data%subflt(n,3) = srcmod_fsp_data%subflt(n,3)*1.0d3

            elseif (linef_sf_cols(i).eq.'RAKE') then
                read(linef(i-1),*) srcmod_fsp_data%subflt(n,6)

            elseif (linef_sf_cols(i).eq.'SLIP') then
                read(linef(i-1),*) srcmod_fsp_data%subflt(n,7)

            elseif (linef_sf_cols(i).eq.'TRUP') then
                read(linef(i-1),*) srcmod_fsp_data%subflt(n,10)

            endif
        enddo
        ! write(stderr,*) srcmod_fsp_data%subflt(n,:)
    endif

1002 if (ios.ne.0) then
        exit ! end of file
    endif
1003 if (ios.ne.0) then
        write(stderr,*) 'read_srcmod_fsp: error reading line ',iline,' from FSP file ', &
                        trim(srcmod_fsp_file)
        ierr = 1
        return
    endif

    iline = iline + 1
enddo

close(35)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine get_ffm_outline(input_ffm,outline)
!----
! Determine the (quadrilateral) outline for a finite fault model when the input sub-faults
! are rectangular patches.
!----

use trig, only: d2r
use geom, only: distaz2lola

implicit none

! Arguments
type(ffm_data) :: input_ffm
double precision :: outline(input_ffm%nseg,4,2)

! Local variables
integer :: i, j, iseg, iflt, ioutline(input_ffm%nseg,4)
logical :: isSegInit(input_ffm%nseg)
double precision :: evlo, evla, str, dip, wid, len, lo, la, corners(4,2)

! Initialize
outline = 0.0d0
isSegInit = .false.

! Determine the outline of each segment
do i = 1,input_ffm%nflt

    ! Define more useful variables
    iseg = input_ffm%seg(i)
    evlo = input_ffm%subflt(i,1)
    evla = input_ffm%subflt(i,2)

    ! Set initial outline value to first sub-fault coordinate
    if (.not.isSegInit(iseg)) then
        outline(iseg,1:2,1) = evlo
        outline(iseg,3:4,2) = evla
        ioutline(iseg,:) = i
        isSegInit(iseg) = .true.
    endif

    ! Minimum longitude of segment
    if (evlo.lt.outline(iseg,1,1)) then
        outline(iseg,1,1) = evlo
        ioutline(iseg,1) = i
    endif

    ! Maximum longitude of segment
    if (evlo.gt.outline(iseg,2,1)) then
        outline(iseg,2,1) = evlo
        ioutline(iseg,2) = i
    endif

    ! Minimum latitude of segment
    if (evla.lt.outline(iseg,3,2)) then
        outline(iseg,3,2) = evla
        ioutline(iseg,3) = i
    endif

    ! Maximum latitude of segment
    if (evla.gt.outline(iseg,4,2)) then
        outline(iseg,4,2) = evla
        ioutline(iseg,4) = i
    endif
enddo

! Move max and min values from fault center to fault corner
do iseg = 1,input_ffm%nseg
    do i = 1,4
        iflt = ioutline(iseg,i)
        if (iflt.lt.1.or.iflt.gt.input_ffm%nflt) then
            cycle
        endif

        evlo = input_ffm%subflt(iflt,1)
        evla = input_ffm%subflt(iflt,2)
        str = input_ffm%subflt(iflt,4)
        dip = input_ffm%subflt(iflt,5)
        wid = input_ffm%subflt(iflt,8)
        len = input_ffm%subflt(iflt,9)
        call distaz2lola(evlo,evla,len/2.0d0/6371.0d3,str,lo,la)
        call distaz2lola(lo,la,wid/2.0d0*cos(dip*d2r)/6371.0d3,str-90.0d0,corners(1,1),corners(1,2))
        call distaz2lola(lo,la,wid/2.0d0*cos(dip*d2r)/6371.0d3,str+90.0d0,corners(2,1),corners(2,2))
        call distaz2lola(evlo,evla,len/2.0d0/6371.0d3,str+180.0d0,lo,la)
        call distaz2lola(lo,la,wid/2.0d0*cos(dip*d2r)/6371.0d3,str-90.0d0,corners(3,1),corners(3,2))
        call distaz2lola(lo,la,wid/2.0d0*cos(dip*d2r)/6371.0d3,str+90.0d0,corners(4,1),corners(4,2))

        if (i.eq.1) then
            outline(iseg,i,1) = minval(corners(:,1))
            do j = 1,4
                if (abs(corners(j,1)-outline(iseg,i,1)).lt.1.0d-10) then
                    outline(iseg,i,2) = corners(j,2)
                endif
            enddo
        elseif (i.eq.2) then
            outline(iseg,i,1) = maxval(corners(:,1))
            do j = 1,4
                if (abs(corners(j,1)-outline(iseg,i,1)).lt.1.0d-10) then
                    outline(iseg,i,2) = corners(j,2)
                endif
            enddo
        elseif (i.eq.3) then
            outline(iseg,i,2) = minval(corners(:,2))
            do j = 1,4
                if (abs(corners(j,2)-outline(iseg,i,2)).lt.1.0d-10) then
                    outline(iseg,i,1) = corners(j,1)
                endif
            enddo
        elseif (i.eq.4) then
            outline(iseg,i,2) = maxval(corners(:,2))
            do j = 1,4
                if (abs(corners(j,2)-outline(iseg,i,2)).lt.1.0d-10) then
                    outline(iseg,i,1) = corners(j,1)
                endif
            enddo
        endif

    enddo
enddo

! do iseg = 1,input_ffm%nseg
!     write(0,*) iseg
!     do j = 1,4
!         do i = 1,2
!             write(0,*) outline(iseg,j,i)
!         enddo
!     enddo
! enddo

return
end subroutine


!--------------------------------------------------------------------------------------------------!

subroutine read_mag(mag_file,mag_data,ierr)

!----
! Read shear dislocations from fault file in GMT psmeca -Sa format:
!     evlo evla evdp(km) str dip rak mag
!----

use io, only: stderr, fileExists, line_count
implicit none

! Arguments
character(len=*) :: mag_file
type(ffm_data) :: mag_data
integer :: ierr

! Local variables
integer :: iflt, j


! In the beginning, it was good
ierr = 0

! Check that the file exists
if (.not.fileExists(mag_file)) then
    write(stderr,*) 'read_mag: no fault file found named ',trim(mag_file)
    ierr = 1
    return
endif

! Allocate memory to sub-fault array
mag_data%nflt = line_count(mag_file)
if (allocated(mag_data%subflt)) then
    deallocate(mag_data%subflt)
endif
allocate(mag_data%subflt(mag_data%nflt,9))

! Read the file
open(unit=31,file=mag_file,status='old')
do iflt = 1,mag_data%nflt
    read(31,*) (mag_data%subflt(iflt,j),j=1,7)
    mag_data%subflt(iflt,3) = mag_data%subflt(iflt,3)*1.0d3 ! Depth km->m
enddo
close(31)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine read_flt(flt_file,flt_data,ierr)
!----
! Read shear dislocations from fault file in GMT psmeca -Sa format, except with slip wid len:
!     evlo evla evdp(km) str dip rak slip(m) wid(km) len(km)
!----

use io, only: stderr, fileExists, line_count
implicit none

! Arguments
character(len=*) :: flt_file
type(ffm_data) :: flt_data
integer :: ierr

! Local variables
integer :: iflt, j


! Yayyy, things are okay! Or are they?
ierr = 0

! Check that the file exists
if (.not.fileExists(flt_file)) then
    write(stderr,*) 'read_flt: no fault file found named ',trim(flt_file)
    ierr = 1
    return
endif

! Allocate memory to sub-fault array
flt_data%nflt = line_count(flt_file)
if (allocated(flt_data%subflt)) then
    deallocate(flt_data%subflt)
endif
allocate(flt_data%subflt(flt_data%nflt,9))

! Read the file
open(unit=32,file=flt_file,status='old')
do iflt = 1,flt_data%nflt
    read(32,*) (flt_data%subflt(iflt,j),j=1,9)
    flt_data%subflt(iflt,3) = flt_data%subflt(iflt,3)*1.0d3 ! Depth km->m
    flt_data%subflt(iflt,8) = flt_data%subflt(iflt,8)*1.0d3 ! Width km->m
    flt_data%subflt(iflt,9) = flt_data%subflt(iflt,9)*1.0d3 ! Length km->m
enddo
close(32)

return
end subroutine



end module

!==================================================================================================!
