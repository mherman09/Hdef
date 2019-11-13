module ff2gmt

character(len=512) :: input_file
character(len=16) :: ffm_type

character(len=512) :: slip_file
character(len=512) :: time_file
character(len=512) :: depth_file
character(len=512) :: area_file

character(len=512) :: clip_file
character(len=16) :: clip_type
character(len=512) :: epi_file

logical :: isOutputDefined

end module

!==================================================================================================!

program main
!----
! Read a finite fault model and convert it to a format for use with GMT
!----

use io, only: stderr, fileExists
use trig, only: d2r
use ffm, only: ffm_data, read_srcmod_fsp, read_usgs_param, get_ffm_outline

use ff2gmt, only: input_file, &
                  ffm_type, &
                  slip_file, &
                  time_file, &
                  depth_file, &
                  area_file, &
                  clip_file, &
                  clip_type, &
                  epi_file, &
                  isOutputDefined

implicit none

! Local variables
integer :: i, ierr
type(ffm_data) :: ff
double precision, allocatable :: clip(:,:,:), evlo, evla, evdp, str, dip, rak, slip, wid, len, time


! Parse command line
call gcmdln()


! Check that input/output are defined
if (input_file.eq.'') then
    call usage('ff2gmt: input finite fault model required')
else
    if (.not.fileExists(input_file)) then
        call usage('ff2gmt: no finite fault model file found named '//trim(input_file))
    endif
endif

if (.not.isOutputDefined) then
    call usage('ff2gmt: no output defined')
endif


! Read the finite fault model
if (ffm_type.eq.'usgs_param') then
    call read_usgs_param(input_file,ff,ierr)
elseif (ffm_type.eq.'srcmod_fsp') then
    call read_srcmod_fsp(input_file,ff,ierr)
else
    call usage('ff2gmt: no ffm_type named '//trim(ffm_type))
endif

if (ierr.ne.0) then
    call usage('ff2gmt: error reading finite fault model')
endif


! Open files for specified outputs
if (slip_file.ne.'') then
    open(unit=21,file=slip_file,status='unknown')
endif
if (time_file.ne.'') then
    if (maxval(ff%subflt(:,10)).lt.1.0d0) then
        write(stderr,*) 'ff2gmt: FFM has no timing information'
        time_file = ''
    else
        open(unit=22,file=time_file,status='unknown')
    endif
endif
if (depth_file.ne.'') then
    open(unit=23,file=depth_file,status='unknown')
endif
if (area_file.ne.'') then
    open(unit=24,file=area_file,status='unknown')
endif
if (epi_file.ne.'') then
    open(unit=31,file=epi_file,status='unknown')
endif
if (clip_file.ne.'') then
    open(unit=41,file=clip_file,status='unknown')
endif


! Write FFM data to files

! Sub-fault data goes in psxy -SJ -C<cpt> format: lon lat <value> strike len wid*cos(dip)
do i = 1,ff%nflt
    evlo = ff%subflt(i,1)
    evla = ff%subflt(i,2)
    evdp = ff%subflt(i,3)/1e3
    str = ff%subflt(i,4)
    dip = ff%subflt(i,5)
    rak = ff%subflt(i,6)
    slip = ff%subflt(i,7)
    wid = ff%subflt(i,8)/1d3
    len = ff%subflt(i,9)/1d3
    time = ff%subflt(i,10)
    if (slip_file.ne.'') then
        write(21,*) evlo,evla,slip,str,len,wid*cos(dip*d2r)
    endif
    if (time_file.ne.'') then
        write(22,*) evlo,evla,time,str,len,wid*cos(dip*d2r)
    endif
    if (depth_file.ne.'') then
        write(23,*) evlo,evla,evdp,str,len,wid*cos(dip*d2r)
    endif
    if (area_file.ne.'') then
        write(24,*) len*wid,len,wid
    endif
enddo

! Epicenter data plots as a single point
if (epi_file.ne.'') then
    write(31,*) ff%hylo, ff%hyla
endif

! Clipping mask is the fault outline
if (clip_file.ne.'') then
    if (clip_type.eq.'seg') then
        ! Everything is initially set to treat each segment separately
    elseif (clip_type.eq.'all') then
        ! Re-define FFM sub-faults to all be in first segment
        ff%seg = 1
    else
        write(stderr,*) 'ff2gmt: no clip_type named ',trim(clip_type)
    endif

    ! Extract the outline to clip array
    allocate(clip(ff%nseg,4,2))
    call get_ffm_outline(ff,clip)

    ! Redefine number of fault segments
    if (clip_type.eq.'all') then
        ff%nseg = 1
    endif

    ! Write fault outline to GMT multi-segment file
    do i = 1,ff%nseg
        write(41,'(A)') '>'
        write(41,*) clip(i,1,1),clip(i,1,2)
        write(41,*) clip(i,3,1),clip(i,3,2)
        write(41,*) clip(i,2,1),clip(i,2,2)
        write(41,*) clip(i,4,1),clip(i,4,2)
        write(41,*) clip(i,1,1),clip(i,1,2)
    enddo
    deallocate(clip)
endif


! All done. Shut it down.
close(21)
close(22)
close(23)
close(24)
close(31)
close(41)

end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use ff2gmt, only: input_file, &
                  ffm_type, &
                  slip_file, &
                  time_file, &
                  depth_file, &
                  area_file, &
                  clip_file, &
                  clip_type, &
                  epi_file, &
                  isOutputDefined

implicit none

! Local variables
integer :: i, narg
character(len=512) :: tag

! Initialize control parameters
input_file = ''
ffm_type = ''
slip_file = ''
time_file = ''
depth_file = ''
area_file = ''
clip_file = ''
epi_file = ''
isOutputDefined = .false.

! Number of arguments
narg = command_argument_count()
if (narg.eq.0) then
  call usage('')
endif

! Read arguments
i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    if (trim(tag).eq.'-ffm') then
        ffm_type = 'usgs_param'
        i = i + 1
        call get_command_argument(i,input_file)
    elseif (trim(tag).eq.'-fsp') then
        ffm_type = 'srcmod_fsp'
        i = i + 1
        call get_command_argument(i,input_file)

    elseif (trim(tag).eq.'-slip') then
        isOutputDefined = .true.
        i = i + 1
        call get_command_argument(i,slip_file)
    elseif (trim(tag).eq.'-trup') then
        isOutputDefined = .true.
        i = i + 1
        call get_command_argument(i,time_file)
    elseif (trim(tag).eq.'-dep') then
        isOutputDefined = .true.
        i = i + 1
        call get_command_argument(i,depth_file)
    elseif (trim(tag).eq.'-area') then
        isOutputDefined = .true.
        i = i + 1
        call get_command_argument(i,area_file)

    elseif (trim(tag).eq.'-clipseg') then
        isOutputDefined = .true.
        clip_type = 'seg'
        i = i + 1
        call get_command_argument(i,clip_file)
    elseif (trim(tag).eq.'-clip') then
        isOutputDefined = .true.
        clip_type = 'all'
        i = i + 1
        call get_command_argument(i,clip_file)
    elseif (trim(tag).eq.'-epi') then
        isOutputDefined = .true.
        i = i + 1
        call get_command_argument(i,epi_file)

    else
        call usage('ff2gmt: no command line option '//trim(tag))
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

write(stderr,*) 'Usage: ff2gmt ...options...'
write(stderr,*)
write(stderr,*) 'Input'
write(stderr,*) '-ffm FFMFILE       Finite fault model in USGS .param format'
write(stderr,*) '-fsp FSPFILE       Finite fault model in SRCMOD FSP format'
write(stderr,*)
write(stderr,*) 'Output (for use with "gmt psxy -SJ -C<cpt>": lon lat <value> str len wid)'
write(stderr,*) '-slip SLIPFILE     Slip magnitude'
write(stderr,*) '-trup TIMEFILE     Rupture time'
write(stderr,*) '-dep DEPFILE       Sub-fault depth'
write(stderr,*)
write(stderr,*) 'Other options'
write(stderr,*) '-clip CLIPFILE     Write outline of FFM to file'
write(stderr,*) '-clipseg CLIPFILE  Write outline of each segment of FFM to file'
write(stderr,*) '-epi EPIFILE       Write epicenter to file'
write(stderr,*) '-area AREAFILE     Sub-fault area, length, and width'
write(stderr,*)

stop
end subroutine
