module eqempirical

character(len=512) :: relation
double precision :: magnitude

end module

!==================================================================================================!

program main
!----
! Print fault characteristics (dimensions, slip) from published empirical relations
!----

use io, only: stdout
use eq, only: empirical
use eqempirical, only: relation, magnitude

implicit none

! Local variables
double precision :: len, wid


! Parse command line inputs
call gcmdln()
if (relation.eq.'') then
    call usage('eqempirical: empirical relation not defined')
endif
if (magnitude.lt.-10.0d0 .and. relation.ne.'list') then
    call usage('eqempirical: magnitude not defined')
endif


if (relation.eq.'WC') then
    call empirical(magnitude,wid,len,'WC:print','')
elseif (relation.eq.'MB') then
    call empirical(magnitude,wid,len,'MB:print','')
elseif (relation.eq.'B') then
    call empirical(magnitude,wid,len,'B:print','')
elseif (relation.eq.'YM') then
    call empirical(magnitude,wid,len,'YM:print','')
elseif (relation.eq.'AH') then
    call empirical(magnitude,wid,len,'AH:print','')
elseif (relation.eq.'list') then
    write(stdout,*) 'Empirical Relations:'
    write(stdout,*) '  WC: Wells and Coppersmith, BSSA (1994)'
    write(stdout,*) '  MB: Mai and Beroza, BSSA (2000)'
    write(stdout,*) '  B: Blaser et al., BSSA (2010)'
    write(stdout,*) '  YM: Yen and Ma, BSSA (2011)'
    write(stdout,*) '  AH: Allen and Hayes, BSSA (2017)'
else
    call usage('eqempirical: no empirical relationship named '//trim(relation))
endif

end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use eqempirical, only: relation, &
                       magnitude
implicit none

! Local variables
integer :: i, ios, narg
character(len=512) :: tag, arg


! Initialize control parameters
ios = 0
relation = ''
magnitude = -1d10


! Number of arguments
narg = command_argument_count()
if (narg.eq.0) then
  call usage('')
endif

! Read arguments
i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    if (tag.eq.'-model') then
        i = i + 1
        call get_command_argument(i,arg,status=ios)
        read (arg,*,err=9032,iostat=ios) relation
    elseif (tag.eq.'-mag') then
        i = i + 1
        call get_command_argument(i,arg,status=ios)
        read (arg,*,err=9032,iostat=ios) magnitude
    elseif (tag.eq.'-list') then
        relation = 'list'
    else
        call usage('')
    endif

    9032 if (ios.ne.0) then
        call usage('eqempirical: error parsing "'//trim(tag)//'" flag arguments')
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

write(stderr,*) 'Usage: eqempirical -model MODEL -mag MAGNITUDE'
write(stderr,*)
write(stderr,*) '-model MODEL      Empirical relation name'
write(stderr,*) '-mag MAGNITUDE    Earthquake magnitude'
write(stderr,*) '-list             List available models'
write(stderr,*)

call error_exit(1)
end subroutine
