module io

integer :: verbosity
character(len=512) :: gcmtFile
character(len=8) :: geoMode
double precision :: lonMin, lonMax, latMin, latMax, lon0, lat0, distMax
double precision :: magMin, magMax
double precision :: jdB, jdE

contains

    subroutine gcmdln()
    implicit none
    integer :: i, narg
    character(len=512) :: tag
    character(len=20) :: timeBeg, timeEnd
    integer :: tB(6), tE(6)
    double precision :: dyB, dyE

    ! Command line variables
    gcmtFile = '/Users/mwh5316/Downloads/jan76_dec17.ndk'
    verbosity = 0
    geoMode = 'rect' ! rect, circle
    lonmin = -360
    lonmax = 360
    latmin = -90
    latmax = 90
    lon0 = 0
    lat0 = 0
    distMax = 111.19*180
    magMin = 0.0
    magMax = 10.0
    timeBeg = '1976-01-01T00:00:00'
    timeEnd = '2076-01-01T00:00:00'

    narg = command_argument_count()
    if (narg.eq.0) then
        call usage('')
    endif
    i = 1
    do while (i.le.narg)
        call get_command_argument(i,tag)
        if (trim(tag).eq.'-f') then
            i = i + 1
            call get_command_argument(i,gcmtFile)
        elseif (trim(tag).eq.'-v') then
            verbosity = 1
        elseif (trim(tag).eq.'-rect') then
            geoMode = 'rect'
            i = i + 1
            call get_command_argument(i,tag)
            read(tag,*) lonMin
            i = i + 1
            call get_command_argument(i,tag)
            read(tag,*) lonMax
            i = i + 1
            call get_command_argument(i,tag)
            read(tag,*) latMin
            i = i + 1
            call get_command_argument(i,tag)
            read(tag,*) latMax
        elseif (trim(tag).eq.'-circle') then
            geoMode = 'circle'
            i = i + 1
            call get_command_argument(i,tag)
            read(tag,*) lon0
            i = i + 1
            call get_command_argument(i,tag)
            read(tag,*) lat0
            i = i + 1
            call get_command_argument(i,tag)
            read(tag,*) distMax
        elseif (trim(tag).eq.'-mag') then
            i = i + 1
            call get_command_argument(i,tag)
            read(tag,*) magMin
            i = i + 1
            call get_command_argument(i,tag)
            read(tag,*) magMax
        elseif (trim(tag).eq.'-date') then
            i = i + 1
            call get_command_argument(i,timeBeg)
            i = i + 1
            call get_command_argument(i,timeEnd)
        else
            call usage('!! Error: No option '//trim(tag))
        endif
        i = i + 1
    enddo

    ! timeBeg: YYYY-MM-DDTHH:MM:SS
    read(timeBeg(1:4),'(I4)') tB(1)
    read(timeBeg(6:7),'(I2)') tB(2)
    read(timeBeg(9:10),'(I2)') tB(3)
    read(timeBeg(12:13),'(I2)') tB(4)
    read(timeBeg(15:16),'(I2)') tB(5)
    read(timeBeg(18:19),'(I2)') tB(6)
    dyB = dble(tB(3)) + dble(tB(4))/2.4d1 + dble(tB(5))/1.44d3 + dble(tB(6))/8.64d4
    call date2jd(jdB,tB(1),tB(2),dyB)
    ! timeEnd: YYYY-MM-DDTHH:MM:SS
    read(timeEnd(1:4),'(I4)') tE(1)
    read(timeEnd(6:7),'(I2)') tE(2)
    read(timeEnd(9:10),'(I2)') tE(3)
    read(timeEnd(12:13),'(I2)') tE(4)
    read(timeEnd(15:16),'(I2)') tE(5)
    read(timeEnd(18:19),'(I2)') tE(6)
    dyE = dble(tE(3)) + dble(tE(4))/2.4d1 + dble(tE(5))/1.44d3 + dble(tE(6))/8.64d4
    call date2jd(jdE,tE(1),tE(2),dyE)

    !write(0,*) 'mag ',magMin,magMax
    !write(0,*) 'date ',timeBeg,timeEnd,dyE,dyB,jdB,jdE

    return
    end subroutine gcmdln

end module io

!--------------------------------------------------------------------------------------------------!

module event_data

! Input lines
character(len=80) :: line(5)

! Hypocenter information
character(len=4) :: hypoCatalog
integer :: origin(7)
double precision :: hypoLat, hypoLon, hypoDep, mb, ms
character(len=32) :: eventLabel

! CMT information
character(len=16) :: gcmtID
integer :: nstaB,ncompB,nstaS,ncompS,nstaM,ncompM
double precision :: minPeriodB,minPeriodS,minPeriodM
integer :: gcmtType
character(len=8) :: stfType
double precision :: stfHalfDuration

! CMT origin and location parameters from inversion
double precision :: centroidTime, centroidTimeError
double precision :: centroidLat, centroidLatError
double precision :: centroidLon, centroidLonError
double precision :: centroidDep, centroidDepError
character(len=4) :: centroidDepType
character(len=16) :: timestamp

! CMT source parameters from inversion
integer :: cmtExponent
double precision :: cmt(6), cmtError(6)

! CMT parameters
character(len=3) :: gcmtVersion
double precision :: cmtEig(3), cmtPlunge(3), cmtAz(3)
double precision :: cmtMoment, cmtMagnitude
double precision :: cmtStr(2), cmtDip(2), cmtRak(2)

end module event_data

!--------------------------------------------------------------------------------------------------!

program main
use io
use event_data
implicit none
integer :: i, ios, luin, keepEvent

call gcmdln()

luin = 11
open(unit=luin,file=gcmtFile,status='old')

ios = 0
do while (ios.eq.0)
    do i = 1,5
        read(luin,'(A)',iostat=ios,end=101) line(i)
    enddo
    call read_event_data()
    call keep_this_event(keepEvent)
    if (keepEvent.eq.1) then
        write(*,1001) origin(1:6),centroidLon, centroidLat, centroidDep, cmtMagnitude, &
                      cmt*10.0d0**cmtExponent/1.0d7
    endif
enddo
101 continue

1001 format(I0.4,"-",I0.2,"-",I0.2,"T",I0.2,":"I0.2,":",I0.2,F10.3,F10.3,F8.1,F8.2,1P6E14.4)

end program main

!--------------------------------------------------------------------------------------------------!

subroutine keep_this_event(keepEvent)
use io, only : geoMode, lonMin, lonMax, latMin, latMax, lon0, lat0, distMax, &
               magMin, magMax, &
               jdB, jdE
use event_data, only : centroidLon, centroidLat, cmtMagnitude, origin
implicit none
integer :: keepEvent
double precision :: dist, az, dyO, jdO
keepEvent = 1

! Geographic test
if (geoMode.eq.'rect') then
    if (lonMin.le.centroidLon.and.centroidLon.le.lonMax .and. &
        latMin.le.centroidLat.and.centroidLat.le.latMax) then
        keepEvent = 1
    else
        keepEvent = 0
        return
    endif
elseif (geoMode.eq.'circle') then
    ! Returns distance in radians
    call ddistaz(dist,az,lon0,lat0,centroidLon,centroidLat)
    ! Compare distance in km
    if (dist*6371.0.le.distMax) then
        keepEvent = 1
    else
        keepEvent = 0
        return
    endif
else
    write(0,*) '!! Error: no geoMode named ',trim(geoMode)
endif

! Magnitude test
if (magMin.le.cmtMagnitude.and.cmtMagnitude.le.magMax) then
    keepEvent = 1
else
    keepEvent = 0
    return
endif

! Date test
dyO = dble(origin(3)) + dble(origin(4))/2.4d1 + dble(origin(5))/1.44d3 + dble(origin(6))/8.64d4
call date2jd(jdO,origin(1),origin(2),dyO)
! Compare number of days in Julian calendar
if (jdB.le.jdO.and.jdO.le.jdE) then
    keepEvent = 1
else
    keepEvent = 0
    return
endif

return
end subroutine keep_this_event

!-------------------------------------------------------------------------------------------------!

SUBROUTINE date2jd(jd,year,month,day)
!----
! Compute the Julian day from the calendar date
! Assumes all dates are after 1582.
! Modified from dateutil.f
!----
IMPLICIT none
integer :: year,month,yrp,mop
double precision :: day,A,B,C,D,jd

if (year.le.1582) then
    write(0,*) 'WARNING: CALCULATION MAY NOT WORK CORRECTLY FOR YEARS BEFORE 1582!'
endif

if (month.eq.1.or.month.eq.2) then
    yrp = year - 1
    mop = month + 12
else
    yrp = year
    mop = month
endif

if (year.lt.1582 .or. &
    year.eq.1582.and.month.lt.10 .or. &
    year.eq.1582.and.month.eq.10.and.day.lt.15.0d0) then
    B = 0.0d0
else
    A = floor(dble(yrp)/100.0d0)
    B = 2.0d0 - A + floor(A*0.25d0)
endif

if (yrp.lt.0) then
    C = floor(365.25d0*dble(yrp) - 0.75d0)
else
    C = floor(365.25d0*dble(yrp))
endif
D = floor(30.6001d0*(dble(mop)+1.0d0))
jd = B+C+D+day+1720994.5d0

RETURN
END

!-------------------------------------------------------------------------------------------------!

subroutine read_event_data()
use io, only : verbosity
use event_data
implicit none
! Local variables
character(len=128) :: fmt_string
character(len=80) :: string(10)
character(len=1) :: c
integer :: i, ios

! First line: Hypocenter line
! [1-4]   Hypocenter reference catalog (e.g., PDE for USGS location, ISC for
!         ISC catalog, SWE for surface-wave location, [Ekstrom, BSSA, 2006])
! [6-15]  Date of reference event
! [17-26] Time of reference event
! [28-33] Latitude
! [35-41] Longitude
! [43-47] Depth
! [49-55] Reported magnitudes, usually mb and MS
! [57-80] Geographical location (24 characters)
fmt_string = '(A4,X,A10,X,A10,X,A6,X,A6,X,A5,X,A8,X,A24)'
read(line(1),fmt_string) (string(i),i=1,8)

read(string(1),'(A)') hypoCatalog
read(string(2),'(I4,X,I2,X,I2)') origin(1),origin(2),origin(3)
read(string(3),'(I2,X,I2,X,I2,X,I1)') origin(4),origin(5),origin(6),origin(7)
read(string(4),*) hypoLat
read(string(5),*) hypoLon
read(string(6),*) hypoDep
read(string(7),*) mb, ms
read(string(8),'(A)') eventLabel

if (verbosity.ge.1) then
    print *, 'PARSED LINE 1: ',line(1)
    print *, 'hypoCatalog: ',hypoCatalog
    print *, 'hypoDate:    ',origin(1), origin(2), origin(3)
    print *, 'hypoTime:    ',origin(4), origin(5), origin(6), origin(7)
    print *, 'hypoLat:     ',hypoLat
    print *, 'hypoLon:     ',hypoLon
    print *, 'hypoDep:     ',hypoDep
    print *, 'mb:          ',mb
    print *, 'ms:          ',ms
    print *, 'eventLabel:  ',eventLabel
    print *
endif

! Second line: CMT info (1)
! [1-16]  CMT event name. This string is a unique CMT-event identifier. Older
!         events have 8-character names, current ones have 14-character names.
!         See note (1) below for the naming conventions used.
! [18-61] Data used in the CMT inversion. Three data types may be used: 
!         Long-period body waves (B), Intermediate-period surface waves (S),
!         and long-period mantle waves (M). For each data type, three values
!         are given: the number of stations used, the number of components 
!         used, and the shortest period used.
! [63-68] Type of source inverted for: "CMT: 0" - general moment tensor; 
!         "CMT: 1" - moment tensor with constraint of zero trace (standard); 
!         "CMT: 2" - double-couple source.
! [70-80] Type and duration of moment-rate function assumed in the inversion. 
!         "TRIHD" indicates a triangular moment-rate function, "BOXHD" indicates
!         a boxcar moment-rate function. The value given is half the duration
!         of the moment-rate function. This value is assumed in the inversion,
!         following a standard scaling relationship (see note (2) below),
!         and is not derived from the analysis.
fmt_string = '(A16,X,A44,X,A6,X,A11)'
read(line(2),fmt_string,iostat=ios,err=9992) (string(i),i=1,4)

read(string(1),'(A)',iostat=ios,err=9992) gcmtID
fmt_string = '(3(A2,I3,I5,F4.0,X))'
read(string(2),fmt_string,iostat=ios,err=9992) c,nstaB,ncompB,minPeriodB, &
                                      c,nstaS,ncompS,minPeriodS, &
                                      c,nstaM,ncompM,minPeriodM
read(string(3),*,iostat=ios,err=9992) c,gcmtType
read(string(4),*,iostat=ios,err=9992) stfType,stfHalfDuration
i = index(stfType,':')
stfType(i:i) = ''

9992 if (ios.ne.0) then
    write(0,*) '                123456789 123456789 123456789 123456789 123456789'
    write(0,*) 'Offending line: ',line(2)
    stop
endif

if (verbosity.ge.1) then
    print *, 'PARSED LINE 2: ',line(2)
    print *,'gcmtID:   ',gcmtID
    print *,'body:     ',nstaB,ncompB,minPeriodB
    print *,'surface:  ',nstaS,ncompS,minPeriodS
    print *,'mantle:   ',nstaM,ncompM,minPeriodM
    print *,'gcmtType: ',gcmtType
    print *,'stfType:  ',stfType
    print *,'stfHalfDuration: ',stfHalfDuration
    print *
endif

! Third line: CMT info (2)
! [1-58]  Centroid parameters determined in the inversion. Centroid time, given
!         with respect to the reference time, centroid latitude, centroid
!         longitude, and centroid depth. The value of each variable is followed
!         by its estimated standard error. See note (3) below for cases in
!         which the hypocentral coordinates are held fixed.
! [60-63] Type of depth. "FREE" indicates that the depth was a result of the
!         inversion; "FIX " that the depth was fixed and not inverted for;
!         "BDY " that the depth was fixed based on modeling of broad-band 
!         P waveforms.
! [65-80] Timestamp. This 16-character string identifies the type of analysis that
!         led to the given CMT results and, for recent events, the date and 
!         time of the analysis. This is useful to distinguish Quick CMTs ("Q-"), 
!         calculated within hours of an event, from Standard CMTs ("S-"), which 
!         are calculated later. The format for this string should not be 
!         considered fixed.
fmt_string = '(A58,X,A4,X,A16)'
read(line(3),fmt_string,iostat=ios,err=9993) (string(i),i=1,3)

fmt_string = '(A10,F8.0,F4.0,F7.0,F5.0,F8.0,F5.0,F6.0,F5.0)'
read(string(1),fmt_string,iostat=ios,err=9993) c, centroidTime, centroidTimeError, &
                                                  centroidLat, centroidLatError, &
                  centroidLon, centroidLonError, centroidDep, centroidDepError
read(string(2),'(A)',iostat=ios,err=9993) centroidDepType
read(string(3),'(A)',iostat=ios,err=9993) timestamp

9993 if (ios.ne.0) then
    write(0,*) 'Offending line: ',line(3)
    stop
endif

if (verbosity.ge.1) then
        print *, 'PARSED LINE 3: ',line(3)
    print *,'centroidTime:    ',centroidTime, centroidTimeError
    print *,'centroidLat:     ',centroidLat, centroidLatError
    print *,'centroidLon:     ',centroidLon, centroidLonError
    print *,'centroidDep:     ',centroidDep, centroidDepError
    print *,'centroidDepType: ',centroidDepType
    print *,'timestamp:       ',timestamp
    print *
endif

! Fourth line: CMT info (3)
! [1-2]   The exponent for all following moment values. For example, if the
!         exponent is given as 24, the moment values that follow, expressed in 
!         dyne-cm, should be multiplied by 10**24.
! [3-80]  The six moment-tensor elements: Mrr, Mtt, Mpp, Mrt, Mrp, Mtp, where r
!         is up, t is south, and p is east. See Aki and Richards for conversions
!         to other coordinate systems. The value of each moment-tensor
!         element is followed by its estimated standard error. See note (4)
!         below for cases in which some elements are constrained in the inversion.
fmt_string = '(A2,A78)'
read(line(4),fmt_string) (string(i),i=1,2)

read(string(1),*) cmtExponent
read(string(2),*) (cmt(i), cmtError(i), i=1,6)

if (verbosity.ge.1) then
    print *, 'PARSED LINE 4: ',line(4)
    print *,'cmtExponent: ',cmtExponent
    print *,'cmtMRR:      ',cmt(1), cmtError(1)
    print *,'cmtMTT:      ',cmt(2), cmtError(2)
    print *,'cmtMPP:      ',cmt(3), cmtError(3)
    print *,'cmtMRT:      ',cmt(4), cmtError(4)
    print *,'cmtMRP:      ',cmt(5), cmtError(5)
    print *,'cmtMTP:      ',cmt(6), cmtError(6)
    print *
endif

! Fifth line: CMT info (4)
! [1-3]   Version code. This three-character string is used to track the version 
!         of the program that generates the "ndk" file.
! [4-48]  Moment tensor expressed in its principal-axis system: eigenvalue, 
!         plunge, and azimuth of the three eigenvectors. The eigenvalue should be
!         multiplied by 10**(exponent) as given on line four.
! [50-56] Scalar moment, to be multiplied by 10**(exponent) as given on line four.
! [58-80] Strike, dip, and rake for first nodal plane of the best-double-couple 
!         mechanism, repeated for the second nodal plane. The angles are defined
!         as in Aki and Richards.
fmt_string = '(A3,X,A45,X,A6,X,A23)'
read(line(5),fmt_string,iostat=ios,err=9995) (string(i),i=1,4)

read(string(1),'(A)',iostat=ios,err=9995) gcmtVersion
read(string(2),*,iostat=ios,err=9995) (cmtEig(i),cmtPlunge(i),cmtAz(i),i=1,3)
cmtEig(1) = cmtEig(1)*10.0**cmtExponent
cmtEig(2) = cmtEig(2)*10.0**cmtExponent
cmtEig(3) = cmtEig(3)*10.0**cmtExponent
read(string(3),*,iostat=ios,err=9995) cmtMoment
cmtMoment = cmtMoment*10.0**cmtExponent
cmtMagnitude = 2.0d0/3.0d0*dlog10(cmtMoment)-10.7d0
read(string(4),*,iostat=ios,err=9995) (cmtStr(i), cmtDip(i), cmtRak(i), i=1,2)

9995 if (ios.ne.0) then
    write(0,*) 'Offending line: ',line(5)
    stop
endif

if (verbosity.ge.1) then
    print *, 'PARSED LINE 5: ',line(5)
    print *,'gcmtVersion: ',gcmtVersion
    print *,'cmtEig1:     ',cmtEig(1), cmtPlunge(1), cmtAz(1)
    print *,'cmtEig2:     ',cmtEig(2), cmtPlunge(2), cmtAz(2)
    print *,'cmtEig3:     ',cmtEig(3), cmtPlunge(3), cmtAz(3)
    print *,'cmtMoment:   ',cmtMoment
    print *,'cmtDC1:      ',cmtStr(1), cmtDip(1), cmtRak(1)
    print *,'cmtDC2:      ',cmtStr(2), cmtDip(2), cmtRak(2)
endif

return
end subroutine read_event_data

!--------------------------------------------------------------------------------------------------!

subroutine usage(string)
!----
! Print program usage statement and exit
!----
implicit none
character(len=*) :: string
integer :: string_length, stderr
stderr = 0
if (string.ne.'') then
    string_length = len(string)
    write(stderr,'(A)') trim(string)
    write(stderr,'(A)')
endif
write(stderr,'(A)') 'Usage: readGCMT ...options...'
write(stderr,'(A)')
write(stderr,'(A)') '-f file                            GCMT NDK file '// &
                                                  '(https://www.globalcmt.org/CMTfiles.html)'
write(stderr,'(A)') '-rect lonMin lonMax latMin latMax  Geographic range'
write(stderr,'(A)') '-circle lon0 lat0 dist             Distance (km) from point'
write(stderr,'(A)') '-mag magMin magMax                 Magnitude range'
write(stderr,'(A)') '-date timeBeg timeEnd              Date range (format: YYYY-MM-DDTHH:MM:SS)'
write(stderr,'(A)') '-v                                 Verbose mode'
write(stderr,*)
stop
end subroutine usage
