module readgcmt

character(len=512) :: gcmt_file
character(len=8) :: geo_mode
double precision :: lon_min
double precision :: lon_max
double precision :: lat_min
double precision :: lat_max
double precision :: lon0
double precision :: lat0
double precision :: dist_max
double precision :: mag_min
double precision :: mag_max
double precision :: jd_beg
double precision :: jd_end

! Hypocenter information
character(len=4) :: hypo_catalog
integer :: origin(7)
double precision :: hypo_lon
double precision :: hypo_lat
double precision :: hypo_dep
double precision :: mb
double precision :: ms
character(len=32) :: event_label

! CMT information
character(len=16) :: gcmt_id
integer :: nsta_B
integer :: ncomp_B
integer :: nsta_S
integer :: ncomp_S
integer :: nsta_M
integer :: ncomp_M
double precision :: min_period_B
double precision :: min_period_S
double precision :: min_period_M
integer :: gcmt_type
character(len=8) :: stf_type
double precision :: stf_half_duration

! CMT origin and location parameters from inversion
double precision :: centroid_time
double precision :: centroid_time_error
double precision :: centroid_lat
double precision :: centroid_lat_error
double precision :: centroid_lon
double precision :: centroid_lon_error
double precision :: centroid_dep
double precision :: centroid_dep_error
character(len=4) :: centroid_dep_type
character(len=16) :: timestamp

! CMT source parameters from inversion
integer :: cmt_exponent
double precision :: cmt(6), cmt_error(6)

! CMT parameters
character(len=3) :: gcmt_version
double precision :: cmt_eig(3)
double precision :: cmt_plunge(3)
double precision :: cmt_az(3)
double precision :: cmt_moment
double precision :: cmt_magnitude
double precision :: cmt_str(2)
double precision :: cmt_dip(2)
double precision :: cmt_rak(2)

end module

!==================================================================================================!

program main

use eq, only: dynecm2nm

use readgcmt, only: gcmt_file, &
                    origin, &
                    centroid_lon, &
                    centroid_lat, &
                    centroid_dep, &
                    cmt_magnitude, &
                    cmt, &
                    cmt_exponent

implicit none

! Local variables
integer :: i, ios, luin
logical :: keepEvent
character(len=128) :: line(5)

call gcmdln()

luin = 11
open(unit=luin,file=gcmt_file,status='old')

ios = 0
do
    do i = 1,5
        read(luin,'(A)',iostat=ios) line(i)
    enddo
    if (ios.ne.0) then
        exit
    endif

    call read_event_data(line)
    call keep_this_event(keepEvent)
    if (keepEvent) then
        write(*,1001) origin(1:6),centroid_lon, centroid_lat, centroid_dep, cmt_magnitude, &
                      cmt*10.0d0**cmt_exponent*dynecm2nm
    endif
enddo
1001 format(I0.4,"-",I0.2,"-",I0.2,"T",I0.2,":"I0.2,":",I0.2,F10.3,F10.3,F8.1,F8.2,1P6E14.4)

end program main

!-------------------------------------------------------------------------------------------------!

subroutine read_event_data(line)

use io, only: stdout, stderr, verbosity
use eq, only: dynecm2nm, mom2mag

use readgcmt, only: hypo_catalog, &
                    origin, &
                    hypo_lon, &
                    hypo_lat, &
                    hypo_dep, &
                    mb, &
                    ms, &
                    event_label, &
                    gcmt_id, &
                    nsta_B, &
                    ncomp_B, &
                    nsta_S, &
                    ncomp_S, &
                    nsta_M, &
                    ncomp_M, &
                    min_period_B, &
                    min_period_S, &
                    min_period_M, &
                    gcmt_type, &
                    stf_type, &
                    stf_half_duration, &
                    centroid_time, &
                    centroid_time_error, &
                    centroid_lat, &
                    centroid_lat_error, &
                    centroid_lon, &
                    centroid_lon_error, &
                    centroid_dep, &
                    centroid_dep_error, &
                    centroid_dep_type, &
                    timestamp, &
                    cmt_exponent, &
                    cmt, &
                    cmt_error, &
                    gcmt_version, &
                    cmt_eig, &
                    cmt_plunge, &
                    cmt_az, &
                    cmt_moment, &
                    cmt_magnitude, &
                    cmt_str, &
                    cmt_dip, &
                    cmt_rak
implicit none

! Arguments
character(len=*) :: line(*)

! Local variables
integer :: i, ios
character(len=128) :: fmt_string, string(10)
character(len=1) :: c


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

read(string(1),'(A)') hypo_catalog
read(string(2),'(I4,X,I2,X,I2)') (origin(i),i=1,3)
read(string(3),'(I2,X,I2,X,I2,X,I1)') (origin(i),i=4,7)
read(string(4),*) hypo_lat
read(string(5),*) hypo_lon
read(string(6),*) hypo_dep
read(string(7),*) mb, ms
read(string(8),'(A)') event_label

if (verbosity.ge.1) then
    write(stdout,*) 'PARSED LINE 1: ',trim(line(1))
    write(stdout,*) 'hypo_catalog: ',hypo_catalog
    write(stdout,*) 'hypo_date:    ',origin(1), origin(2), origin(3)
    write(stdout,*) 'hypo_time:    ',origin(4), origin(5), origin(6), origin(7)
    write(stdout,*) 'hypo_lat:     ',hypo_lat
    write(stdout,*) 'hypo_lon:     ',hypo_lon
    write(stdout,*) 'hypo_dep:     ',hypo_dep
    write(stdout,*) 'mb:          ',mb
    write(stdout,*) 'ms:          ',ms
    write(stdout,*) 'event_label:  ',event_label
    write(stdout,*)
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

read(string(1),'(A)',iostat=ios,err=9992) gcmt_id
fmt_string = '(3(A2,I3,I5,F4.0,X))'
read(string(2),fmt_string,iostat=ios,err=9992) c,nsta_B,ncomp_B,min_period_B, &
                                               c,nsta_S,ncomp_S,min_period_S, &
                                               c,nsta_M,ncomp_M,min_period_M
read(string(3),*,iostat=ios,err=9992) c, gcmt_type
read(string(4),*,iostat=ios,err=9992) stf_type, stf_half_duration
i = index(stf_type,':')
stf_type(i:i) = ''

9992 if (ios.ne.0) then
    write(stderr,*) '                123456789 123456789 123456789 123456789 123456789'
    write(stderr,*) 'Offending line: ',trim(line(2))
    stop
endif

if (verbosity.ge.1) then
    write(stdout,*) 'PARSED LINE 2: ',trim(line(2))
    write(stdout,*)'gcmt_id:   ',gcmt_id
    write(stdout,*)'body:     ',nsta_B,ncomp_B,min_period_B
    write(stdout,*)'surface:  ',nsta_S,ncomp_S,min_period_S
    write(stdout,*)'mantle:   ',nsta_M,ncomp_M,min_period_M
    write(stdout,*)'gcmt_type: ',gcmt_type
    write(stdout,*)'stf_type:  ',stf_type
    write(stdout,*)'stf_half_duration: ',stf_half_duration
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
read(string(1),fmt_string,iostat=ios,err=9993) c, centroid_time, centroid_time_error, &
                                                  centroid_lat, centroid_lat_error, &
                                                  centroid_lon, centroid_lon_error, &
                                                  centroid_dep, centroid_dep_error
read(string(2),'(A)',iostat=ios,err=9993) centroid_dep_type
read(string(3),'(A)',iostat=ios,err=9993) timestamp

9993 if (ios.ne.0) then
    write(stderr,*) 'Offending line: ',trim(line(3))
    stop
endif

if (verbosity.ge.1) then
    write(stdout,*) 'PARSED LINE 3: ',trim(line(3))
    write(stdout,*)'centroid_time:    ',centroid_time, centroid_time_error
    write(stdout,*)'centroid_lat:     ',centroid_lat, centroid_lat_error
    write(stdout,*)'centroid_lon:     ',centroid_lon, centroid_lon_error
    write(stdout,*)'centroid_dep:     ',centroid_dep, centroid_dep_error
    write(stdout,*)'centroid_dep_type: ',centroid_dep_type
    write(stdout,*)'timestamp:       ',timestamp
    write(stdout,*)
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

read(string(1),*) cmt_exponent
read(string(2),*) (cmt(i), cmt_error(i), i=1,6)

if (verbosity.ge.1) then
    write(stdout,*) 'PARSED LINE 4: ',trim(line(4))
    write(stdout,*)'cmt_exponent: ',cmt_exponent
    write(stdout,*)'cmtMRR:      ',cmt(1), cmt_error(1)
    write(stdout,*)'cmtMTT:      ',cmt(2), cmt_error(2)
    write(stdout,*)'cmtMPP:      ',cmt(3), cmt_error(3)
    write(stdout,*)'cmtMRT:      ',cmt(4), cmt_error(4)
    write(stdout,*)'cmtMRP:      ',cmt(5), cmt_error(5)
    write(stdout,*)'cmtMTP:      ',cmt(6), cmt_error(6)
    write(stdout,*)
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

read(string(1),'(A)',iostat=ios,err=9995) gcmt_version
read(string(2),*,iostat=ios,err=9995) (cmt_eig(i),cmt_plunge(i),cmt_az(i),i=1,3)
cmt_eig(1) = cmt_eig(1)*10.0d0**cmt_exponent
cmt_eig(2) = cmt_eig(2)*10.0d0**cmt_exponent
cmt_eig(3) = cmt_eig(3)*10.0d0**cmt_exponent
read(string(3),*,iostat=ios,err=9995) cmt_moment
cmt_moment = cmt_moment*10.0**cmt_exponent*dynecm2nm
call mom2mag(cmt_moment,cmt_magnitude)
read(string(4),*,iostat=ios,err=9995) (cmt_str(i), cmt_dip(i), cmt_rak(i), i=1,2)

9995 if (ios.ne.0) then
    write(stderr,*) 'Offending line: ',trim(line(5))
    stop
endif

if (verbosity.ge.1) then
    write(stdout,*) 'PARSED LINE 5: ',trim(line(5))
    write(stdout,*)'gcmt_version: ',gcmt_version
    write(stdout,*)'cmt_eig:     ',cmt_eig(1), cmt_plunge(1), cmt_az(1)
    write(stdout,*)'cmt_eig:     ',cmt_eig(2), cmt_plunge(2), cmt_az(2)
    write(stdout,*)'cmt_eig:     ',cmt_eig(3), cmt_plunge(3), cmt_az(3)
    write(stdout,*)'cmt_moment:   ',cmt_moment
    write(stdout,*)'cmtDC1:      ',cmt_str(1), cmt_dip(1), cmt_rak(1)
    write(stdout,*)'cmtDC2:      ',cmt_str(2), cmt_dip(2), cmt_rak(2)
    write(stdout,*)
endif

return
end subroutine read_event_data

!--------------------------------------------------------------------------------------------------!

subroutine keep_this_event(keepEvent)

use io, only: stderr
use calendar, only: date2jd
use earth, only: radius_earth_km
use geom, only: lola2distaz

use readgcmt, only: geo_mode, &
                    lon_min, &
                    lon_max, &
                    lat_min, &
                    lat_max, &
                    lon0, &
                    lat0, &
                    dist_max, &
                    jd_beg, &
                    jd_end, &
                    origin, &
                    mag_min, &
                    mag_max, &
                    jd_beg, &
                    jd_end, &
                    centroid_lat, &
                    centroid_lon, &
                    cmt_magnitude

! use io, only : geoMode, lonMin, lonMax, latMin, latMax, lon0, lat0, distMax, &
!                magMin, magMax, &
!                jdB, jdE
! use event_data, only : centroidLon, centroidLat, cmtMagnitude, origin

implicit none

! Arguments
logical :: keepEvent

! Local variables
integer :: ierr
double precision :: dist, az, jd


keepEvent = .true.

! Geographic test
if (geo_mode.eq.'rect') then
    if (lon_min.lt.0.0d0) then
        lon_min = lon_min + 360.0d0
    endif
    if (lon_max.lt.0.0d0) then
        lon_max = lon_max + 360.0d0
    endif
    if (centroid_lon.lt.0.0d0) then
        centroid_lon = centroid_lon + 360.0d0
    endif
    if (lon_min.le.centroid_lon.and.centroid_lon.le.lon_max .and. &
        lat_min.le.centroid_lat.and.centroid_lat.le.lat_max) then
        ! Keep event based on location
    else
        keepEvent = .false.
        return
    endif
elseif (geo_mode.eq.'circle') then
    call lola2distaz(lon0,lat0,centroid_lon,centroid_lat,dist,az, &
                     'radians','radians',ierr)
    if (ierr.ne.0) then
        call usage('readGCMT: error computing distance')
    endif
    if (dist*radius_earth_km.le.dist_max) then
        ! Keep event based on location
    else
        keepEvent = .false.
        return
    endif
else
    write(stderr,*) 'keep_this_event: no geo_mode named ',trim(geo_mode)
    keepEvent = .false.
    return
endif

! Magnitude test
if (mag_min.le.cmt_magnitude.and.cmt_magnitude.le.mag_max) then
    ! Keep event based on magnitude
else
    keepEvent = .false.
    return
endif

! Date test
call date2jd(origin,jd)
if (jd_beg.le.jd.and.jd.le.jd_end) then
    ! Keep event based on origin time
else
    keepEvent = .false.
    return
endif

return
end subroutine keep_this_event

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use io, only: verbosity
use trig, only: pi
use earth, only: radius_earth_km
use calendar, only: date2jd

use readgcmt, only: gcmt_file, &
                    geo_mode, &
                    lon_min, &
                    lon_max, &
                    lat_min, &
                    lat_max, &
                    lon0, &
                    lat0, &
                    dist_max, &
                    mag_min, &
                    mag_max, &
                    jd_beg, &
                    jd_end

implicit none

! Local variables
integer :: i, narg
character(len=512) :: tag

character(len=32) :: date_beg, date_end
integer :: date_beg_int(7), date_end_int(7)


! Initialize variables
gcmt_file = '/Users/mwh5316/Downloads/jan76_dec17.ndk'
geo_mode = 'rect' ! rect, circle
lon_min = -360.0d0
lon_max = 360.0d0
lat_min = -90.0d0
lat_max = 90.0d0
lon0 = 0.0d0
lat0 = 0.0d0
dist_max = radius_earth_km*pi
mag_min = 0.0
mag_max = 20.0
date_beg = '1976-01-01T00:00:00'
date_end = '2076-01-01T00:00:00'
date_beg_int = 0
date_end_int = 0

verbosity = 0

narg = command_argument_count()
if (narg.eq.0) then
    call usage('')
endif

i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    if (trim(tag).eq.'-f') then
        i = i + 1
        call get_command_argument(i,gcmt_file)

    elseif (trim(tag).eq.'-rect') then
        geo_mode = 'rect'
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) lon_min
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) lon_max
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) lat_min
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) lat_max

    elseif (trim(tag).eq.'-circle') then
        geo_mode = 'circle'
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) lon0
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) lat0
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) dist_max

    elseif (trim(tag).eq.'-mag') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) mag_min
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) mag_max

    elseif (trim(tag).eq.'-date') then
        i = i + 1
        call get_command_argument(i,date_beg)
        i = i + 1
        call get_command_argument(i,date_end)

    elseif (trim(tag).eq.'-v') then
        verbosity = 1

    else
        call usage('readGCMT: no option '//trim(tag))
    endif

    i = i + 1
enddo

! timeBeg: YYYY-MM-DDTHH:MM:SS
read(date_beg(1:4),  '(I4)') date_beg_int(1)
read(date_beg(6:7),  '(I2)') date_beg_int(2)
read(date_beg(9:10), '(I2)') date_beg_int(3)
read(date_beg(12:13),'(I2)') date_beg_int(4)
read(date_beg(15:16),'(I2)') date_beg_int(5)
read(date_beg(18:19),'(I2)') date_beg_int(6)
call date2jd(date_beg_int,jd_beg)

! timeEnd: YYYY-MM-DDTHH:MM:SS
read(date_end(1:4),  '(I4)') date_end_int(1)
read(date_end(6:7),  '(I2)') date_end_int(2)
read(date_end(9:10), '(I2)') date_end_int(3)
read(date_end(12:13),'(I2)') date_end_int(4)
read(date_end(15:16),'(I2)') date_end_int(5)
read(date_end(18:19),'(I2)') date_end_int(6)
call date2jd(date_end_int,jd_end)

return
end subroutine gcmdln

!--------------------------------------------------------------------------------------------------!

subroutine usage(string)

use io, only: stderr

implicit none

character(len=*) :: string

if (string.ne.'') then
    write(stderr,*) trim(string)
    write(stderr,*)
endif
write(stderr,*) 'Usage: readGCMT ...options...'
write(stderr,*)
write(stderr,*) '-f file                                GCMT NDK file '// &
                                                      '(https://www.globalcmt.org/CMTfiles.html)'
write(stderr,*) '-rect lon_min lon_max lat_min lat_max  Geographic range'
write(stderr,*) '-circle lon0 lat0 dist                 Distance (km) from point'
write(stderr,*) '-mag mag_min mag_max                   Magnitude range'
write(stderr,*) '-date date_beg date_end                Date range (format: YYYY-MM-DDTHH:MM:SS)'
write(stderr,*)

stop
end subroutine usage
