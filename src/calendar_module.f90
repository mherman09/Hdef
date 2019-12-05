!--------------------------------------------------------------------------------------------------!
! Module calendar
!
! Routines for manipulating date and time information.
!
!--------------------------------------------------------------------------------------------------!

module calendar

double precision, parameter, public :: yr2sec = 365.25d0*24.0d0*60.0d0*60.0d0
double precision, parameter, public :: sec2yr = 1.0d0/yr2sec

public :: parse_date
public :: date2jd
public :: jd2date

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

subroutine parse_date(date_string,date,date_format,ierr)
!----
! Read a date string into an integer array. The following date_format strings are accepted:
!     YYYY-MM-DD[THH:MM:SS[.MMM]]
!     YYYY MM DD[ HH MM SS[ MMM]]
!     YYYYMMDD[HHMMSS[MMM]]
!----

use io, only: stderr

implicit none

! Arguments
character(len=*) :: date_string, date_format
integer :: date(7), ierr

! Local variables
integer :: j, ios
double precision :: dp

! Initialize
date = 0
ierr = 0

! Read date string in specified format
if (date_format.eq.'YYYY-MM-DD') then
    if (len(trim(date_string)).eq.10) then
        read(date_string,'(I4,X,I2,X,I2)',iostat=ios) (date(j),j=1,3)
    else
        ios = 1
    endif
elseif (date_format.eq.'YYYY MM DD') then
    read(date_string,*,iostat=ios) (date(j),j=1,3)
elseif (date_format.eq.'YYYYMMDD') then
    if (len(trim(date_string)).eq.8) then
        read(date_string,'(I4,I2,I2)',iostat=ios) (date(j),j=1,3)
    else
        ios = 1
    endif
elseif (date_format.eq.'YYYY-MM-DDTHH:MM:SS') then
    if (len(trim(date_string)).ge.19) then
        read(date_string,'(I4,X,I2,X,I2,X,I2,X,I2,X,I2)',iostat=ios) (date(j),j=1,6)
    else
        ios = 1
    endif
elseif (date_format.eq.'YYYY MM DD HH MM SS') then
    read(date_string,*,iostat=ios) (date(j),j=1,5),dp
    date(6) = int(dp)
elseif (date_format.eq.'YYYYMMDDHHMMSS') then
    if (len(trim(date_string)).eq.14) then
        read(date_string,'(I4,I2,I2,I2,I2,I2)',iostat=ios) (date(j),j=1,6)
    else
        ios = 1
    endif
elseif (date_format.eq.'YYYY-MM-DDTHH:MM:SS.MMM') then
    if (len(trim(date_string)).eq.23) then
        read(date_string,'(I4,X,I2,X,I2,X,I2,X,I2,X,I2,X,I3)',iostat=ios) (date(j),j=1,7)
    else
        ios = 1
    endif
elseif (date_format.eq.'YYYY MM DD HH MM SS MMM') then
    read(date_string,*,iostat=ios) (date(j),j=1,7)
elseif (date_format.eq.'YYYYMMDDHHMMSSMMM') then
    if (len(trim(date_string)).eq.17) then
        read(date_string,'(I4,I2,I2,I2,I2,I2,I3)',iostat=ios) (date(j),j=1,7)
    else
        ios = 1
    endif
else
    write(stderr,*) 'parse_date: date_format "',trim(date_format),'" not recognized'
    ierr = 1
    return
endif

if (ios.ne.0) then
    write(stderr,*) 'parse_date: error parsing date_string "',trim(date_string),'" with format ', &
                    trim(date_format)
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine date2jd(date,julian_day)
!----
! Compute the Julian day from the calendar date. Only works for dates after 1582.
!
! Note: I am not sure where I got this algorithm from, so it does not have a reference or good
! documentation. Sorry!
!----

use io, only: stderr
implicit none

! Arguments
integer :: date(7)
double precision :: julian_day

! Local variables
integer :: year, month
double precision :: A, B, C, D, day_fraction

if (date(1).lt.1582) then
    write(stderr,*) 'date2jd: calculation may not work correctly for years before 1582'
endif

if (date(2).le.2) then
    year = date(1) - 1
    month = date(2) + 12
else
    year = date(1)
    month = date(2)
endif

if (date(1).lt.1582 .or. &
    date(1).eq.1582.and.date(2).lt.10 .or. &
    date(1).eq.1582.and.date(2).eq.10.and.date(3).lt.15) then
    B = 0.0d0
else
    A = floor(dble(year)/100.0d0)
    B = 2.0d0 - A + floor(A*0.25d0)
endif

if (year.lt.0) then
    C = floor(365.25d0*dble(year) - 0.75d0)
else
    C = floor(365.25d0*dble(year))
endif

D = floor(30.6001d0*(dble(month)+1.0d0))

day_fraction = dble(date(4))/24.0d0 + &
               dble(date(5))/24.0d0/60.0d0 + &
               dble(date(6))/24.0d0/60.0d0/60.0d0 + &
               dble(date(7))/24.0d0/60.0d0/60.0d0/1.0d3

julian_day = B + C + D + dble(date(3)) + day_fraction + 1720994.5d0

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine jd2date(julian_day,date)
!----
! Compute the calendar date corresponding to a Julian day.
!
! Note: I am not sure where I got this algorithm from, so it does not have a reference or good
! documentation. Sorry!
!----

implicit none

! Arguments
double precision :: julian_day
integer :: date(7)

! Local variables
double precision :: A, B, C, D, E, F, G, I
double precision :: real_year, real_month, real_day, real_hour, real_min, real_sec, real_ms

julian_day = julian_day + 0.5d0

F = julian_day - floor(julian_day)
I = julian_day - F
A = floor((I-1867216.25d0)/36524.25d0)

if (I.gt.2299160d0) then
    B = I + 1.0d0 + A - floor(A*0.25d0)
else
    B = I
endif

C = B + 1524.0d0
D = floor((C-122.1d0)/365.25d0)
E = floor(365.25d0*D)
G = floor((C-E)/30.6001d0)

real_day =  C-E+F-floor(30.6001d0*G)
real_hour = (real_day-floor(real_day))*24.0d0
real_min = (real_hour-floor(real_hour))*60.0d0
real_sec = (real_min-floor(real_min))*60.0d0
real_ms = (real_sec-floor(real_sec))*1000.0d0

if (G.lt.13.5d0) then
    real_month = G-1.0d0
else
    real_month = G-13.0d0
endif

if (real_month.gt.2.5d0) then
    real_year = D - 4716.0d0
else
    real_year = D - 4715.0d0
endif

! Convert real date to integer date
date(1) = int(real_year)
date(2) = int(real_month)
date(3) = int(floor(real_day))
date(4) = int(floor(real_hour))
date(5) = int(floor(real_min))
date(6) = int(floor(real_sec))
date(7) = int(floor(real_ms))

return
end subroutine

end module
