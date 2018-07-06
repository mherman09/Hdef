program main
implicit none
character(len=128) :: line
integer :: ios, luin, luout
character(len=1024) :: input_file, output_file
character(len=64) :: mode
double precision :: lon0, lat0, lon, lat, radius
double precision :: k, x, y, rho, c, arg1
double precision, parameter :: d2r = datan(1.0d0)*4.0d0/1.8d2
logical :: ex

! Parse command line
call gcmdln(mode,lon0,lat0,radius,input_file,output_file)
!write(0,*) trim(input_file)
!write(0,*) trim(output_file)
!write(0,*) trim(mode)
!write(0,*) lon0
!write(0,*) lat0
!write(0,*) radius

if (trim(input_file).eq.'stdin') then
    luin = 5
else
    inquire(file=input_file,exist=ex)
    if (.not.ex) then
        call usage('!! Error: no input file found named '//trim(input_file))
    endif
    luin = 11
    open(unit=11,file=input_file,status='old')
endif
if (trim(output_file).eq.'stdout') then
    luout = 6
else
    luout = 12
    open(unit=12,file=output_file,status='unknown')
endif

! Stereographic projection or inverse
do
    read(luin,'(A)',iostat=ios) line
    if (ios.ne.0) then
        exit
    endif
    if (mode.eq.'project') then
        read(line,*,iostat=ios) lon,lat
        k = 2.0d0*radius/(1.0d0 + dsin(lat0*d2r)*dsin(lat*d2r) + &
                              dcos(lat0*d2r)*dcos(lat*d2r)*dcos(lon*d2r-lon0*d2r))
        x = k*dcos(lat*d2r)*dsin(lon*d2r-lon0*d2r)
        y = k*(dcos(lat0*d2r)*dsin(lat*d2r)-dsin(lat0*d2r)*dcos(lat*d2r)*dcos(lon*d2r-lon0*d2r))
        write(luout,*) x, y
    elseif (mode.eq.'inverse') then
        read(line,*,iostat=ios) x,y
        rho = dsqrt(x*x+y*y)
        c = 2.0d0*datan2(rho,2.0d0*radius)
        arg1 = dcos(c)*dsin(lat0*d2r) + y*dsin(c)*dcos(lat0*d2r)/rho
        lat = datan2(arg1, dsqrt(1.0d0-arg1*arg1))
        lon = lon0 + datan2(x*dsin(c),rho*dcos(lat0*d2r)*dcos(c)-y*dsin(lat0*d2r)*dsin(c))
        write(luout,*) lon/d2r,lat/d2r
    else
        call usage('!! Error: mode must be "project" or "inverse"')
    endif
enddo

end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln(mode,lon0,lat0,radius,input_file,output_file)
implicit none
! I/O variables
character(len=*) :: mode, input_file, output_file
double precision :: lon0, lat0, radius
! Local variables
integer :: i, narg
character(len=1024) :: tag
input_file = 'stdin'
output_file = 'stdout'
mode = ''
lon0 = 1.0d3
lat0 = 1.0d3
radius = 6371.0d0
narg = command_argument_count()
if (narg.eq.0) then
    call usage('')
endif
i = 1
do while (i.le.narg)
    call get_command_argument(i,tag)
    if (trim(tag).eq.'-mode') then
        i = i + 1
        call get_command_argument(i,mode)
    elseif (trim(tag).eq.'-lon0') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) lon0
    elseif (trim(tag).eq.'-lat0') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) lat0
    elseif (trim(tag).eq.'-radius') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) radius
    elseif (trim(tag).eq.'-f') then
        i = i + 1
        call get_command_argument(i,input_file)
    elseif (trim(tag).eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file)
    else
        call usage('!! Error: no option '//trim(tag))
    endif
    i = i + 1
enddo
return
end

!--------------------------------------------------------------------------------------------------!

subroutine usage(string)
implicit none
! I/O variables
character(len=*) :: string
! Local variables
integer :: string_length
integer, parameter :: stderr = 0
if (string.ne.'') then
    string_length = len(string)
    write(stderr,*) trim(string)
    write(stderr,*)
endif
write(stderr,*) 'Usage: stereo_project ...options...'
write(stderr,*)
write(stderr,*) '-f      input_file (default: stdin)'
write(stderr,*) '-o      output_file (default: stdout)'
write(stderr,*) '-mode   project|inverse'
write(stderr,*) '-lon0   center_longitude'
write(stderr,*) '-lat0   center_latitude'
write(stderr,*) '-radius radius'
stop 1
end
