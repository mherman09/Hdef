module okada92
!----
! Routines to compute deformation from point and finite sources in an elastic half-space.
!
! Okada, Y. (1992). Internal deformation due to shear and tensile faults in a half-space.
! Bulletin of the Seismological Society of America 82, 1018-1040.
!
! - All variables are in SI units
! - Angles input into subroutines are in degrees
! - Coordinate system is relative to the fault orientation:
!   - X points in the along-strike direction
!   - Y points in the horizontal up-dip direction
!   - Depths are defined positive down
!----

! Half-space variables
double precision :: a, CA1, CA2, CB, CC

! Fault dip variables
double precision :: sd, cd, s2d, c2d, cdcd, sdsd, cdsd

! Basic geometric variables
double precision :: x, y, z, c, d

! Derived geometric variables
double precision :: p, q, s, t, xx, xy, yy, dd, xd, xc, xq, xs, xt, yq, ys, yt, dq, pq, qq
double precision :: R, R2, R3, R4, R5, R7, Rd
double precision :: A3, A5, A7, B3, B5, B7, C3, C5, C7
double precision :: I1, I2, I3, I4, I5, J1, J2, J3, J4, J5, J6, K1, K2, K3, K4
double precision :: U2, V2, W2, U3, V3, W3
double precision :: eta_vec(2), ksi_vec(2), ksi, eta
double precision :: cbar, dbar, ybar
double precision :: X11, X32, X53, Y11, Y32, Y53, Y0, Z32, Z53, Z0
double precision :: E2, F2, G2, H2, P2, Q2, E3, F3, G3, H3, P3, Q3
double precision :: D11, Re, Rk, logRe, logRk, TH

! Integration constants for rectangular faults, stored as an array:
!    chinnery_factor = [  1 -1 ]
!                      [ -1  1 ]
double precision :: chinnery_factor(2,2)=reshape((/1.0d0,-1.0d0,-1.0d0,1.0d0/),(/2,2/))

! Singular flag
logical :: isSingular

! Private variables and routines are only public for unit testing
#ifdef UNIT_TEST
  public :: a, CA1, CA2, CB, CC
  public :: sd, cd, s2d, c2d, cdcd, sdsd, cdsd
  public :: x, y, z, c, d
  public :: p, q, s, t, xx, xy, yy, dd, xd, xc, xq, xs, xt, yq, ys, yt, dq, pq, qq
  public :: R, R2, R3, R4, R5, R7, Rd
  public :: A3, A5, A7, B3, B5, B7, C3, C5, C7
  public :: I1, I2, I3, I4, I5, J1, J2, J3, J4, J5, J6, K1, K2, K3, K4
  public :: U2, V2, W2, U3, V3, W3
  public :: eta_vec, ksi_vec, ksi, eta
  public :: cbar, dbar, ybar
  public :: X11, X32, X53, Y11, Y32, Y53, Y0, Z32, Z53, Z0
  public :: E2, F2, G2, H2, P2, Q2, E3, F3, G3, H3, P3, Q3
  public :: D11, Re, Rk, logRe, logRk, TH
  public :: chinnery_factor
  public :: isSingular
  public :: halfspace_vars
  public :: dip_vars
  public :: pt_src_vars
  public :: rect_src_coords
  public :: rect_src_vars
  public :: check_singular_pt
  public :: check_singular_rect
#else
  private :: a, CA1, CA2, CB, CC
  private :: sd, cd, s2d, c2d, cdcd, sdsd, cdsd
  private :: x, y, z, c, d
  private :: p, q, s, t, xx, xy, yy, dd, xd, xc, xq, xs, xt, yq, ys, yt, dq, pq, qq
  private :: R, R2, R3, R4, R5, R7, Rd
  private :: A3, A5, A7, B3, B5, B7, C3, C5, C7
  private :: I1, I2, I3, I4, I5, J1, J2, J3, J4, J5, J6, K1, K2, K3, K4
  private :: U2, V2, W2, U3, V3, W3
  private :: eta_vec, ksi_vec, ksi, eta
  private :: cbar, dbar, ybar
  private :: X11, X32, X53, Y11, Y32, Y53, Y0, Z32, Z53, Z0
  private :: E2, F2, G2, H2, P2, Q2, E3, F3, G3, H3, P3, Q3
  private :: D11, Re, Rk, logRe, logRk, TH
  private :: chinnery_factor
  private :: isSingular
  private :: halfspace_vars
  private :: dip_vars
  private :: pt_src_vars
  private :: rect_src_coords
  private :: rect_src_vars
  private :: check_singular_pt
  private :: check_singular_rect
#endif

! Public routines
public :: o92_pt_disp
public :: o92_pt_strain
public :: o92_pt_partials
public :: o92_rect_disp
public :: o92_rect_strain
public :: o92_rect_partials
public :: test_halfspace_vars
public :: test_dip_vars
public :: test_pt_src_vars
public :: test_rect_src_vars
public :: test_rect_disp_components
public :: test_rect_partial_components

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!


subroutine halfspace_vars(lambda,mu)
!----
! Coefficients in displacement and strain equations, where lambda is the Lame parameter and mu is
! the shear modulus, both in Pa.
!----

implicit none

! Arguments
double precision :: lambda, mu

a   = (lambda+mu)/(lambda+2.0d0*mu)
CA1 = (1.0d0-a)*0.5d0
CA2 = a*0.5d0
CB  = (1.0d0-a)/a
CC  = 1.0d0-a

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine dip_vars(dip)
!----
! Trigonometric functions of the fault dip (in degrees) found in the deformation equations.
!----

use trig, only: d2r

implicit none

! Arguments
double precision :: dip

! Local variables
double precision, parameter :: dip_max = 89.997d0

! To avoid divide by zero problems, avoid dip of 90 degress
if (dip.gt.dip_max) then
  dip = dip_max
endif

sd = dsin(dip*d2r)
cd = dcos(dip*d2r)
s2d = dsin(2.0d0*dip*d2r)
c2d = dcos(2.0d0*dip*d2r)
cdcd = cd*cd
sdsd = sd*sd
cdsd = sd*cd

return
end subroutine

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!------------------------------------- POINT SOURCE -----------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine o92_pt_disp(disp,sta_coord,evdp,dip,moment,lambda,shear_modulus)
!----
! Compute displacement vector at a point in an elastic half-space due to a point source dislocation
! Subroutine arguments:
!     disp(3): output displacement vector (m)
!     sta_coord(3): input station location(m); x=along-strike, y=hor-up-dip, z=depth
!     evdp: source depth (m), positive down
!     dip: source fault dip (degrees)
!     moment(4): strike-slip, dip-slip, tensile, and volume moments (N-m = shear_modulus*area*slip)
!     lambda, shear_modulus: half-space parameters (Pa)
!----

use io, only: stderr
use trig, only: pi

implicit none

! Arguments
double precision :: disp(3), sta_coord(3), evdp, dip, moment(4), lambda, shear_modulus

! Local variables
double precision :: Css, Cds, Cts, Cvol, u_ss(3), u_ds(3), u_ts(3), u_vol(3)

! write(0,*) 'o92_pt_disp: starting'

! Initialize displacements
disp = 0.0d0

! If station is above surface, exit with warning
if (sta_coord(3).lt.0.0d0) then
    write(stderr,*) 'o92_pt_disp: station depth less than zero, setting displacement to zero'
    return
endif

! Calculate some parameters for fault geometry and elastic half-space
call dip_vars(dip)
call halfspace_vars(lambda,shear_modulus)

! Scaling factors
Css = moment(1)/(2.0d0*pi*shear_modulus)
Cds = moment(2)/(2.0d0*pi*shear_modulus)
Cts = moment(3)/(2.0d0*pi*shear_modulus)
Cvol = moment(4)/(2.0d0*pi*shear_modulus)

! Start calculation with values that depend on station depth (-sta_coord(3) saved as z in subroutine)
call pt_src_vars(sta_coord(1),sta_coord(2),-sta_coord(3),evdp)

! Check for singular solution when R=0
call check_singular_pt(isSingular,R)
if (isSingular) then
    write(stderr,*) 'o92_pt_disp: solution is singular, setting displacement to zero'
    return
endif

u_ss = 0.0d0
u_ds = 0.0d0
u_ts = 0.0d0
u_vol = 0.0d0
if (dabs(Css).gt.1.0d-6) then
    u_ss = uA_ss() + uB_ss() + z*uC_ss()
endif
if (dabs(Cds).gt.1.0d-6) then
    u_ds = uA_ds() + uB_ds() + z*uC_ds()
endif
if (dabs(Cts).gt.1.0d-6) then
    u_ts = uA_ts() + uB_ts() + z*uC_ts()
endif
if (dabs(Cvol).gt.1.0d-6) then
    u_vol = uA_vol() + uB_vol() + z*uC_vol()
endif
disp = disp + Css*u_ss + Cds*u_ds + Cts*u_ts + Cvol*u_vol

! Add term that depends on negative station depth
call pt_src_vars(sta_coord(1),sta_coord(2),sta_coord(3),evdp)

! Check for singular solution when R=0
call check_singular_pt(isSingular,R)
if (isSingular) then
    write(stderr,*) 'o92_pt_disp: solution is singular, setting displacement to zero'
    return
endif

u_ss = 0.0d0
u_ds = 0.0d0
u_ts = 0.0d0
u_vol = 0.0d0
if (dabs(Css).gt.1.0d-6) then
    u_ss = uA_ss()
endif
if (dabs(Cds).gt.1.0d-6) then
    u_ds = uA_ds()
endif
if (dabs(Cts).gt.1.0d-6) then
    u_ts = uA_ts()
endif
if (dabs(Cvol).gt.1.0d-6) then
    u_vol = uA_vol()
endif
disp = disp - (Css*u_ss + Cds*u_ds + Cts*u_ts + Cvol*u_vol)

! write(0,*) 'o92_pt_disp: finished'

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine o92_pt_strain(strain,sta_coord,evdp,dip,moment,lambda,shear_modulus)
!----
! Compute strain tensor at a point in an elastic half-space due to a point source dislocation
! Subroutine arguments:
!     strain(3,3): output strain tensor
!     sta_coord(3): input station location(m); x=along-strike, y=hor-up-dip, z=depth
!     evdp: source depth (m), positive down
!     dip: source fault dip (degrees)
!     moment(4): strike-slip, dip-slip, tensile, and volume moments (N-m = shear_modulus*area*slip)
!     lambda, shear_modulus: half-space parameters (Pa)
!----

use io, only: stderr

implicit none

! Arguments
double precision :: strain(3,3), sta_coord(3), evdp, dip, moment(4), lambda, shear_modulus

! Local variables
double precision :: uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz

! write(0,*) 'o92_pt_strain: starting'

! Initialize strains
strain = 0.0d0

! If station is above surface, exit with warning
if (sta_coord(3).lt.0.0d0) then
    write(stderr,*) 'o92_pt_strain: station depth less than zero; setting strain to zero'
    return
endif

! Partial derivatives of displacement equations
call o92_pt_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz, &
                     sta_coord,evdp,dip,moment,lambda,shear_modulus)

! Strains
strain(1,1) = uxx
strain(1,2) = 1.0d0/2.0d0*(uxy+uyx)
strain(1,3) = 1.0d0/2.0d0*(uxz+uzx)
strain(2,1) = strain(1,2)
strain(2,2) = uyy
strain(2,3) = 1.0d0/2.0d0*(uyz+uzy)
strain(3,1) = strain(1,3)
strain(3,2) = strain(2,3)
strain(3,3) = uzz

! write(0,*) 'o92_pt_strain: finished'

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine o92_pt_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz, &
                           sta_coord,evdp,dip,moment,lambda,shear_modulus)
!----
! Compute partial derivatives at a point in an elastic half-space due to a point source dislocation
!----

use io, only: stderr
use trig, only: pi

implicit none

! Arguments
double precision :: uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz, &
                    sta_coord(3), evdp, dip, moment(4), lambda, shear_modulus

! Local variables
double precision :: Css, Cds, Cts, Cvol
double precision :: dudx_ss(3), dudy_ss(3), dudz_ss(3), u_ss(3)
double precision :: dudx_ds(3), dudy_ds(3), dudz_ds(3), u_ds(3)
double precision :: dudx_ts(3), dudy_ts(3), dudz_ts(3), u_ts(3)
double precision :: dudx_vol(3), dudy_vol(3), dudz_vol(3), u_vol(3)

! write(0,*) 'o92_pt_partials: starting'

! Initialize partial derivatives
uxx = 0.0d0
uyx = 0.0d0
uzx = 0.0d0
uxy = 0.0d0
uyy = 0.0d0
uzy = 0.0d0
uxz = 0.0d0
uyz = 0.0d0
uzz = 0.0d0

! Calculate some parameters for fault geometry and elastic half-space
call dip_vars(dip)
call halfspace_vars(lambda,shear_modulus)

! Scaling factors
Css = moment(1)/(2.0d0*pi*shear_modulus)
Cds = moment(2)/(2.0d0*pi*shear_modulus)
Cts = moment(3)/(2.0d0*pi*shear_modulus)
Cvol = moment(4)/(2.0d0*pi*shear_modulus)

! Start calculation with values that depend on station depth (-sta_coord(3) saved as z in subroutine)
call pt_src_vars(sta_coord(1),sta_coord(2),-sta_coord(3),evdp)

! Check for singular solution when R=0
call check_singular_pt(isSingular,R)
if (isSingular) then
    write(stderr,*) 'o92_pt_partials: solution is singular, setting partial derivatives to zero'
    return
endif

dudx_ss = 0.0d0
dudy_ss = 0.0d0
dudz_ss = 0.0d0
dudx_ds = 0.0d0
dudy_ds = 0.0d0
dudz_ds = 0.0d0
dudx_ts = 0.0d0
dudy_ts = 0.0d0
dudz_ts = 0.0d0
dudx_vol = 0.0d0
dudy_vol = 0.0d0
dudz_vol = 0.0d0
u_ss = 0.0d0
u_ds = 0.0d0
u_ts = 0.0d0
u_vol = 0.0d0
if (dabs(Css).gt.1.0d-6) then
    dudx_ss = duAdx_ss() + duBdx_ss() + z*duCdx_ss()
    dudy_ss = duAdy_ss() + duBdy_ss() + z*duCdy_ss()
    dudz_ss = duAdz_ss() + duBdz_ss() + z*duCdz_ss()
    u_ss = uC_ss()
endif
if (dabs(Cds).gt.1.0d-6) then
    dudx_ds = duAdx_ds() + duBdx_ds() + z*duCdx_ds()
    dudy_ds = duAdy_ds() + duBdy_ds() + z*duCdy_ds()
    dudz_ds = duAdz_ds() + duBdz_ds() + z*duCdz_ds()
    u_ds = uC_ds()
endif
if (dabs(Cts).gt.1.0d-6) then
    dudx_ts = duAdx_ts() + duBdx_ts() + z*duCdx_ts()
    dudy_ts = duAdy_ts() + duBdy_ts() + z*duCdy_ts()
    dudz_ts = duAdz_ts() + duBdz_ts() + z*duCdz_ts()
    u_ts = uC_ts()
endif
if (dabs(Cvol).gt.1.0d-6) then
    dudx_vol = duAdx_vol() + duBdx_vol() + z*duCdx_vol()
    dudy_vol = duAdy_vol() + duBdy_vol() + z*duCdy_vol()
    dudz_vol = duAdz_vol() + duBdz_vol() + z*duCdz_vol()
    u_vol = uC_vol()
endif

! Add term that depends on negative station depth
call pt_src_vars(sta_coord(1),sta_coord(2),sta_coord(3),evdp)

! Check for singular solution when R=0
call check_singular_pt(isSingular,R)
if (isSingular) then
    write(stderr,*) 'o92_pt_partials: solution is singular, setting partial derivatives to zero'
    return
endif

if (dabs(Css).gt.1.0d-6) then
    dudx_ss = dudx_ss - duAdx_ss()
    dudy_ss = dudy_ss - duAdy_ss()
    dudz_ss = dudz_ss + duAdz_ss()
endif
if (dabs(Cds).gt.1.0d-6) then
    dudx_ds = dudx_ds - duAdx_ds()
    dudy_ds = dudy_ds - duAdy_ds()
    dudz_ds = dudz_ds + duAdz_ds()
endif
if (dabs(Cts).gt.1.0d-6) then
    dudx_ts = dudx_ts - duAdx_ts()
    dudy_ts = dudy_ts - duAdy_ts()
    dudz_ts = dudz_ts + duAdz_ts()
endif
if (dabs(Cvol).gt.1.0d-6) then
    dudx_vol = dudx_vol - duAdx_vol()
    dudy_vol = dudy_vol - duAdy_vol()
    dudz_vol = dudz_vol + duAdz_vol()
endif

! Total partial derivatives
uxx = Css*dudx_ss(1) + Cds*dudx_ds(1) + Cts*dudx_ts(1) + Cvol*dudx_vol(1)
uxy = Css*dudy_ss(1) + Cds*dudy_ds(1) + Cts*dudy_ts(1) + Cvol*dudy_vol(1)
uxz = Css*dudz_ss(1) + Cds*dudz_ds(1) + Cts*dudz_ts(1) + Cvol*dudz_vol(1) + &
      Css*u_ss(1) + Cds*u_ds(1) + Cts*u_ts(1) + Cvol*u_vol(1)
uyx = Css*dudx_ss(2) + Cds*dudx_ds(2) + Cts*dudx_ts(2) + Cvol*dudx_vol(2)
uyy = Css*dudy_ss(2) + Cds*dudy_ds(2) + Cts*dudy_ts(2) + Cvol*dudy_vol(2)
uyz = Css*dudz_ss(2) + Cds*dudz_ds(2) + Cts*dudz_ts(2) + Cvol*dudz_vol(2) + &
      Css*u_ss(2) + Cds*u_ds(2) + Cts*u_ts(2) + Cvol*u_vol(2)
uzx = Css*dudx_ss(3) + Cds*dudx_ds(3) + Cts*dudx_ts(3) + Cvol*dudx_vol(3)
uzy = Css*dudy_ss(3) + Cds*dudy_ds(3) + Cts*dudy_ts(3) + Cvol*dudy_vol(3)
uzz = Css*dudz_ss(3) + Cds*dudz_ds(3) + Cts*dudz_ts(3) + Cvol*dudz_vol(3) + &
      Css*u_ss(3) + Cds*u_ds(3) + Cts*u_ts(3) + Cvol*u_vol(3)

! write(0,*) 'o92_pt_partials: finished'

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine pt_src_vars(xin,yin,zin,cin)
!----
! Calculate variables in point source equations from x, y, z (receiver location), and c (source
! depth), plus fault dip
!----

implicit none

! Arguments
double precision :: xin, yin, zin, cin

! Local variables
double precision :: Rd2, Rd3, Rd4

! write(0,*) 'pt_src_vars: starting'

! Save passed arguments to module variables
x = xin
y = yin
z = zin
c = cin
d = c - z

xx = x*x
xy = x*y
yy = y*y
dd = d*d
xd = x*d
xc = x*c

! R = distance between source and station
R  = dsqrt(xx+yy+dd)
R2 = R*R
R3 = R2*R
R4 = R3*R
R5 = R4*R
R7 = R5*R2
Rd = R + d
Rd2 = Rd*Rd
Rd3 = Rd2*Rd
Rd4 = Rd3*Rd

! p, q, s, t = distances in rotated reference frame
! p = fault dip-parallel distance
! q = fault dip-perpendicular distance
p = y*cd + d*sd
q = y*sd - d*cd
s = p*sd + q*cd
t = p*cd - q*sd
xq = x*q
xs = x*s
xt = x*t
yq = y*q
ys = y*s
yt = y*t
dq = d*q
pq = p*q
qq = q*q

! Other variables
A3 = 1.0d0 - 3.0d0*xx/R2
A5 = 1.0d0 - 5.0d0*xx/R2
A7 = 1.0d0 - 7.0d0*xx/R2
B3 = 1.0d0 - 3.0d0*yy/R2
B5 = 1.0d0 - 5.0d0*yy/R2
B7 = 1.0d0 - 7.0d0*yy/R2
C3 = 1.0d0 - 3.0d0*dd/R2
C5 = 1.0d0 - 5.0d0*dd/R2
C7 = 1.0d0 - 7.0d0*dd/R2

I1 =  y*(1.0d0/(R*Rd2) - xx*(3.0d0*R+d)/(R3*Rd3))
I2 =  x*(1.0d0/(R*Rd2) - yy*(3.0d0*R+d)/(R3*Rd3))
I3 =  x/R3 - I2
I4 = -xy*(2.0d0*R+d)/(R3*Rd2)
I5 =  1.0d0/(R*Rd) - xx*(2.0d0*R+d)/(R3*Rd2)

J1 = -3.0d0*xy*((3.0d0*R+d)/(R3*Rd3) - xx*(5.0d0*R2+4.0d0*R*d+dd)/(R5*Rd4))
J2 = 1.0d0/R3 - 3.0d0/(R*Rd2) + 3.0d0*xx*yy*(5.0d0*R2+4.0d0*R*d+dd)/(R5*Rd4)
J3 = A3/R3 - J2
J4 = -3.0d0*xy/R5 - J1

K1 = -y*((2.0d0*R+d)/(R3*Rd2) - xx*(8.0d0*R2+9.0d0*R*d+3.0d0*dd)/(R5*Rd3))
K2 = -x*((2.0d0*R+d)/(R3*Rd2) - yy*(8.0d0*R2+9.0d0*R*d+3.0d0*dd)/(R5*Rd3))
K3 = -3.0d0*x*d/R5 - K2

U2 = sd - 5.0d0*yq/R2
V2 = s  - 5.0d0*y*pq/R2
W2 = sd + U2
U3 = cd + 5.0d0*dq/R2
V3 = t  + 5.0d0*d*pq/R2
W3 = cd + U3

! write(0,*) 'pt_src_vars: finished'

return
end subroutine

!--------------------------------------------------------------------------------------------------!
! COMPONENTS OF DISPLACEMENT FOR A POINT SOURCE DISLOCATION - EDIT AT YOUR OWN RISK!               !
!--------------------------------------------------------------------------------------------------!

function uA_ss()
implicit none
double precision :: uA_ss(3)
uA_ss(1) =  CA1*q/R3    + CA2*3.0d0*xx*q/R5
uA_ss(2) =  CA1*x*sd/R3 + CA2*3.0d0*xy*q/R5
uA_ss(3) = -CA1*x*cd/R3 + CA2*3.0d0*xd*q/R5
return
end function

function uA_ds()
implicit none
double precision :: uA_ds(3)
uA_ds(1) =             CA2*3.0d0*x*pq/R5
uA_ds(2) =  CA1*s/R3 + CA2*3.0d0*y*pq/R5
uA_ds(3) = -CA1*t/R3 + CA2*3.0d0*d*pq/R5
return
end function

function uA_ts()
implicit none
double precision :: uA_ts(3)
uA_ts(1) = CA1*x/R3 - CA2*3.0d0*x*qq/R5
uA_ts(2) = CA1*t/R3 - CA2*3.0d0*y*qq/R5
uA_ts(3) = CA1*s/R3 - CA2*3.0d0*d*qq/R5
return
end function

function uA_vol()
implicit none
double precision :: uA_vol(3)
uA_vol(1) = -CA1*x/R3
uA_vol(2) = -CA1*y/R3
uA_vol(3) = -CA1*d/R3
return
end function

function uB_ss()
implicit none
double precision :: uB_ss(3)
uB_ss(1) = -3.0d0*xx*q/R5 - CB*I1*sd
uB_ss(2) = -3.0d0*xy*q/R5 - CB*I2*sd
uB_ss(3) = -3.0d0*xc*q/R5 - CB*I4*sd
return
end function

function uB_ds()
implicit none
double precision :: uB_ds(3)
uB_ds(1) = -3.0d0*x*pq/R5 + CB*I3*cdsd
uB_ds(2) = -3.0d0*y*pq/R5 + CB*I1*cdsd
uB_ds(3) = -3.0d0*c*pq/R5 + CB*I5*cdsd
return
end function

function uB_ts()
implicit none
double precision :: uB_ts(3)
uB_ts(1) = 3.0d0*x*qq/R5 - CB*I3*sdsd
uB_ts(2) = 3.0d0*y*qq/R5 - CB*I1*sdsd
uB_ts(3) = 3.0d0*c*qq/R5 - CB*I5*sdsd
return
end function

function uB_vol()
implicit none
double precision :: uB_vol(3)
uB_vol(1) = CB*x/R3
uB_vol(2) = CB*y/R3
uB_vol(3) = CB*d/R3
return
end function

function uC_ss()
implicit none
double precision :: uC_ss(3)
uC_ss(1) = -CC*A3*cd/R3 + a*3.0d0*c*q*A5/R5
uC_ss(2) =  CC*3.0d0*xy*cd/R5 + a*3.0d0*xc*(sd-5.0d0*y*q/R2)/R5
uC_ss(3) = -CC*3.0d0*xy*sd/R5 + a*3.0d0*xc*(cd+5.0d0*d*q/R2)/R5
return
end function

function uC_ds()
implicit none
double precision :: uC_ds(3)
uC_ds(1) =  CC*3.0d0*x*t/R5 - a*15.0d0*xc*pq/R7
uC_ds(2) = -CC*(c2d-3.0d0*yt/R2)/R3 + a*3.0d0*c*(s-5.0d0*y*pq/R2)/R5
uC_ds(3) = -CC*A3*cdsd/R3 + a*3.0d0*c*(t+5.0d0*d*pq/R2)/R5
return
end function

function uC_ts()
implicit none
double precision :: uC_ts(3)
uC_ts(1) = -CC*3.0d0*x*s/R5 + a*15.0d0*c*x*qq/R7 - a*3.0d0*x*z/R5
uC_ts(2) =  CC*(s2d-3.0d0*ys/R2)/R3 + a*3.0d0*c*(t-y+5.0d0*y*qq/R2)/R5 - a*3.0d0*y*z/R5
uC_ts(3) = -CC*(1.0d0-A3*sdsd)/R3 - a*3.0d0*c*(s-d+5.0d0*d*qq/R2)/R5 + a*3.0d0*d*z/R5
return
end function

function uC_vol()
implicit none
double precision :: uC_vol(3)
uC_vol(1) = CC*3.0d0*x*d/R5
uC_vol(2) = CC*3.0d0*y*d/R5
uC_vol(3) = CC*C3/R3
return
end function

!--------------------------------------------------------------------------------------------------!
! COMPONENTS OF PARTIAL DERIVATIVES FOR A POINT SOURCE DISLOCATION - EDIT AT YOUR OWN RISK!        !
!--------------------------------------------------------------------------------------------------!

function duAdx_ss()
implicit none
double precision :: duAdx_ss(3)
duAdx_ss(1) = -CA1*3.0d0*xq/R5 + CA2*3.0d0*xq*(1.0d0+A5)/R5
duAdx_ss(2) =  CA1*A3*sd/R3     + CA2*3.0d0*yq*A5/R5
duAdx_ss(3) = -CA1*A3*cd/R3     + CA2*3.0d0*dq*A5/R5
return
end function

function duAdx_ds()
implicit none
double precision :: duAdx_ds(3)
duAdx_ds(1) =                     CA2*3.0d0*pq*A5/R5
duAdx_ds(2) = -CA1*3.0d0*xs/R5 - CA2*15.0d0*xy*pq/R7
duAdx_ds(3) =  CA1*3.0d0*xt/R5 - CA2*15.0d0*xd*pq/R7
return
end function

function duAdx_ts()
implicit none
double precision :: duAdx_ts(3)
duAdx_ts(1) =  CA1*A3/R3       - CA2*3.0d0*qq*A5/R5
duAdx_ts(2) = -CA1*3.0d0*xt/R5 + CA2*15.0d0*xy*qq/R7
duAdx_ts(3) = -CA1*3.0d0*xs/R5 + CA2*15.0d0*xd*qq/R7
return
end function

function duAdx_vol()
implicit none
double precision :: duAdx_vol(3)
duAdx_vol(1) = -CA1*A3/R3
duAdx_vol(2) =  CA1*3.0d0*xy/R5
duAdx_vol(3) =  CA1*3.0d0*xd/R5
return
end function

function duBdx_ss()
implicit none
double precision :: duBdx_ss(3)
duBdx_ss(1) = -3.0d0*xq*(1.0d0+A5)/R5 - CB*J1*sd
duBdx_ss(2) = -3.0d0*yq*A5/R5         - CB*J2*sd
duBdx_ss(3) = -3.0d0*c*q*A5/R5        - CB*K1*sd
return
end function

function duBdx_ds()
implicit none
double precision :: duBdx_ds(3)
duBdx_ds(1) = -3.0d0 *pq*A5/R5 + CB*J3*cdsd
duBdx_ds(2) =  15.0d0*xy*pq/R7 + CB*J1*cdsd
duBdx_ds(3) =  15.0d0*xc*pq/R7 + CB*K3*cdsd
return
end function

function duBdx_ts()
implicit none
double precision :: duBdx_ts(3)
duBdx_ts(1) =   3.0d0*qq*A5/R5 - CB*J3*sdsd
duBdx_ts(2) = -15.0d0*xy*qq/R7 - CB*J1*sdsd
duBdx_ts(3) = -15.0d0*xc*qq/R7 - CB*K3*sdsd
return
end function

function duBdx_vol()
implicit none
double precision :: duBdx_vol(3)
duBdx_vol(1) =  CB*A3/R3
duBdx_vol(2) = -CB*3.0d0*xy/R5
duBdx_vol(3) = -CB*3.0d0*xd/R5
return
end function

function duCdx_ss()
implicit none
double precision :: duCdx_ss(3)
duCdx_ss(1) =  CC*3.0d0*x*(2.0d0+A5)*cd/R5 - a*15.0d0*c*xq*(2.0d0+A7)/R7
duCdx_ss(2) =  CC*3.0d0*y*A5*cd/R5 + a*3.0d0*c*(A5*sd-5.0d0*yq*A7/R2)/R5
duCdx_ss(3) = -CC*3.0d0*y*A5*sd/R5 + a*3.0d0*c*(A5*cd+5.0d0*dq*A7/R2)/R5
return
end function

function duCdx_ds()
implicit none
double precision :: duCdx_ds(3)
duCdx_ds(1) =  CC*3.0d0*t*A5/R5 - a*15.0d0*c*pq*A7/R7
duCdx_ds(2) =  CC*3.0d0*x*(c2d-5.0d0*yt/R2)/R5 - a*15.0d0*xc*(s-7.0d0*y*pq/R2)/R7
duCdx_ds(3) =  CC*3.0d0*x*(2.0d0+A5)*sd*cd/R5 - a*15.0d0*xc*(t+7.0d0*d*pq/R2)/R7
return
end function

function duCdx_ts()
implicit none
double precision :: duCdx_ts(3)
duCdx_ts(1) = -CC*3.0d0*s*A5/R5 + a*15.0d0*c*qq*A7/R7 - a*3.0d0*z*A5/R5
duCdx_ts(2) = -CC*3.0d0*x*(s2d-5.0d0*ys/R2)/R5 - a*15.0d0*xc*(t-y+7.0d0*y*qq/R2)/R7 + &
                 a*15.0d0*xy*z/R7
duCdx_ts(3) =  CC*3.0d0*x*(1.0d0-(2.0d0+A5)*sdsd)/R5 + a*15.0d0*xc*(s-d+7.0d0*d*qq/R2)/R7 - &
                 a*15.0d0*xd*z/R7
return
end function

function duCdx_vol()
implicit none
double precision :: duCdx_vol(3)
duCdx_vol(1) =  CC*3.0d0*d*A5/R5
duCdx_vol(2) = -CC*15.0d0*xy*d/R7
duCdx_vol(3) = -CC*3.0d0*x*C5/R5
return
end function


function duAdy_ss()
implicit none
double precision :: duAdy_ss(3)
duAdy_ss(1) =  CA1*(sd-3.0d0*yq/R2)/R3 + CA2*3.0d0*xx*U2/R5
duAdy_ss(2) = -CA1*3.0d0*xy*sd/R5      + CA2*3.0d0*x*(y*U2+q)/R5
duAdy_ss(3) =  CA1*3.0d0*xy*cd/R5      + CA2*3.0d0*xd*U2/R5
return
end function

function duAdy_ds()
implicit none
double precision :: duAdy_ds(3)
duAdy_ds(1) =                              CA2*3.0d0*x*V2/R5
duAdy_ds(2) =  CA1*(s2d-3.0d0*ys/R2)/R3 + CA2*3.0d0*(y*V2+pq)/R5
duAdy_ds(3) = -CA1*(c2d-3.0d0*yt/R2)/R3 + CA2*3.0d0*d*V2/R5
return
end function

function duAdy_ts()
implicit none
double precision :: duAdy_ts(3)
duAdy_ts(1) = -CA1*3.0d0*xy/R5 - CA2*3.0d0*xq*W2/R5
duAdy_ts(2) =  CA1*(c2d-3.0d0*yt/R2)/R3 - CA2*3.0d0*(yq*W2+qq)/R5
duAdy_ts(3) =  CA1*(s2d-3.0d0*ys/R2)/R3 - CA2*3.0d0*dq*W2/R5
return
end function

function duAdy_vol()
implicit none
double precision :: duAdy_vol(3)
duAdy_vol(1) =  CA1*3.0d0*xy/R5
duAdy_vol(2) = -CA1*B3/R3
duAdy_vol(3) =  CA1*3.0d0*y*d/R5
return
end function

function duBdy_ss()
implicit none
double precision :: duBdy_ss(3)
duBdy_ss(1) = -3.0d0*xx*U2/R5                - CB*J2*sd
duBdy_ss(2) = -3.0d0*xy*U2/R5 - 3.0d0*x*q/R5 - CB*J4*sd
duBdy_ss(3) = -3.0d0*xc*U2/R5                - CB*K2*sd
return
end function

function duBdy_ds()
implicit none
double precision :: duBdy_ds(3)
duBdy_ds(1) = -3.0d0*x*V2/R5                  + CB*J1*cdsd
duBdy_ds(2) = -3.0d0*y*V2/R5   - 3.0d0*p*q/R5 + CB*J2*cdsd
duBdy_ds(3) = -3.0d0*c*V2/R5                  + CB*K1*cdsd
return
end function

function duBdy_ts()
implicit none
double precision :: duBdy_ts(3)
duBdy_ts(1) = 3.0d0*xq*W2/R5                - CB*J1*sdsd
duBdy_ts(2) = 3.0d0*yq*W2/R5 + 3.0d0*q*q/R5 - CB*J2*sdsd
duBdy_ts(3) = 3.0d0*c*q*W2/R5               - CB*K1*sdsd
return
end function

function duBdy_vol()
implicit none
double precision :: duBdy_vol(3)
duBdy_vol(1) = -CB*3.0d0*xy/R5
duBdy_vol(2) =  CB*B3/R3
duBdy_vol(3) = -CB*3.0d0*y*d/R5
return
end function

function duCdy_ss()
implicit none
double precision :: duCdy_ss(3)
duCdy_ss(1) =  CC*3.0d0*y*A5*cd/R5 + a*3.0d0*c*(A5*sd-5.0d0*yq*A7/R2)/R5
duCdy_ss(2) =  CC*3.0d0*x*B5*cd/R5 - a*15.0d0*xc*(2.0d0*y*sd+q*B7)/R7
duCdy_ss(3) = -CC*3.0d0*x*B5*sd/R5 + a*15.0d0*xc*(d*B7*sd-y*C7*cd)/R7
return
end function

function duCdy_ds()
implicit none
double precision :: duCdy_ds(3)
duCdy_ds(1) =  CC*3.0d0*x*(c2d-5.0d0*yt/R2)/R5 - a*15.0d0*xc*(s-7.0d0*y*pq/R2)/R7
duCdy_ds(2) =  CC*3.0d0*(2.0d0*y*c2d+t*B5)/R5 + a*3.0d0*c*(s2d-10.0d0*ys/R2-5.0d0*pq*B7/R2)/R5
duCdy_ds(3) =  CC*3.0d0*y*A5*cdsd/R5 - a*3.0d0*c*((3.0d0+A5)*c2d+35.0d0*y*d*pq/R4)/R5
return
end function

function duCdy_ts()
implicit none
double precision :: duCdy_ts(3)
duCdy_ts(1) = -CC*3.0d0*x*(s2d-5.0d0*ys/R2)/R5 - a*15.0d0*xc*(t-y+7.0d0*y*qq/R2)/R7 + &
              a*15.0d0*xy*z/R7
duCdy_ts(2) = -CC*3.0d0*(2.0d0*y*s2d+s*B5)/R5 - &
              a*3.0d0*c*(2.0d0*sdsd+10.0d0*y*(t-y)/R2-5.0d0*qq*B7/R2)/R5 - a*3.0d0*z*B5/R5
duCdy_ts(3) =  CC*3.0d0*y*(1.0d0-A5*sdsd)/R5 + &
              a*3.0d0*c*((3.0d0+A5)*s2d-5.0d0*y*d*(2.0d0-7.0d0*q*q/R2)/R2)/R5 - a*15.0d0*y*d*z/R7
return
end function

function duCdy_vol()
implicit none
double precision :: duCdy_vol(3)
duCdy_vol(1) = -CC*15.0d0*xy*d/R7
duCdy_vol(2) =  CC*3.0d0*d*B5/R5
duCdy_vol(3) = -CC*3.0d0*y*C5/R5
return
end function


function duAdz_ss()
implicit none
double precision :: duAdz_ss(3)
duAdz_ss(1) =  CA1*(cd+3.0d0*d*q/R2)/R3 + CA2*3.0d0*xx*U3/R5
duAdz_ss(2) =  CA1*3.0d0*xd*sd/R5       + CA2*3.0d0*xy*U3/R5
duAdz_ss(3) = -CA1*3.0d0*xd*cd/R5       + CA2*3.0d0*x*(d*U3-q)/R5
return
end function

function duAdz_ds()
implicit none
double precision :: duAdz_ds(3)
duAdz_ds(1) =                              CA2*3.0d0*x*V3/R5
duAdz_ds(2) =  CA1*(c2d+3.0d0*d*s/R2)/R3 + CA2*3.0d0*y*V3/R5
duAdz_ds(3) =  CA1*(s2d-3.0d0*d*t/R2)/R3 + CA2*3.0d0*(d*V3-pq)/R5
return
end function

function duAdz_ts()
implicit none
double precision :: duAdz_ts(3)
duAdz_ts(1) =  CA1*3.0d0*xd/R5 - CA2*3.0d0*xq*W3/R5
duAdz_ts(2) = -CA1*(s2d-3.0d0*d*t/R2)/R3 - CA2*3.0d0*yq*W3/R5
duAdz_ts(3) =  CA1*(c2d+3.0d0*d*s/R2)/R3 - CA2*3.0d0*q*(d*W3-q)/R5
return
end function

function duAdz_vol()
implicit none
double precision :: duAdz_vol(3)
duAdz_vol(1) = -CA1*3.0d0*xd/R5
duAdz_vol(2) = -CA1*3.0d0*y*d/R5
duAdz_vol(3) =  CA1*C3/R3
return
end function

function duBdz_ss()
implicit none
double precision :: duBdz_ss(3)
duBdz_ss(1) = -3.0d0*xx*U3/R5 + CB*K1*sd
duBdz_ss(2) = -3.0d0*xy*U3/R5 + CB*K2*sd
duBdz_ss(3) = -3.0d0*xc*U3/R5 + CB*3.d0*xy*sd/R5
return
end function

function duBdz_ds()
implicit none
double precision :: duBdz_ds(3)
duBdz_ds(1) = -3.0d0*x*V3/R5   - CB*K3*sd*cd
duBdz_ds(2) = -3.0d0*y*V3/R5   - CB*K1*sd*cd
duBdz_ds(3) = -3.0d0*c*V3/R5   + CB*A3*sd*cd/R3
return
end function

function duBdz_ts()
implicit none
double precision :: duBdz_ts(3)
duBdz_ts(1) = 3.0d0*xq*W3/R5 + CB*K3*sdsd
duBdz_ts(2) = 3.0d0*yq*W3/R5 + CB*K1*sdsd
duBdz_ts(3) = 3.0d0*c*q*W3/R5 - CB*A3*sdsd/R3
return
end function

function duBdz_vol()
implicit none
double precision :: duBdz_vol(3)
duBdz_vol(1) =  CB*3.0d0*xd/R5
duBdz_vol(2) =  CB*3.0d0*y*d/R5
duBdz_vol(3) = -CB*C3/R3
return
end function

function duCdz_ss()
implicit none
double precision :: duCdz_ss(3)
duCdz_ss(1) = -CC*3.0d0*d*A5*cd/R5  + a*3.0d0*c*(A5*cd+5.0d0*d*q*A7/R2)/R5
duCdz_ss(2) =  CC*15.0d0*xy*d*cd/R7 + a*15.0d0*xc*(d*B7*sd-y*C7*cd)/R7
duCdz_ss(3) = -CC*15.0d0*xy*d*sd/R7 + a*15.0d0*xc*(2.0d0*d*cd-q*C7)/R7
return
end function

function duCdz_ds()
implicit none
double precision :: duCdz_ds(3)
duCdz_ds(1) = -CC*3.0d0*x*(s2d-5.0d0*d*t/R2)/R5 - a*15.0d0*c*x*(t+7.0d0*d*pq/R2)/R7
duCdz_ds(2) = -CC*3.0d0*(d*B5*c2d+y*C5*s2d)/R5 - a*3.0d0*c*((3.0d0+A5)*c2d+35.0d0*y*d*pq/R4)/R5
duCdz_ds(3) = -CC*3.0d0*d*A5*sd*cd/R5 - a*3.0d0*c*(s2d-(10.0d0*d*t-5.0d0*pq*C7)/R2)/R5
return
end function

function duCdz_ts()
implicit none
double precision :: duCdz_ts(3)
duCdz_ts(1) = -CC*3.0d0*x*(c2d+5.0d0*d*s/R2)/R5 + a*15.0d0*xc*(s-d+7.0d0*d*qq/R2)/R7 - &
              a*3.0d0*x*(1.0d0+5.0d0*d*z/R2)/R5
duCdz_ts(2) =  CC*3.0d0*(d*B5*s2d-y*C5*c2d)/R5 + &
               a*3.0d0*c*((3.0d0+A5)*s2d-5.0d0*y*d*(2.0d0-7.0d0*qq/R2)/R2)/R5 - &
               a*3.0d0*y*(1.0d0+5.0d0*d*z/R2)/R5
duCdz_ts(3) = -CC*3.0d0*d*(1.0d0-A5*sdsd)/R5 - &
              a*3.0d0*c*(c2d+10.0d0*d*(s-d)/R2-5.0d0*qq*C7/R2)/R5 - a*3.0d0*z*(1.0d0+C5)/R5
return
end function

function duCdz_vol()
implicit none
double precision :: duCdz_vol(3)
duCdz_vol(1) = -CC*3.0d0*x*C5/R5
duCdz_vol(2) = -CC*3.0d0*y*C5/R5
duCdz_vol(3) =  CC*3.0d0*d*(2.0d0+C5)/R5
return
end function

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------- RECTANGULAR SOURCE ---------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


subroutine o92_rect_disp(disp,sta_coord,evdp,dip,slip,wid,len,lambda,shear_modulus)
!----
! Compute displacement vector at a point in an elastic half-space due to a rectangular dislocation
! Subroutine arguments:
!     disp(3): output displacement vector (m)
!     sta_coord(3): input station location(m) relative to center of rectangle;
!                   x=along-strike, y=hor-up-dip, z=depth
!     evdp: depth to center of rectangle (m), positive down
!     dip: source fault dip (degrees)
!     wid: down-dip width of fault (m)
!     len: along-strike length of fault (m)
!     slip(3): strike-slip, dip-slip, tensile-slip (m)
!     lambda, shear_modulus: half-space parameters (Pa)
!----

use io, only: stderr
use trig, only: pi

implicit none

! Arguments
double precision :: disp(3), sta_coord(3), evdp, dip, wid, len, slip(3), lambda, shear_modulus

! Local variables
integer :: i, j
double precision :: Css, Cds, Cts, u_ss(3), u_ds(3), u_ts(3), uA(3), uB(3), uC(3)

! write(0,*) 'o92_rect_disp: starting'
! write(0,*) 'o92_rect_disp: sta_coord: ',sta_coord
! write(0,*) 'o92_rect_disp: evdp     : ',evdp
! write(0,*) 'o92_rect_disp: dip      : ',dip
! write(0,*) 'o92_rect_disp: slip      : ',slip
! write(0,*) 'o92_rect_disp: wid      : ',wid
! write(0,*) 'o92_rect_disp: len      : ',len
! write(0,*) 'o92_rect_disp: lambda   : ',lambda
! write(0,*) 'o92_rect_disp: shear    : ',shear_modulus

! Initialize displacements
disp = 0.0d0

! If station is above surface, exit with warning
if (sta_coord(3).lt.0.0d0) then
    write(stderr,*) 'o92_rect_disp: station depth less than zero; setting displacement to zero'
    return
endif

! Calculate some parameters for fault geometry and elastic half-space
call dip_vars(dip)
! write(0,*) 'o92_rect_disp: cd       : ',cd
! write(0,*) 'o92_rect_disp: sd       : ',sd
call halfspace_vars(lambda,shear_modulus)

Css = slip(1)/(2.0d0*pi)
Cds = slip(2)/(2.0d0*pi)
Cts = slip(3)/(2.0d0*pi)

! Calculate the transformed coordinates (p, q, ksi, eta) for the terms that depend on the station
! at its actual depth below the surface
call rect_src_coords(sta_coord(1),sta_coord(2),-sta_coord(3),evdp,wid,len)

! Check for singular solution when observation point lies on fault edge
call check_singular_rect(isSingular,ksi_vec,eta_vec,q)
if (isSingular) then
    write(stderr,*) 'o92_rect_disp: solution is singular, setting displacements to zero'
    return
endif

! Calculate the components of displacement for z negative
u_ss = 0.0d0
u_ds = 0.0d0
u_ts = 0.0d0
do i = 1,2
    do j = 1,2
        call rect_src_vars(ksi_vec(i),eta_vec(j))
        if (dabs(Css).gt.1.0d-6) then
            uA = chinnery_factor(i,j)*uA_ss_rect()
            uB = chinnery_factor(i,j)*uB_ss_rect()
            uC = chinnery_factor(i,j)*z*uC_ss_rect()
            u_ss(1) = u_ss(1) + uA(1) + uB(1) + uC(1)
            u_ss(2) = u_ss(2) + (uA(2)+uB(2)+uC(2))*cd - (uA(3)+uB(3)+uC(3))*sd
            u_ss(3) = u_ss(3) + (uA(2)+uB(2)-uC(2))*sd + (uA(3)+uB(3)-uC(3))*cd
        endif
        if (dabs(Cds).gt.1.0d-6) then
            uA = chinnery_factor(i,j)*uA_ds_rect()
            uB = chinnery_factor(i,j)*uB_ds_rect()
            uC = chinnery_factor(i,j)*z*uC_ds_rect()
            u_ds(1) = u_ds(1) + uA(1) + uB(1) + uC(1)
            u_ds(2) = u_ds(2) + (uA(2)+uB(2)+uC(2))*cd - (uA(3)+uB(3)+uC(3))*sd
            u_ds(3) = u_ds(3) + (uA(2)+uB(2)-uC(2))*sd + (uA(3)+uB(3)-uC(3))*cd
        endif
        if (dabs(Cts).gt.1.0d-6) then
            uA = chinnery_factor(i,j)*uA_ts_rect()
            uB = chinnery_factor(i,j)*uB_ts_rect()
            uC = chinnery_factor(i,j)*z*uC_ts_rect()
            u_ts(1) = u_ts(1) + uA(1) + uB(1) + uC(1)
            u_ts(2) = u_ts(2) + (uA(2)+uB(2)+uC(2))*cd - (uA(3)+uB(3)+uC(3))*sd
            u_ts(3) = u_ts(3) + (uA(2)+uB(2)-uC(2))*sd + (uA(3)+uB(3)-uC(3))*cd
        endif
    enddo
enddo

! Calculate the transformed coordinates (p, q, ksi, eta) for the terms that depend on the station
! at the opposite depth
call rect_src_coords(sta_coord(1),sta_coord(2),sta_coord(3),evdp,wid,len)

! Check for singular solution when observation point lies on fault edge
call check_singular_rect(isSingular,ksi_vec,eta_vec,q)
if (isSingular) then
    write(stderr,*) 'o92_rect_disp: solution is singular, setting displacements to zero'
    return
endif

! Displacement components for z positive
do i = 1,2
    do j = 1,2
        call rect_src_vars(ksi_vec(i),eta_vec(j))
        if (dabs(Css).gt.1.0d-6) then
            uA = chinnery_factor(i,j)*uA_ss_rect()
            u_ss(1) = u_ss(1) - uA(1)
            u_ss(2) = u_ss(2) - uA(2)*cd + uA(3)*sd
            u_ss(3) = u_ss(3) - uA(2)*sd - uA(3)*cd
        endif
        if (dabs(Cds).gt.1.0d-6) then
            uA = chinnery_factor(i,j)*uA_ds_rect()
            u_ds(1) = u_ds(1) - uA(1)
            u_ds(2) = u_ds(2) - uA(2)*cd + uA(3)*sd
            u_ds(3) = u_ds(3) - uA(2)*sd - uA(3)*cd
        endif
        if (dabs(Cts).gt.1.0d-6) then
            uA = chinnery_factor(i,j)*uA_ts_rect()
            u_ts(1) = u_ts(1) - uA(1)
            u_ts(2) = u_ts(2) - uA(2)*cd + uA(3)*sd
            u_ts(3) = u_ts(3) - uA(2)*sd - uA(3)*cd
        endif
    enddo
enddo

disp(1) = Css*u_ss(1) + Cds*u_ds(1) + Cts*u_ts(1)
disp(2) = Css*u_ss(2) + Cds*u_ds(2) + Cts*u_ts(2)
disp(3) = Css*u_ss(3) + Cds*u_ds(3) + Cts*u_ts(3)

! write(0,*) 'o92_rect_disp: finished'

return
end subroutine o92_rect_disp

!--------------------------------------------------------------------------------------------------!

subroutine o92_rect_strain(strain,sta_coord,evdp,dip,slip,wid,len,lambda,shear_modulus)
!----
! Compute strain tensor at a point in an elastic half-space due to a rectangular dislocation
! Subroutine arguments:
!     strain(3,3): output strain tensor
!     sta_coord(3): input station location(m) relative to center of rectangle;
!                   x=along-strike, y=hor-up-dip, z=depth
!     evdp: depth to center of rectangle (m), positive down
!     dip: source fault dip (degrees)
!     wid: down-dip width of fault (m)
!     len: along-strike length of fault (m)
!     slip(3): strike-slip, dip-slip, tensile-slip (m)
!     lambda, shear_modulus: half-space parameters (Pa)
!----

use io, only: stderr

implicit none

! Arguments
double precision :: strain(3,3), sta_coord(3), evdp, dip, wid, len, slip(3), lambda, shear_modulus

! Local variables
double precision :: uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz

! Initialize strain tensor
strain = 0.0d0

! If station is above surface, exit with warning
if (sta_coord(3).lt.0.0d0) then
    write(stderr,*) 'o92_rect_strain: station depth less than zero; setting strain to zero'
    return
endif

! Partial derivatives of displacement equations
call o92_rect_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz, &
                       sta_coord,evdp,dip,slip,wid,len,lambda,shear_modulus)

! Strains
strain(1,1) = uxx
strain(1,2) = 1.0d0/2.0d0*(uxy+uyx)
strain(1,3) = 1.0d0/2.0d0*(uxz+uzx)
strain(2,1) = strain(1,2)
strain(2,2) = uyy
strain(2,3) = 1.0d0/2.0d0*(uyz+uzy)
strain(3,1) = strain(1,3)
strain(3,2) = strain(2,3)
strain(3,3) = uzz

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine o92_rect_partials(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz, &
                             sta_coord,evdp,dip,slip,wid,len,lambda,shear_modulus)
!----
! Compute partial derivatives of displacement at a point in an elastic half-space due to a
! rectangular source dislocation.
!----

use io, only: stderr
use trig, only: pi

implicit none

! Arguments
double precision :: uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz, &
                    sta_coord(3), evdp, dip, slip(3), wid, len, lambda, shear_modulus

! Local variables
double precision :: Css, Cds, Cts, dudx_ss(3), dudy_ss(3), dudz_ss(3), dudx_ds(3), dudy_ds(3), &
                    dudz_ds(3), dudx_ts(3), dudy_ts(3), dudz_ts(3), duA(3), duB(3), duC(3)
integer :: i, j

! Initialize partial derivatives
uxx = 0.0d0
uyx = 0.0d0
uzx = 0.0d0
uxy = 0.0d0
uyy = 0.0d0
uzy = 0.0d0
uxz = 0.0d0
uyz = 0.0d0
uzz = 0.0d0

! Calculate some parameters for fault geometry and elastic half-space
call dip_vars(dip)
call halfspace_vars(lambda,shear_modulus)

! Scaling factors
Css = slip(1)/(2.0d0*pi)
Cds = slip(2)/(2.0d0*pi)
Cts = slip(3)/(2.0d0*pi)

! Calculate the transformed coordinates (p, q, ksi, eta) for the terms that depend on the station
! at its actual depth below the surface
call rect_src_coords(sta_coord(1),sta_coord(2),-sta_coord(3),evdp,wid,len)

! Check for singular solution when observation point lies on fault edge
call check_singular_rect(isSingular,ksi_vec,eta_vec,q)
if (isSingular) then
    write(stderr,*) 'o92_rect_strain: solution is singular, setting strains to zero'
    return
endif

! Calculate the components of partial derivatives for z negative
dudx_ss = 0.0d0
dudy_ss = 0.0d0
dudz_ss = 0.0d0
dudx_ds = 0.0d0
dudy_ds = 0.0d0
dudz_ds = 0.0d0
dudx_ts = 0.0d0
dudy_ts = 0.0d0
dudz_ts = 0.0d0
do i = 1,2
    do j = 1,2
        call rect_src_vars(ksi_vec(i),eta_vec(j))

        if (dabs(Css).gt.1.0d-6) then
            ! X-derivatives
            duA = chinnery_factor(i,j)*duAdx_ss_rect()
            duB = chinnery_factor(i,j)*duBdx_ss_rect()
            duC = chinnery_factor(i,j)*z*duCdx_ss_rect()
            dudx_ss(1) = dudx_ss(1) + duA(1) + duB(1) + duC(1)
            dudx_ss(2) = dudx_ss(2) + (duA(2)+duB(2)+duC(2))*cd - (duA(3)+duB(3)+duC(3))*sd
            dudx_ss(3) = dudx_ss(3) + (duA(2)+duB(2)-duC(2))*sd + (duA(3)+duB(3)-duC(3))*cd

            ! Y-derivatives
            duA = chinnery_factor(i,j)*duAdy_ss_rect()
            duB = chinnery_factor(i,j)*duBdy_ss_rect()
            duC = chinnery_factor(i,j)*z*duCdy_ss_rect()
            dudy_ss(1) = dudy_ss(1) + duA(1) + duB(1) + duC(1)
            dudy_ss(2) = dudy_ss(2) + (duA(2)+duB(2)+duC(2))*cd - (duA(3)+duB(3)+duC(3))*sd
            dudy_ss(3) = dudy_ss(3) + (duA(2)+duB(2)-duC(2))*sd + (duA(3)+duB(3)-duC(3))*cd

            ! Z-derivatives
            duA = chinnery_factor(i,j)*duAdz_ss_rect()
            duB = chinnery_factor(i,j)*duBdz_ss_rect()
            duC = chinnery_factor(i,j)*(uC_ss_rect()+z*duCdz_ss_rect())
            dudz_ss(1) = dudz_ss(1) + duA(1) + duB(1) + duC(1)
            dudz_ss(2) = dudz_ss(2) + (duA(2)+duB(2)+duC(2))*cd - (duA(3)+duB(3)+duC(3))*sd
            dudz_ss(3) = dudz_ss(3) + (duA(2)+duB(2)-duC(2))*sd + (duA(3)+duB(3)-duC(3))*cd
        endif

        if (dabs(Cds).gt.1.0d-6) then
            ! X-derivatives
            duA = chinnery_factor(i,j)*duAdx_ds_rect()
            duB = chinnery_factor(i,j)*duBdx_ds_rect()
            duC = chinnery_factor(i,j)*z*duCdx_ds_rect()
            dudx_ds(1) = dudx_ds(1) + duA(1) + duB(1) + duC(1)
            dudx_ds(2) = dudx_ds(2) + (duA(2)+duB(2)+duC(2))*cd - (duA(3)+duB(3)+duC(3))*sd
            dudx_ds(3) = dudx_ds(3) + (duA(2)+duB(2)-duC(2))*sd + (duA(3)+duB(3)-duC(3))*cd

            ! Y-derivatives
            duA = chinnery_factor(i,j)*duAdy_ds_rect()
            duB = chinnery_factor(i,j)*duBdy_ds_rect()
            duC = chinnery_factor(i,j)*z*duCdy_ds_rect()
            dudy_ds(1) = dudy_ds(1) + duA(1) + duB(1) + duC(1)
            dudy_ds(2) = dudy_ds(2) + (duA(2)+duB(2)+duC(2))*cd - (duA(3)+duB(3)+duC(3))*sd
            dudy_ds(3) = dudy_ds(3) + (duA(2)+duB(2)-duC(2))*sd + (duA(3)+duB(3)-duC(3))*cd

            ! Z-derivatives
            duA = chinnery_factor(i,j)*duAdz_ds_rect()
            duB = chinnery_factor(i,j)*duBdz_ds_rect()
            duC = chinnery_factor(i,j)*(uC_ds_rect()+z*duCdz_ds_rect())
            dudz_ds(1) = dudz_ds(1) + duA(1) + duB(1) + duC(1)
            dudz_ds(2) = dudz_ds(2) + (duA(2)+duB(2)+duC(2))*cd - (duA(3)+duB(3)+duC(3))*sd
            dudz_ds(3) = dudz_ds(3) + (duA(2)+duB(2)-duC(2))*sd + (duA(3)+duB(3)-duC(3))*cd
        endif

        if (dabs(Cts).gt.1.0d-6) then
            ! X-derivatives
            duA = chinnery_factor(i,j)*duAdx_ts_rect()
            duB = chinnery_factor(i,j)*duBdx_ts_rect()
            duC = chinnery_factor(i,j)*z*duCdx_ts_rect()
            dudx_ts(1) = dudx_ts(1) + duA(1) + duB(1) + duC(1)
            dudx_ts(2) = dudx_ts(2) + (duA(2)+duB(2)+duC(2))*cd - (duA(3)+duB(3)+duC(3))*sd
            dudx_ts(3) = dudx_ts(3) + (duA(2)+duB(2)-duC(2))*sd + (duA(3)+duB(3)-duC(3))*cd

            ! Y-derivatives
            duA = chinnery_factor(i,j)*duAdy_ts_rect()
            duB = chinnery_factor(i,j)*duBdy_ts_rect()
            duC = chinnery_factor(i,j)*z*duCdy_ts_rect()
            dudy_ts(1) = dudy_ts(1) + duA(1) + duB(1) + duC(1)
            dudy_ts(2) = dudy_ts(2) + (duA(2)+duB(2)+duC(2))*cd - (duA(3)+duB(3)+duC(3))*sd
            dudy_ts(3) = dudy_ts(3) + (duA(2)+duB(2)-duC(2))*sd + (duA(3)+duB(3)-duC(3))*cd

            ! Z-derivatives
            duA = chinnery_factor(i,j)*duAdz_ts_rect()
            duB = chinnery_factor(i,j)*duBdz_ts_rect()
            duC = chinnery_factor(i,j)*(uC_ts_rect()+z*duCdz_ts_rect())
            dudz_ts(1) = dudz_ts(1) + duA(1) + duB(1) + duC(1)
            dudz_ts(2) = dudz_ts(2) + (duA(2)+duB(2)+duC(2))*cd - (duA(3)+duB(3)+duC(3))*sd
            dudz_ts(3) = dudz_ts(3) + (duA(2)+duB(2)-duC(2))*sd + (duA(3)+duB(3)-duC(3))*cd
        endif

    enddo
enddo

! Calculate the transformed coordinates (p, q, ksi, eta) for the terms that depend on the station
! at the opposite depth
call rect_src_coords(sta_coord(1),sta_coord(2),sta_coord(3),evdp,wid,len)

! Check for singular solution when observation point lies on fault edge
call check_singular_rect(isSingular,ksi_vec,eta_vec,q)
if (isSingular) then
    write(stderr,*) 'o92_rect_strain: solution is singular, setting strains to zero'
    return
endif

! Calculate the components of partial derivatives for z positive
do i = 1,2
    do j = 1,2
        call rect_src_vars(ksi_vec(i),eta_vec(j))

        if (dabs(Css).gt.1.0d-6) then
            ! X-derivatives
            duA = chinnery_factor(i,j)*duAdx_ss_rect()
            dudx_ss(1) = dudx_ss(1) - duA(1)
            dudx_ss(2) = dudx_ss(2) - duA(2)*cd + duA(3)*sd
            dudx_ss(3) = dudx_ss(3) - duA(2)*sd - duA(3)*cd

            ! Y-derivatives
            duA = chinnery_factor(i,j)*duAdy_ss_rect()
            dudy_ss(1) = dudy_ss(1) - duA(1)
            dudy_ss(2) = dudy_ss(2) - duA(2)*cd + duA(3)*sd
            dudy_ss(3) = dudy_ss(3) - duA(2)*sd - duA(3)*cd

            ! Z-derivatives
            duA = chinnery_factor(i,j)*duAdz_ss_rect()
            dudz_ss(1) = dudz_ss(1) + duA(1)
            dudz_ss(2) = dudz_ss(2) + duA(2)*cd - duA(3)*sd
            dudz_ss(3) = dudz_ss(3) + duA(2)*sd + duA(3)*cd
        endif

        if (dabs(Cds).gt.1.0d-6) then
            ! X-derivatives
            duA = chinnery_factor(i,j)*duAdx_ds_rect()
            dudx_ds(1) = dudx_ds(1) - duA(1)
            dudx_ds(2) = dudx_ds(2) - duA(2)*cd + duA(3)*sd
            dudx_ds(3) = dudx_ds(3) - duA(2)*sd - duA(3)*cd

            ! Y-derivatives
            duA = chinnery_factor(i,j)*duAdy_ds_rect()
            dudy_ds(1) = dudy_ds(1) - duA(1)
            dudy_ds(2) = dudy_ds(2) - duA(2)*cd + duA(3)*sd
            dudy_ds(3) = dudy_ds(3) - duA(2)*sd - duA(3)*cd

            ! Z-derivatives
            duA = chinnery_factor(i,j)*duAdz_ds_rect()
            dudz_ds(1) = dudz_ds(1) + duA(1)
            dudz_ds(2) = dudz_ds(2) + duA(2)*cd - duA(3)*sd
            dudz_ds(3) = dudz_ds(3) + duA(2)*sd + duA(3)*cd
        endif

        if (dabs(Cts).gt.1.0d-6) then
            ! X-derivatives
            duA = chinnery_factor(i,j)*duAdx_ts_rect()
            dudx_ts(1) = dudx_ts(1) - duA(1)
            dudx_ts(2) = dudx_ts(2) - duA(2)*cd + duA(3)*sd
            dudx_ts(3) = dudx_ts(3) - duA(2)*sd - duA(3)*cd

            ! Y-derivatives
            duA = chinnery_factor(i,j)*duAdy_ts_rect()
            dudy_ts(1) = dudy_ts(1) - duA(1)
            dudy_ts(2) = dudy_ts(2) - duA(2)*cd + duA(3)*sd
            dudy_ts(3) = dudy_ts(3) - duA(2)*sd - duA(3)*cd

            ! Z-derivatives
            duA = chinnery_factor(i,j)*duAdz_ts_rect()
            dudz_ts(1) = dudz_ts(1) + duA(1)
            dudz_ts(2) = dudz_ts(2) + duA(2)*cd - duA(3)*sd
            dudz_ts(3) = dudz_ts(3) + duA(2)*sd + duA(3)*cd
        endif

    enddo
enddo

uxx = Css*dudx_ss(1) + Cds*dudx_ds(1) + Cts*dudx_ts(1)
uxy = Css*dudy_ss(1) + Cds*dudy_ds(1) + Cts*dudy_ts(1)
uxz = Css*dudz_ss(1) + Cds*dudz_ds(1) + Cts*dudz_ts(1)
uyx = Css*dudx_ss(2) + Cds*dudx_ds(2) + Cts*dudx_ts(2)
uyy = Css*dudy_ss(2) + Cds*dudy_ds(2) + Cts*dudy_ts(2)
uyz = Css*dudz_ss(2) + Cds*dudz_ds(2) + Cts*dudz_ts(2)
uzx = Css*dudx_ss(3) + Cds*dudx_ds(3) + Cts*dudx_ts(3)
uzy = Css*dudy_ss(3) + Cds*dudy_ds(3) + Cts*dudy_ts(3)
uzz = Css*dudz_ss(3) + Cds*dudz_ds(3) + Cts*dudz_ts(3)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine rect_src_coords(xin,yin,zin,cin,wid,len)
!----
! Calculate rotated coordinates p, q, eta, and ksi from x, y, z (station location), c (source depth),
! fault width and length
!----

implicit none

! Arguments
double precision :: xin, yin, zin, cin, wid, len

! write(0,*) 'rect_src_coords: starting'

! Save passed arguments to module variables
x = xin
y = yin
z = zin
c = cin
d =  c - z

p = y*cd + d*sd
q = y*sd - d*cd

eta_vec(1) = p - 0.5d0*wid
eta_vec(2) = p + 0.5d0*wid
ksi_vec(1) = x - 0.5d0*len
ksi_vec(2) = x + 0.5d0*len

! write(0,*) 'rect_src_coords: finished'

return
end subroutine rect_src_coords

!--------------------------------------------------------------------------------------------------!

subroutine rect_src_vars(ksiin,etain)
!----
! Calculate variables from ksi, eta (receiver location relative to fault edge)
!----

implicit none

! Arguments
double precision :: ksiin, etain

! Local variables
double precision :: h, xxx
double precision, parameter :: eps = 1.0d-4

! Save module variables from arguments (other variables already set)
ksi = ksiin
eta = etain

! Some more rotated coordinates
ybar = eta*cd + q*sd
dbar = eta*sd - q*cd
cbar = dbar + z

R2 = ksi*ksi + eta*eta + q*q
R  = dsqrt(R2)
R3 = R2*R
R5 = R3*R2
Rd = R + dbar
h = q*cd - z

qq = q*q

! Singular case (i) from Okada (1992)
if (dabs(q).lt.eps) then
    TH = 0.0d0
else
    TH = datan(ksi*eta/(q*R))
endif

! Singular case (iii) from Okada (1992)
if (ksi.lt.0.0d0.and.dabs(R+ksi).lt.eps) then
    X11 = 0.0d0
    X32 = 0.0d0
    X53 = 0.0d0
    logRk = -dlog(R-ksi)
else
    Rk = R + ksi
    X11 = 1.0d0/(R*Rk)
    X32 = (2.0d0*R+ksi)/(R3*Rk*Rk)
    X53 = (8.0d0*R2+9.0d0*R*ksi+3.0d0*ksi*ksi)/(R5*Rk*Rk*Rk)
    logRk = dlog(R+ksi)
endif

! Singular case (iv) from Okada (1992)
if (eta.lt.0.0d0.and.dabs(R+eta).lt.eps) then
    Y11 = 0.0d0
    Y32 = 0.0d0
    Y53 = 0.0d0
    logRe = -dlog(R-eta)
else
    Re = R+eta
    Y11 = 1.0d0/(R*Re)
    Y32 = (2.0d0*R+eta)/(R3*Re*Re)
    Y53 = (8.0d0*R2+9.0d0*R*eta+3.0d0*eta*eta)/(R5*Re*Re*Re)
    logRe = dlog(R+eta)
endif
Y0  = Y11 - ksi*ksi*Y32

! Singular case (ii) from Okada (1992)
if (dabs(ksi).lt.eps) then
    I4 = 0.0d0
else
    xxx = sqrt(ksi*ksi+q*q)
    I4 = sd*ksi/(Rd*cd) + 2.0d0*datan((eta*(xxx+q*cd)+xxx*(R+xxx)*sd)/(ksi*(R+xxx)*cd))/(cdcd)
endif
I3 = ybar/(cd*Rd) - (logRe-sd*dlog(Rd))/(cdcd)
I1 = -ksi*cd/Rd - I4*sd
I2 = dlog(Rd) + I3*sd

Z32 = sd/R3 - h*Y32
Z53 = 3.0d0*sd/R5 - h*Y53
Z0  = Z32 - ksi*ksi*Z53

D11 = 1.0d0/(R*Rd)

K1 = ksi*(D11-Y11*sd)/cd
K3 = (q*Y11-ybar*D11)/cd
K2 = 1.0d0/R+K3*sd
K4 = ksi*Y11*cd-K1*sd

J2 = ksi*ybar*D11/Rd
J5 = -(dbar+ybar*ybar/Rd)*D11
J3 = (K1-J2*sd)/cd
J6 = (K3-J5*sd)/cd
J1 = J5*cd - J6*sd
J4 = -ksi*Y11 - J2*cd + J3*sd

E2 = sd/R - ybar*q/R3
F2 = dbar/R3 + ksi*ksi*Y32*sd
G2 = 2.0d0*X11*sd - ybar*q*X32
H2 = dbar*q*X32 + ksi*q*Y32*sd
P2 = cd/R3 + q*Y32*sd
Q2 = 3.0d0*cbar*dbar/R5 - (z*Y32+Z32+Z0)*sd

E3 = cd/R + dbar*q/R3
F3 = ybar/R3 + ksi*ksi*Y32*cd
G3 = 2.0d0*X11*cd + dbar*q*X32
H3 = ybar*q*X32 + ksi*q*Y32*cd
P3 = sd/R3 - q*Y32*cd
Q3 = 3.0d0*cbar*ybar/R5 + q*Y32 - (z*Y32+Z32+Z0)*cd

return
end subroutine rect_src_vars

!--------------------------------------------------------------------------------------------------!
! COMPONENTS OF DISPLACEMENT FOR A RECTANGULAR SOURCE DISLOCATION - EDIT AT YOUR OWN RISK!         !
!--------------------------------------------------------------------------------------------------!

function uA_ss_rect()
implicit none
double precision :: uA_ss_rect(3)
uA_ss_rect(1) = 0.5d0*TH + CA2*ksi*q*Y11
uA_ss_rect(2) =            CA2*q/R
uA_ss_rect(3) = CA1*logRe   - CA2*qq*Y11
return
end function

function uA_ds_rect()
implicit none
double precision :: uA_ds_rect(3)
uA_ds_rect(1) =             CA2*q/R
uA_ds_rect(2) = TH*0.5d0  + CA2*eta*q*X11
uA_ds_rect(3) = CA1*logRk - CA2*qq*X11
return
end function

function uA_ts_rect()
implicit none
double precision :: uA_ts_rect(3)
uA_ts_rect(1) = -CA1*logRe - CA2*qq*Y11
uA_ts_rect(2) = -CA1*logRk - CA2*qq*X11
uA_ts_rect(3) =   TH*0.5d0 - CA2*q*(eta*X11+ksi*Y11)
return
end function

function uB_ss_rect()
implicit none
double precision :: uB_ss_rect(3)
uB_ss_rect(1) = -ksi*q*Y11 - TH - CB*I1*sd
uB_ss_rect(2) = -q/R            + CB*ybar*sd/Rd
uB_ss_rect(3) =  qq*Y11         - CB*I2*sd
return
end function

function uB_ds_rect()
implicit none
double precision :: uB_ds_rect(3)
uB_ds_rect(1) = -q/R            + CB*I3*cdsd
uB_ds_rect(2) = -eta*q*X11 - TH - CB*ksi*cdsd/Rd
uB_ds_rect(3) =  qq*X11         + CB*I4*cdsd
return
end function

function uB_ts_rect()
implicit none
double precision :: uB_ts_rect(3)
uB_ts_rect(1) = qq*Y11                   - CB*I3*sdsd
uB_ts_rect(2) = qq*X11                   + CB*ksi*sdsd/Rd
uB_ts_rect(3) = q*(eta*X11+ksi*Y11) - TH - CB*I4*sdsd
return
end function

function uC_ss_rect()
implicit none
double precision :: uC_ss_rect(3)
uC_ss_rect(1) = CC*ksi*Y11*cd - a*ksi*q*Z32
uC_ss_rect(2) = CC*(cd/R+2.0d0*q*Y11*sd) - a*cbar*q/R3
uC_ss_rect(3) = CC*q*Y11*cd - a*(cbar*eta/R3-z*Y11+ksi*ksi*Z32)
return
end function

function uC_ds_rect()
implicit none
double precision :: uC_ds_rect(3)
uC_ds_rect(1) = CC*cd/R - q*Y11*sd - a*cbar*q/R3
uC_ds_rect(2) = CC*ybar*X11 - a*cbar*eta*q*X32
uC_ds_rect(3) = -dbar*X11 - ksi*Y11*sd - a*cbar*(X11 - qq*X32)
return
end function

function uC_ts_rect()
implicit none
double precision :: uC_ts_rect(3)
uC_ts_rect(1) = -CC*(sd/R+q*Y11*cd) - a*(z*Y11-qq*Z32)
uC_ts_rect(2) =  CC*2.0d0*ksi*Y11*sd + dbar*X11 - a*cbar*(X11-qq*X32)
uC_ts_rect(3) =  CC*(ybar*X11+ksi*Y11*cd) + a*q*(cbar*eta*X32+ksi*Z32)
return
end function

!--------------------------------------------------------------------------------------------------!
! COMPONENTS OF PARTIAL DERIVATIVES FOR A RECTANGULAR SOURCE DISLOCATION - EDIT AT YOUR OWN RISK!  !
!--------------------------------------------------------------------------------------------------!

function duAdx_ss_rect()
implicit none
double precision :: duAdx_ss_rect(3)
duAdx_ss_rect(1) = -CA1*q*Y11   - CA2*ksi*ksi*q*Y32
duAdx_ss_rect(2) =              - CA2*ksi*q/R3
duAdx_ss_rect(3) =  CA1*ksi*Y11 + CA2*ksi*qq*Y32
return
end function

function duAdx_ds_rect()
implicit none
double precision :: duAdx_ds_rect(3)
duAdx_ds_rect(1) =              - CA2*ksi*q/R3
duAdx_ds_rect(2) = -q*Y11*0.5d0 - CA2*eta*q/R3
duAdx_ds_rect(3) =  CA1/R       + CA2*qq/R3
return
end function

function duAdx_ts_rect()
implicit none
double precision :: duAdx_ts_rect(3)
duAdx_ts_rect(1) = -CA1*ksi*Y11 + CA2*ksi*qq*Y32
duAdx_ts_rect(2) = -CA1/R       + CA2*qq/R3
duAdx_ts_rect(3) = -CA1*q*Y11   - CA2*qq*q*Y32
return
end function

function duBdx_ss_rect()
implicit none
double precision :: duBdx_ss_rect(3)
duBdx_ss_rect(1) =  ksi*ksi*q*Y32     - CB*J1*sd
duBdx_ss_rect(2) =  ksi*q/R3          - CB*J2*sd
duBdx_ss_rect(3) = -ksi*qq*Y32       - CB*J3*sd
return
end function

function duBdx_ds_rect()
implicit none
double precision :: duBdx_ds_rect(3)
duBdx_ds_rect(1) =  ksi*q/R3          + CB*J4*cdsd
duBdx_ds_rect(2) =  eta*q/R3 + q*Y11  + CB*J5*cdsd
duBdx_ds_rect(3) = -qq/R3            + CB*J6*cdsd
return
end function

function duBdx_ts_rect()
implicit none
double precision :: duBdx_ts_rect(3)
duBdx_ts_rect(1) = -ksi*qq*Y32      - CB*J4*sdsd
duBdx_ts_rect(2) = -qq/R3           - CB*J5*sdsd
duBdx_ts_rect(3) = qq*q*Y32         - CB*J6*sdsd
return
end function

function duCdx_ss_rect()
implicit none
double precision :: duCdx_ss_rect(3)
duCdx_ss_rect(1) =  CC*Y0*cd - a*q*Z0
duCdx_ss_rect(2) = -CC*ksi*(cd/R3+2.0d0*q*Y32*sd) + a*3.0d0*cbar*ksi*q/R5
duCdx_ss_rect(3) = -CC*ksi*q*Y32*cd + a*ksi*(3.0d0*cbar*eta/R5-z*Y32-Z32-Z0)
return
end function

function duCdx_ds_rect()
implicit none
double precision :: duCdx_ds_rect(3)
duCdx_ds_rect(1) = -CC*ksi*cd/R3 + ksi*q*Y32*sd + a*3.0d0*cbar*ksi*q/R5
duCdx_ds_rect(2) = -CC*ybar/R3 + a*3.0d0*cbar*eta*q/R5
duCdx_ds_rect(3) = dbar/R3 - Y0*sd + a*cbar*(1.0d0-3.0d0*qq/R2)/R3
return
end function

function duCdx_ts_rect()
implicit none
double precision :: duCdx_ts_rect(3)
duCdx_ts_rect(1) = CC*ksi*sd/R3 + ksi*q*Y32*cd + a*ksi*(3.0d0*cbar*eta/R5 - 2.0d0*Z32 - Z0)
duCdx_ts_rect(2) =  CC*2.0d0*Y0*sd - dbar/R3 + a*cbar*(1.0d0-3.0d0*qq/R2)/R3
duCdx_ts_rect(3) = -CC*(ybar/R3-Y0*cd) - a*(3.0d0*cbar*eta*q/R5-q*Z0)
return
end function

function duAdy_ss_rect()
implicit none
double precision :: duAdy_ss_rect(3)
duAdy_ss_rect(1) =  CA1*ksi*Y11*sd + dbar*X11*0.5d0 + CA2*ksi*F2
duAdy_ss_rect(2) =                                    CA2*E2
duAdy_ss_rect(3) =  CA1*(cd/R+q*Y11*sd)             - CA2*q*F2
return
end function

function duAdy_ds_rect()
implicit none
double precision :: duAdy_ds_rect(3)
duAdy_ds_rect(1) =                                     CA2*E2
duAdy_ds_rect(2) =   CA1*dbar*X11 + ksi*Y11*sd*0.5d0 + CA2*eta*G2
duAdy_ds_rect(3) =   CA1*ybar*X11                    - CA2*q*G2
return
end function

function duAdy_ts_rect()
implicit none
double precision :: duAdy_ts_rect(3)
duAdy_ts_rect(1) =  -CA1*(cd/R+q*Y11*sd)             - CA2*q*F2
duAdy_ts_rect(2) =  -CA1*ybar*X11                    - CA2*q*G2
duAdy_ts_rect(3) =   CA1*(dbar*X11+ksi*Y11*sd)       + CA2*q*H2
return
end function

function duBdy_ss_rect()
implicit none
double precision :: duBdy_ss_rect(3)
duBdy_ss_rect(1) =  -ksi*F2 - dbar*X11     + CB*(ksi*Y11+J4)*sd
duBdy_ss_rect(2) =  -E2                    + CB*(1.0d0/R+J5)*sd
duBdy_ss_rect(3) =   q*F2                  - CB*(q*Y11-J6)*sd
return
end function

function duBdy_ds_rect()
implicit none
double precision :: duBdy_ds_rect(3)
duBdy_ds_rect(1) =  -E2                    + CB*J1*cdsd
duBdy_ds_rect(2) =  -eta*G2 - ksi*Y11*sd + CB*J2*cdsd
duBdy_ds_rect(3) =   q*G2                  + CB*J3*cdsd
return
end function

function duBdy_ts_rect()
implicit none
double precision :: duBdy_ts_rect(3)
duBdy_ts_rect(1) =  q*F2 - CB*J1*sdsd
duBdy_ts_rect(2) =  q*G2 - CB*J2*sdsd
duBdy_ts_rect(3) = -q*H2 - CB*J3*sdsd
return
end function

function duCdy_ss_rect()
implicit none
double precision :: duCdy_ss_rect(3)
duCdy_ss_rect(1) = -CC*ksi*P2*cd - a*ksi*Q2
duCdy_ss_rect(2) = 2.0d0*CC*(dbar/R3-Y0*sd)*sd - ybar*cd/R3 - &
                   a*((cbar+dbar)*sd/R3-eta/R3-3.0d0*cbar*ybar*q/R5)
duCdy_ss_rect(3) = -CC*q/R3 + (ybar/R3-Y0*cd)*sd + &
                   a*((cbar+dbar)*cd/R3+3.0d0*cbar*dbar*q/R5-(Y0*cd+q*Z0)*sd)
return
end function

function duCdy_ds_rect()
implicit none
double precision :: duCdy_ds_rect(3)
duCdy_ds_rect(1) = -CC*eta/R3 + Y0*sdsd - a*((cbar+dbar)*sd/R3-3.0d0*cbar*ybar*q/R5)
duCdy_ds_rect(2) =  CC*(X11-ybar*ybar*X32) - a*cbar*((dbar+2.0d0*q*cd)*X32-ybar*eta*q*X53)
duCdy_ds_rect(3) =  ksi*P2*sd + ybar*dbar*X32 + a*cbar*((ybar+2.0d0*q*sd)*X32-ybar*qq*X53)
return
end function

function duCdy_ts_rect()
implicit none
double precision :: duCdy_ts_rect(3)
duCdy_ts_rect(1) =  CC*(q/R3+Y0*cdsd) + a*(z*cd/R3+3.0d0*cbar*dbar*q/R5-q*Z0*sd)
duCdy_ts_rect(2) = -CC*2.0d0*ksi*P2*sd - ybar*dbar*X32 + a*cbar*((ybar+2.0d0*q*sd)*X32-ybar*qq*X53)
duCdy_ts_rect(3) = -CC*(ksi*P2*cd-X11+ybar*ybar*X32) + &
                   a*cbar*((dbar+2.0d0*q*cd)*X32-ybar*eta*q*X53) + a*ksi*Q2
return
end function

function duAdz_ss_rect()
implicit none
double precision :: duAdz_ss_rect(3)
duAdz_ss_rect(1) =  CA1*ksi*Y11*cd + ybar*X11*0.5d0 + CA2*ksi*F3
duAdz_ss_rect(2) =                                    CA2*E3
duAdz_ss_rect(3) = -CA1*(sd/R-q*Y11*cd)             - CA2*q*F3
return
end function

function duAdz_ds_rect()
implicit none
double precision :: duAdz_ds_rect(3)
duAdz_ds_rect(1) =                                    CA2*E3
duAdz_ds_rect(2) =  CA1*ybar*X11 + ksi*Y11*cd*0.5d0 + CA2*eta*G3
duAdz_ds_rect(3) = -CA1*dbar*X11                    - CA2*q*G3
return
end function

function duAdz_ts_rect()
implicit none
double precision :: duAdz_ts_rect(3)
duAdz_ts_rect(1) =   CA1*(sd/R-q*Y11*cd)             - CA2*q*F3
duAdz_ts_rect(2) =  CA1*dbar*X11                    - CA2*q*G3
duAdz_ts_rect(3) =   CA1*(ybar*X11+ksi*Y11*cd)       - CA2*q*H3
return
end function

function duBdz_ss_rect()
implicit none
double precision :: duBdz_ss_rect(3)
duBdz_ss_rect(1) = -ksi*F3 - ybar*X11     + CB*K1*sd
duBdz_ss_rect(2) = -E3                    + CB*ybar*D11*sd
duBdz_ss_rect(3) =  q*F3                  + CB*K2*sd
return
end function

function duBdz_ds_rect()
implicit none
double precision :: duBdz_ds_rect(3)
duBdz_ds_rect(1) = -E3                    - CB*K3*cdsd
duBdz_ds_rect(2) = -eta*G3 - ksi*Y11*cd - CB*ksi*D11*cdsd
duBdz_ds_rect(3) =  q*G3                  - CB*K4*cdsd
return
end function

function duBdz_ts_rect()
implicit none
double precision :: duBdz_ts_rect(3)
duBdz_ts_rect(1) =  q*F3 + CB*K3*sdsd
duBdz_ts_rect(2) =  q*G3 + CB*ksi*D11*sdsd
duBdz_ts_rect(3) = -q*H3 + CB*K4*sdsd
return
end function

function duCdz_ss_rect()
implicit none
double precision :: duCdz_ss_rect(3)
duCdz_ss_rect(1) = CC*ksi*P3*cd - a*ksi*Q3
duCdz_ss_rect(2) = 2.0d0*CC*(ybar/R3-Y0*cd)*sd + dbar*cd/R3 - &
                   a*((cbar+dbar)*cd/R3+3.0d0*cbar*dbar*q/R5)
duCdz_ss_rect(3) = (ybar/R3-Y0*cd)*cd - &
                   a*((cbar+dbar)*sd/R3-3.0d0*cbar*ybar*q/R5-Y0*sd*sd+q*Z0*cd)
return
end function

function duCdz_ds_rect()
implicit none
double precision :: duCdz_ds_rect(3)
duCdz_ds_rect(1) = -q/R3 + Y0*sd*cd - a*((cbar+dbar)*cd/R3+3.0d0*cbar*dbar*q/R5)
duCdz_ds_rect(2) = CC*ybar*dbar*X32 - a*cbar*((ybar-2.0d0*q*sd)*X32+dbar*eta*q*X53)
duCdz_ds_rect(3) = -ksi*P3*sd + X11 - dbar*dbar*X32 - a*cbar*((dbar-2.0*q*cd)*X32-dbar*qq*X53)
return
end function

function duCdz_ts_rect()
implicit none
double precision :: duCdz_ts_rect(3)
duCdz_ts_rect(1) = -eta/R3 + Y0*cdcd - a*(z*sd/R3-3.0d0*cbar*ybar*q/R5-Y0*sdsd+q*Z0*cd)
duCdz_ts_rect(2) =  CC*2.0d0*ksi*P3*sd - X11 + dbar*dbar*X32 - &
                    a*cbar*((dbar-2.0d0*q*cd)*X32-dbar*qq*X53)
duCdz_ts_rect(3) =  CC*(ksi*P3*cd+ybar*dbar*X32) + a*cbar*((ybar-2.0d0*q*sd)*X32+dbar*eta*q*X53) + &
                    a*ksi*Q3
return
end function

!--------------------------------------------------------------------------------------------------!

subroutine check_singular_pt(isPointOnSource,Rin)
!----
! Check to see if observation point lies on the source point
!----

implicit none
! Arguments
logical :: isPointOnSource
double precision :: Rin

! write(0,*) 'check_singular_pt: starting'

isPointOnSource = .false.

 ! Check if point lies on fault edge (q=0 and either ksi=0 or eta=0)
if (dabs(Rin).lt.1.0d-6) then
    isPointOnSource = .false.
endif

! write(0,*) 'check_singular_pt: finished'

return
end subroutine check_singular_pt

!--------------------------------------------------------------------------------------------------!

subroutine check_singular_rect(isPointOnFaultEdge,ksiin,etain,qin)
!----
! Check to see if observation point lies on one of singular locations
!----

implicit none

! Arguments
logical :: isPointOnFaultEdge
double precision :: ksiin(2),etain(2),qin

! Local variables
double precision :: k1k2, e1e2

! write(0,*) 'check_singular_rect: starting'

isPointOnFaultEdge = .false.

k1k2 = ksiin(1)*ksiin(2)
e1e2 = etain(1)*etain(2)

 ! Check if point lies on fault edge (q=0 and either ksi=0 or eta=0)
if (dabs(qin).lt.1.0d-6) then

    ! Point lies on edge with one of the etas=0
    if (k1k2.le.0.0d0.and.dabs(e1e2).lt.1.0d-6) then
        isPointOnFaultEdge = .true.
        return

    ! Point lies on edge with one of the ksis=0
    elseif (e1e2.le.0.0d0.and.dabs(k1k2).lt.1.0d-6) then
        isPointOnFaultEdge = .true.
        return

    endif

endif

! write(0,*) 'check_singular_rect: finished'

return
end subroutine check_singular_rect

end module
