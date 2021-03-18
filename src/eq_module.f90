!--------------------------------------------------------------------------------------------------!
! Module eq
!
! Routines for manipulating earthquake source parameters, including moment tensors, principal
! axes, fault geometries, slip kinematics, magnitudes, and more.
!
! The algorithms for most of these subroutines come from two sources:
!     Aki, K., Richards, P.G. (2002). Quantitative Seismology, Sausalito, CA: University Science
!         Books.
!     Jost, M.L., Herrmann, R.B. (1989). A student's guide to and review of moment tensors.
!         Seismological Research Letters 60, 37-57.
!
!--------------------------------------------------------------------------------------------------!
module eq

double precision, parameter :: dynecm2nm = 1.0d-7
double precision, parameter :: nm2dynecm = 1.0d7

public :: dynecm2nm
public :: nm2dynecm

public :: mag2mom
public :: mom2mag

public :: mij2dcp
public :: mij2mag
public :: mij2mom
public :: mij2pnt
public :: mij2sdr
public :: mij2sv
public :: mij2ter

public :: pnt2dcp
public :: pnt2mag
public :: pnt2mij
public :: pnt2mom
public :: pnt2mom_nondc
public :: pnt2sdr
public :: pnt2sv
public :: pnt2ter

public :: sdr2mij
public :: sdr2pnt
public :: sdr2sdr
public :: sdr2sv
public :: sdr2ter

public :: sdr2kagan
public :: mij2kagan
public :: kaganangle

public :: empirical

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!---------------------------------- mij2xxx -------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine mij2dcp(mrr,mpp,mtt,mrp,mrt,mtp,dcp)
!----
! Compute the percent double couple of a moment tensor.
!----

implicit none

! Arguments
double precision :: mrr, mpp, mtt, mrp, mrt, mtp, dcp

! Local variables
double precision :: pnt(12)

! Compute eigenvectors and eigenvalues
call mij2pnt(mrr,mpp,mtt,mrp,mrt,mtp,pnt)

! Compute double couple percentage
call pnt2dcp(pnt,dcp)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine mij2mag(mrr,mpp,mtt,mrp,mrt,mtp,mag,mode)
!----
! Calculate the moment magnitude of the moment tensor
!----

use io, only: stderr
implicit none

! Arguments
double precision :: mrr, mpp, mtt, mrp, mrt, mtp, mag
character(len=*) :: mode

! Local variables
double precision :: pnt(12), mom

call mij2pnt(mrr,mpp,mtt,mrp,mrt,mtp,pnt)
if (mode.eq.'dc') then
    call pnt2mom(pnt,mom)
elseif (mode.eq.'nondc') then
    call pnt2mom_nondc(pnt,mom)
else
    write(stderr,*) 'mij2mag: no mode named "',trim(mode),'"'
    write(stderr,*) 'returning with magnitude of -100'
    mag = -100.0d0
    return
endif
call mom2mag(mom,mag)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine mij2mom(mrr,mpp,mtt,mrp,mrt,mtp,mom,mode)
!----
! Calculate the scalar moment of the moment tensor
!----

use io, only: stderr
implicit none

! Arguments
double precision :: mrr, mpp, mtt, mrp, mrt, mtp, mom
character(len=*) :: mode

! Local variables
double precision :: pnt(12)

call mij2pnt(mrr,mpp,mtt,mrp,mrt,mtp,pnt)
if (mode.eq.'dc') then
    call pnt2mom(pnt,mom)
elseif (mode.eq.'nondc') then
    call pnt2mom_nondc(pnt,mom)
else
    write(stderr,*) 'mij2mag: no mode named "',trim(mode),'"'
    write(stderr,*) 'returning with moment of 0'
    mom = 0.0d0
    return
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine mij2pnt(mrr,mtt,mpp,mrt,mrp,mtp,pnt)
!----
! Compute the P, N, and T axis orientations (eigenvectors) in east-north-vertical coordinates and
! magnitudes (eigenvalues) for a moment tensor given in spherical coordinates.
!----

use io, only: stderr

#ifdef USE_LAPACK
! With LAPACK, we only need the cross_product and dot_product routines from the algebra module
use algebra, only: cross_product, dot_product
#else
! Without LAPACK, we also need the eigen routines from the algebra module
use algebra, only: cross_product, dot_product, jacobi, eig_sort
#endif


implicit none

! Arguments
double precision :: mrr, mtt, mpp, mrt, mrp, mtp, pnt(12)

! Local variables
double precision :: m_enz(3,3), eigvals(3), xp(3), dp


#ifdef USE_LAPACK

! Variables for running LAPACK dysev (double precision symmetric matrix eigen analysis)
integer :: info, lwork
double precision, allocatable :: work(:)

#else

! Variables for running Numerical Recipes jacobi (symmetric matrix eigen analysis) and eig_sort
integer :: nrot, ierr
double precision :: eigvecs(3,3)

#endif


! Initialize PNT
pnt = 0.0d0


! Convert moment tensor from spherical coordinates to local east-north-vertical coordinates, where
! r=up, t=south, and p=east
m_enz(3,3) =  mrr     ! mrr =  mzz
m_enz(2,2) =  mtt     ! mtt =  mnn
m_enz(1,1) =  mpp     ! mpp =  mee
m_enz(2,3) = -mrt     ! mrt = -mnz
m_enz(3,2) = -mrt
m_enz(1,3) =  mrp     ! mrp =  mez
m_enz(3,1) =  mrp
m_enz(1,2) = -mtp     ! mtp = -men
m_enz(2,1) = -mtp


! Use LAPACK eigensystem solver if it is available. Otherwise, use Numerical Recipes routine.

#ifdef USE_LAPACK

! Query the optimal workspace
allocate(work(1))
lwork = -1
call dsyev('Vectors','Upper',3,m_enz,3,eigvals,work,lwork,info)
lwork = int(work(1))
deallocate(work)
allocate(work(lwork))

! Solve eigenproblem
call dsyev('Vectors','Upper',3,m_enz,3,eigvals,work,lwork,info)

if (info.ne.0) then
    write(stderr,*) 'mij2pnt: error in dsyev; returning pnt with zeros'
    return
endif

! Eigenvectors overwrite the input array
pnt(1) = m_enz(1,1)
pnt(2) = m_enz(2,1)
pnt(3) = m_enz(3,1)
pnt(4) = m_enz(1,2)
pnt(5) = m_enz(2,2)
pnt(6) = m_enz(3,2)
pnt(7) = m_enz(1,3)
pnt(8) = m_enz(2,3)
pnt(9) = m_enz(3,3)

! Eigenvalues are in the w array
pnt(10) = eigvals(1)
pnt(11) = eigvals(2)
pnt(12) = eigvals(3)

#else

! Calculate eigenvalues and eigenvectors with Numerical Recipes routine
call jacobi(m_enz,3,3,eigvals,eigvecs,nrot,ierr)
if (ierr.ne.0) then
    write(stderr,*) 'mij2pnt: error in jacobi; returning pnt with zeros'
    return
endif

call eig_sort(eigvals,eigvecs,3,3)

! Eigenvectors are in the mat array
pnt(1) = eigvecs(1,1)
pnt(2) = eigvecs(2,1)
pnt(3) = eigvecs(3,1)
pnt(4) = eigvecs(1,2)
pnt(5) = eigvecs(2,2)
pnt(6) = eigvecs(3,2)
pnt(7) = eigvecs(1,3)
pnt(8) = eigvecs(2,3)
pnt(9) = eigvecs(3,3)

! Eigenvalues are in the vec array
pnt(10) = eigvals(1)
pnt(11) = eigvals(2)
pnt(12) = eigvals(3)

#endif

! Verify that eigenvectors are right handed
call cross_product(pnt(1:3),pnt(4:6),xp)
call dot_product(xp,pnt(7:9),dp)
if (dp.lt.0.0d0) then
    pnt(7) = -pnt(7)
    pnt(8) = -pnt(8)
    pnt(9) = -pnt(9)
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine mij2sdr(mrr,mtt,mpp,mrt,mrp,mtp,str1,dip1,rak1,str2,dip2,rak2)
!----
! Compute the strike, dip, and rake angles of both nodal planes of the best-fitting double couple
! source for a moment tensor source provided in spherical coordinates used by the Global Centroid
! Moment Tensor project.
!----

implicit none

! Arguments
double precision :: mrr, mtt, mpp, mrt, mrp, mtp, str1, dip1, rak1, str2, dip2, rak2

! Local variables
double precision :: pnt(12)

! Compute the P, N, and T axes for this moment tensor
call mij2pnt(mrr,mtt,mpp,mrt,mrp,mtp,pnt)

! Compute the strike, dip, and rake of one nodal plane
call pnt2sdr(pnt,str1,dip1,rak1,str2,dip2,rak2)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine mij2sv(mrr,mtt,mpp,mrt,mrp,mtp,slip_vec1,slip_vec2)
!----
! Compute the slip vectors of the two nodal planes for the best-fitting double couple source
! of a moment tensor.
!----

implicit none

! Arguments
double precision :: mrr, mtt, mpp, mrt, mrp, mtp, slip_vec1(3), slip_vec2(3)

! Local variables
double precision :: str1, dip1, rak1, str2, dip2, rak2

! Calculate best fitting double couple
call mij2sdr(mrr,mtt,mpp,mrt,mrp,mtp,str1,dip1,rak1,str2,dip2,rak2)

! Compute slip vector for both nodal planes
call sdr2sv(str1,dip1,rak1,slip_vec1)
call sdr2sv(str2,dip2,rak2,slip_vec2)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine mij2ter(mrr,mtt,mpp,mrt,mrp,mtp,fth,fss,fno)
!----
! Compute the fraction thrust, strike-slip, and normal for a moment tensor.
!----

implicit none

double precision :: mrr, mtt, mpp, mrt, mrp, mtp, fth, fss, fno
double precision :: pnt(12)

! Find eigenvectors and eigenvalues
call mij2pnt(mrr,mtt,mpp,mrt,mrp,mtp,pnt)

! Compute fraction thrust, strike-slip, and normal
call pnt2ter(pnt,fth,fss,fno)

return
end subroutine

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!---------------------------------- pnt2xxx -------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine pnt2dcp(pnt,dcp)
!----
! Compute the percent double couple from the P, N, and T values.
!----

implicit none

! Arguments
double precision :: pnt(12), dcp

if (abs(pnt(10)).gt.abs(pnt(12))) then
    dcp = 1.0d0 - 2.0d0*abs(pnt(11)/pnt(10))
else
    dcp = 1.0d0 - 2.0d0*abs(pnt(11)/pnt(12))
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine pnt2mij(pnt,mrr,mtt,mpp,mrt,mrp,mtp)
!----
! Transform the eigenvectors and eigenvalues of the moment tensor (in east-north-vertical
! coordinates) to the full moment tensor (in Global Centroid Moment Tensor project spherical
! coordinates).
!----

implicit none

! Arguments
double precision :: pnt(12), mrr, mtt, mpp, mrt, mrp, mtp

! Local variables
double precision :: mxx, myy, mzz, mxy, mxz, myz

! Calculate the moment tensor in x=east, y=north, z=vertical coordinates
mxx = pnt(10)*pnt(1)*pnt(1) + pnt(11)*pnt(4)*pnt(4) + pnt(12)*pnt(7)*pnt(7)
myy = pnt(10)*pnt(2)*pnt(2) + pnt(11)*pnt(5)*pnt(5) + pnt(12)*pnt(8)*pnt(8)
mzz = pnt(10)*pnt(3)*pnt(3) + pnt(11)*pnt(6)*pnt(6) + pnt(12)*pnt(9)*pnt(9)
mxy = pnt(10)*pnt(1)*pnt(2) + pnt(11)*pnt(4)*pnt(5) + pnt(12)*pnt(7)*pnt(8)
mxz = pnt(10)*pnt(1)*pnt(3) + pnt(11)*pnt(4)*pnt(6) + pnt(12)*pnt(7)*pnt(9)
myz = pnt(10)*pnt(2)*pnt(3) + pnt(11)*pnt(5)*pnt(6) + pnt(12)*pnt(8)*pnt(9)

! Convert to spherical coordinates
mrr = mzz
mtt = myy
mpp = mxx
mrt = -myz
mrp = mxz
mtp = -mxy

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine pnt2mag(pnt,mag)
!----
! Calculate the moment magnitude from the P and T magnitudes
!----

implicit none

! Arguments
double precision :: pnt(12), mag

! Local variables
double precision :: mom

call pnt2mom(pnt,mom)
call mom2mag(mom,mag)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine pnt2mag_nondc(pnt,mag)
!----
! Calculate the moment magnitude from the P and T magnitudes
!----

implicit none

! Arguments
double precision :: pnt(12), mag

! Local variables
double precision :: mom

call pnt2mom_nondc(pnt,mom)
call mom2mag(mom,mag)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine pnt2mom(pnt,mom)
!----
! Calculate the seismic moment from the P and T magnitudes by adding the absolute values of the
! two largest eigenvalues and dividing by two (Jost and Herrmann, 1989: Equation 19). This is the
! scalar moment defined by the Global Centroid Moment Tensor project.
!----

implicit none

! Arguments
double precision :: pnt(12), mom

mom = 0.5d0*(abs(pnt(12))+abs(pnt(10)))

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine pnt2mom_nondc(pnt,mom)
!----
! Calculate the seismic moment from the P and T magnitudes by adding the sum of squares of the
! eigenvalues, dividing by two, and taking the square root (Jost and Herrmann, 1989: Equation 20;
! from Silver and Jordan, 1982). This is equivalent to pnt2mom for a double couple source. This is
! the scalar moment defined by the USGS, and the value used for scaling beachballs in GMT psmeca.
!----

implicit none

! Arguments
double precision :: pnt(12), mom

mom = sqrt(0.5d0*(pnt(10)*pnt(10)+pnt(11)*pnt(11)+pnt(12)*pnt(12)))

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine pnt2sdr(pnt,str1,dip1,rak1,str2,dip2,rak2)
!----
! Compute the strike, dip, and rake of the best-fitting double couple source corresponding to a
! set of P, N, and T axes. Calculate for both nodal planes.
!----

use trig, only: r2d, d2r

implicit none

! Arguments
double precision :: pnt(12), str1, dip1, rak1, str2, dip2, rak2

! Local variables
double precision :: slip_vec(3), normal_vec(3), str_vec(3)


! Compute the unit normal and unit slip vectors corresponding to the P and T axes
slip_vec = (pnt(7:9)+pnt(1:3))/sqrt(2.0d0)
normal_vec = (pnt(7:9)-pnt(1:3))/sqrt(2.0d0)

! Convention is that normal vector must point up
if (normal_vec(3).lt.0.0d0) then
    slip_vec = -slip_vec
    normal_vec = -normal_vec
endif

! Compute strike from normal vector
str1 = atan2(normal_vec(1),normal_vec(2))*r2d - 90.0d0
if (str1.lt.0.0d0) then
    str1 = str1 + 360.0d0
endif

! Compute the dip from the normal vector
dip1 = atan2(sqrt(normal_vec(1)**2+normal_vec(2)**2),normal_vec(3))*r2d

! Vector parallel to strike
str_vec(1) = sin(str1*d2r)
str_vec(2) = cos(str1*d2r)
str_vec(3) = 0.0d0

! Compute the rake
rak1 = acos(str_vec(1)*slip_vec(1)+str_vec(2)*slip_vec(2))*r2d
if (str_vec(1)*slip_vec(2)-str_vec(2)*slip_vec(1).lt.0.0d0) then
    rak1 = -rak1
endif

! Compute strike, dip, and rake of other nodal plane
call sdr2sdr(str1,dip1,rak1,str2,dip2,rak2)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine pnt2sv(pnt,slip_vec1,slip_vec2)
!----
! Compute the slip vectors for the best-fitting double couple to P, N, and T axes
!----

implicit none

! Arguments
double precision :: pnt(12), slip_vec1(3), slip_vec2(3)

! Local variables
double precision :: str1, dip1, rak1, str2, dip2, rak2

call pnt2sdr(pnt,str1,dip1,rak1,str2,dip2,rak2)
call sdr2sv(str1,dip1,rak1,slip_vec1)
call sdr2sv(str2,dip2,rak2,slip_vec2)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine pnt2ter(pnt,fth,fss,fno)
!----
! Compute the fraction thrust, strike-slip, and normal component from P, N, and T axes.
!
! Frohlich, C. (1992). Triangle diagrams: ternary graphs to display similarity and diversity of
!     earthquake focal mechanisms. Physics of the Earth and Planetary Interiors, 75, 193-198.
!----

implicit none

! Arguments
double precision :: pnt(12), fth, fss, fno

! Local variables
double precision :: p_pl, n_pl, t_pl

! Compute the plunge of each axis
p_pl = atan(pnt(3)/dsqrt(pnt(1)*pnt(1)+pnt(2)*pnt(2)))
n_pl = atan(pnt(6)/dsqrt(pnt(4)*pnt(4)+pnt(5)*pnt(5)))
t_pl = atan(pnt(9)/dsqrt(pnt(7)*pnt(7)+pnt(8)*pnt(8)))

! Compute the fraction thrust, strike-slip, and normal
fno = dsin(p_pl)*dsin(p_pl)
fss = dsin(n_pl)*dsin(n_pl)
fth = dsin(t_pl)*dsin(t_pl)

return
end

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!---------------------------------- sdr2xxx -------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine sdr2mij(str,dip,rak,mrr,mtt,mpp,mrt,mrp,mtp)
!----
! Compute the unit moment tensor for a double couple source, given the strike, dip, and rake angles.
! Strike is in degrees clockwise from north, dip is in degrees from horizontal, and rake is in
! degrees counter-clockwise from horizontal along the fault plane. The moment tensor is in the
! spherical coordinate system used by the Global Centroid Moment Tensor project.
!----

use trig, only: d2r

implicit none

! Arguments
double precision :: str, dip, rak, mrr, mtt, mpp, mrt, mrp, mtp

! Local variables
double precision :: cs, ss, s2s, c2s, cd, sd, s2d, c2d, cr, sr

! Trigonometric functions of strike, dip, and rake
cs  = dcos(str*d2r)
ss  = dsin(str*d2r)
s2s = dsin(2.0d0*str*d2r)
c2s = dcos(2.0d0*str*d2r)
cd  = dcos(dip*d2r)
sd  = dsin(dip*d2r)
s2d = dsin(2.0d0*dip*d2r)
c2d = dcos(2.0d0*dip*d2r)
cr  = dcos(rak*d2r)
sr  = dsin(rak*d2r)

! Moment tensor components in GCMT convention (r,t,p)
mrr =  (            s2d*sr      )
mtt = -(sd*cr*s2s + s2d*sr*ss*ss)
mpp =  (sd*cr*s2s - s2d*sr*cs*cs)
mrt = -(cd*cr*cs  + c2d*sr*ss   )
mrp =  (cd*cr*ss  - c2d*sr*cs   )
mtp = -(sd*cr*c2s + s2d*sr*ss*cs)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine sdr2pnt(str,dip,rak,pnt)
!----
! Compute the unit P, N, and T axes from a pure double couple source defined by strike, dip, and
! rake angles. There is no moment magnitude information, so the P and T axes should have a value
! of one, and the N axis has a value of zero.
!----

implicit none

! Arguments
double precision :: str, dip, rak, pnt(12)

! Local variables
double precision :: mrr, mtt, mpp, mrt, mrp, mtp

! Convert strike, dip, and rake angles to a moment tensor
call sdr2mij(str,dip,rak,mrr,mtt,mpp,mrt,mrp,mtp)

! Compute eigenvectors and eigenvalues of moment tensor
call mij2pnt(mrr,mtt,mpp,mrt,mrp,mtp,pnt)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine sdr2sdr(str1,dip1,rak1,str2,dip2,rak2)
!----
! Compute the second nodal plane of a double couple source
!----

use trig, only: d2r, r2d

implicit none

! Arguments
double precision :: str1, dip1, rak1, str2, dip2, rak2

! Local variables
double precision :: normal_vec(3), rake_vec(3)


! Compute normal vector of input plane
normal_vec(1) =  dsin(dip1*d2r)*dcos(str1*d2r)
normal_vec(2) = -dsin(dip1*d2r)*dsin(str1*d2r)
normal_vec(3) =  dcos(dip1*d2r)

! Compute rake vector
rake_vec(1) = dcos(rak1*d2r)*dsin(str1*d2r) - dcos(dip1*d2r)*dsin(rak1*d2r)*dcos(str1*d2r)
rake_vec(2) = dcos(rak1*d2r)*dcos(str1*d2r) + dcos(dip1*d2r)*dsin(rak1*d2r)*dsin(str1*d2r)
rake_vec(3) =                                 dsin(dip1*d2r)*dsin(rak1*d2r)

! If rake vector points down, flip it
if (rake_vec(3).lt.0.0d0) then
    rake_vec = -rake_vec
    normal_vec = -normal_vec
endif

! Compute strike of second nodal plane
str2 = datan2(rake_vec(1),rake_vec(2))*r2d - 90.0d0

! Make sure strike is in domain [0,360]
if (str2.gt.360.0d0) then
    str2 = str2 - 360.0d0
endif
if (str2.lt.0.0d0) then
    str2 = str2 + 360.0d0
endif

! Compute dip of second nodal plane
dip2 = datan(dsqrt(rake_vec(1)*rake_vec(1)+rake_vec(2)*rake_vec(2))/rake_vec(3))*r2d

! Compute rake of second nodal plane
rak2 = dacos(normal_vec(1)*dsin(str2*d2r)+normal_vec(2)*dcos(str2*d2r))*r2d

! Check whether rake is positive or negative
if (normal_vec(2)*dsin(str2*d2r)-normal_vec(1)*dcos(str2*d2r).lt.0.0d0) then
    rak2 = -rak2
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine sdr2sv(str,dip,rak,slip_vec)
!----
! Compute the slip vector (east, north, vertical) of the hanging wall relative to the footwall
! corresponding to the strike, dip, and rake.
!----

use trig, only: d2r

implicit none

! Arguments
double precision :: str, dip, rak, slip_vec(3)

! Local variables
double precision :: vec(3)

! Compute slip vector in fault reference frame (x=strike, y=horizontal updip, z=up)
vec(1) = cos(rak*d2r)
vec(2) = sin(rak*d2r)*cos(dip*d2r)
vec(3) = sin(rak*d2r)*sin(dip*d2r)

! Rotate horizontal components to east, north
slip_vec(1) = vec(1)*sin(str*d2r) - vec(2)*cos(str*d2r)
slip_vec(2) = vec(1)*cos(str*d2r) + vec(2)*sin(str*d2r)
slip_vec(3) = vec(3)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine sdr2ter(str,dip,rak,fth,fss,fno)
!----
! Compute the fraction thrust, strike-slip, and normal component from strike, dip, and rake.
!
! Frohlich, C. (1992). Triangle diagrams: ternary graphs to display similarity and diversity of
!     earthquake focal mechanisms. Physics of the Earth and Planetary Interiors, 75, 193-198.
!----

implicit none

! Arguments
double precision :: str, dip, rak, fth, fss, fno

! Local variables
double precision :: pnt(12)

! Compute the P, N, and T axes for this strike, dip, and rake
call sdr2pnt(str,dip,rak,pnt)

! Compute the fraction thrust, strike-slip, and normal
call pnt2ter(pnt,fth,fss,fno)

return
end subroutine

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----------------------------- Moment <-> Magnitude -----------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine mom2mag(mom,mag)
!----
! Convert the scalar seismic moment (in Nm) to moment magnitude
! From Kanamori (1977), p. 2984
!----
implicit none
! Arguments
double precision :: mom, mag
! mag = 2.0d0/3.0d0*log10(mom*nm2dynecm)-10.7d0
mag = 2.0d0/3.0d0*(log10(mom*nm2dynecm)-16.101d0)
return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine mag2mom(mag,mom)
!----
! Convert the moment magnitude to the scalar seismic moment (in Nm)
! From Kanamori (1977), p. 2984
!----
implicit none
! Arguments
double precision :: mag, mom
! mom = 10.0d0**((mag+10.7d0)*1.5d0)
mom = 10.0d0**(1.5d0*mag+16.101d0)
mom = mom*dynecm2nm
return
end subroutine


!--------------------------------------------------------------------------------------------------!
!-------------------------------- ANGLE BETWEEN MOMENT TENSORS ------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine sdr2kagan(str1,dip1,rak1,str2,dip2,rak2,angle)
implicit none
! Arguments
double precision :: str1, dip1, rak1, str2, dip2, rak2, angle
! Local variables
double precision :: mrr1, mtt1, mpp1, mrt1, mrp1, mtp1, mrr2, mtt2, mpp2, mrt2, mrp2, mtp2
double precision :: mt1(6), mt2(6)
call sdr2mij(str1,dip1,rak1,mrr1,mtt1,mpp1,mrt1,mrp1,mtp1)
call sdr2mij(str2,dip2,rak2,mrr2,mtt2,mpp2,mrt2,mrp2,mtp2)
mt1(1) = mrr1
mt1(2) = mtt1
mt1(3) = mpp1
mt1(4) = mrt1
mt1(5) = mrp1
mt1(6) = mtp1
mt2(1) = mrr2
mt2(2) = mtt2
mt2(3) = mpp2
mt2(4) = mrt2
mt2(5) = mrp2
mt2(6) = mtp2
call kaganangle(mt1,mt2,angle)
return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine mij2kagan(mrr1,mtt1,mpp1,mrt1,mrp1,mtp1,mrr2,mtt2,mpp2,mrt2,mrp2,mtp2,angle)
implicit none
! Arguments
double precision :: mrr1, mtt1, mpp1, mrt1, mrp1, mtp1, mrr2, mtt2, mpp2, mrt2, mrp2, mtp2, angle
! Local variables
double precision :: mt1(6), mt2(6)
mt1(1) = mrr1
mt1(2) = mtt1
mt1(3) = mpp1
mt1(4) = mrt1
mt1(5) = mrp1
mt1(6) = mtp1
mt2(1) = mrr2
mt2(2) = mtt2
mt2(3) = mpp2
mt2(4) = mrt2
mt2(5) = mrp2
mt2(6) = mtp2
call kaganangle(mt1,mt2,angle)
return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine kaganangle(mt1,mt2,angle)

use trig, only: r2d

implicit none
double precision :: mt1(6), mt2(6), angle
double precision :: pnt1(12), pnt2(12)
double precision :: mat1(3,3), mat2(3,3), u(3,3), w, qw(4), qq
integer :: i

call mij2pnt(mt1(1),mt1(2),mt1(3),mt1(4),mt1(5),mt1(6),pnt1)
call mij2pnt(mt2(1),mt2(2),mt2(3),mt2(4),mt2(5),mt2(6),pnt2)

! Transpose of first eigenvector matrix
mat1(1,1) = pnt1(1)
mat1(1,2) = pnt1(2)
mat1(1,3) = pnt1(3)
mat1(2,1) = pnt1(4)
mat1(2,2) = pnt1(5)
mat1(2,3) = pnt1(6)
mat1(3,1) = pnt1(7)
mat1(3,2) = pnt1(8)
mat1(3,3) = pnt1(9)

! Second eigenvector matrix
mat2(1,1) = pnt2(1)
mat2(2,1) = pnt2(2)
mat2(3,1) = pnt2(3)
mat2(1,2) = pnt2(4)
mat2(2,2) = pnt2(5)
mat2(3,2) = pnt2(6)
mat2(1,3) = pnt2(7)
mat2(2,3) = pnt2(8)
mat2(3,3) = pnt2(9)

! Multiply matrices
u(1,1) = mat1(1,1)*mat2(1,1) + mat1(1,2)*mat2(2,1) + mat1(1,3)*mat2(3,1)
u(1,2) = mat1(1,1)*mat2(1,2) + mat1(1,2)*mat2(2,2) + mat1(1,3)*mat2(3,2)
u(1,3) = mat1(1,1)*mat2(1,3) + mat1(1,2)*mat2(2,3) + mat1(1,3)*mat2(3,3)
u(2,1) = mat1(2,1)*mat2(1,1) + mat1(2,2)*mat2(2,1) + mat1(2,3)*mat2(3,1)
u(2,2) = mat1(2,1)*mat2(1,2) + mat1(2,2)*mat2(2,2) + mat1(2,3)*mat2(3,2)
u(2,3) = mat1(2,1)*mat2(1,3) + mat1(2,2)*mat2(2,3) + mat1(2,3)*mat2(3,3)
u(3,1) = mat1(3,1)*mat2(1,1) + mat1(3,2)*mat2(2,1) + mat1(3,3)*mat2(3,1)
u(3,2) = mat1(3,1)*mat2(1,2) + mat1(3,2)*mat2(2,2) + mat1(3,3)*mat2(3,2)
u(3,3) = mat1(3,1)*mat2(1,3) + mat1(3,2)*mat2(2,3) + mat1(3,3)*mat2(3,3)

! Calculate quaternion
w = 0.5d0*sqrt(1.0d0 + u(1,1) + u(2,2) + u(3,3))
qw(1) = w
qw(2) = (u(3,2)-u(2,3))/(4.0d0*w)
qw(3) = (u(1,3)-u(3,1))/(4.0d0*w)
qw(4) = (u(2,1)-u(1,2))/(4.0d0*w)

! Take maximum value and use this to calculate the angle
qq = qw(1)
do i = 2,4
    if (abs(qw(i)).gt.abs(qq)) then
        qq = qw(i)
    endif
enddo

angle = 2.0d0*acos(abs(qq))*r2d

return
end subroutine



!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!-------------------------------- EMPIRICAL RELATIONS ---------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine empirical(mag,wid,len,relation_input,fault_type)
!----
! Compute the fault width (down-dip) and length (along-strike) from the magnitude for a particular
! empirical relation.
!----

use io, only: stderr, stdout

implicit none

! Arguments
double precision :: mag, wid, len
character(len=*) :: relation_input, fault_type

! Local variables
integer :: i
logical :: printRelation
double precision :: mom, log10mom
double precision :: len_surf_rupt(4), len_total(4), wid_total(4), area(4), &
                    slip_max(4), slip_avg(4)
character(len=32) :: relation_name


! Initialize return values
wid = 0.0d0
len = 0.0d0

! Initialize values
wid_total = 0.0d0
len_total = 0.0d0


! If relation_name ends with ":print", then print all parameter information for empirical relation
! Otherwise, extract the length and width corresponding to the relation and fault type only
i = index(relation_input,':print')
if (i.eq.0) then
    relation_name = relation_input
    printRelation = .false.
else
    relation_name = relation_input(1:i-1)
    printRelation = .true.
endif


! Determine earthquake parameters and print if requested
if (relation_name.eq.'wells_coppersmith_94'.or.relation_name.eq.'WC') then
    ! Wells, D.L., Coppersmith, K.J. (1994). New empirical relationships among magnitude, rupture
    !   length, rupture width, rupture area, and surface displacement. Bulletin of the Seismological
    !   Society of America, 84, 974–1002.

    ! 1=ss  2=th  3=no  4=all
    ! Total subsurface rupture length
    len_total(1) = 10.0d0 ** (-2.57d0 + 0.62d0*mag)         ! a: +/-0.12,  b: +/-0.02
    len_total(2) = 10.0d0 ** (-2.42d0 + 0.58d0*mag)         ! a: +/-0.21,  b: +/-0.03
    len_total(3) = 10.0d0 ** (-1.88d0 + 0.50d0*mag)         ! a: +/-0.37,  b: +/-0.06
    len_total(4) = 10.0d0 ** (-2.44d0 + 0.59d0*mag)         ! a: +/-0.11,  b: +/-0.02

    ! Downdip rupture width
    wid_total(1) = 10.0d0 ** (-0.76d0 + 0.27d0*mag)         ! a: +/-0.12,  b: +/-0.02
    wid_total(2) = 10.0d0 ** (-1.61d0 + 0.41d0*mag)         ! a: +/-0.20,  b: +/-0.03
    wid_total(3) = 10.0d0 ** (-1.14d0 + 0.35d0*mag)         ! a: +/-0.28,  b: +/-0.05
    wid_total(4) = 10.0d0 ** (-1.01d0 + 0.32d0*mag)         ! a: +/-0.10,  b: +/-0.02

    if (printRelation) then
        ! Surface rupture length
        len_surf_rupt(1) = 10.0d0 ** (-3.55d0 + 0.74d0*mag) ! a: +/-0.37,  b: +/-0.05
        len_surf_rupt(2) = 10.0d0 ** (-2.86d0 + 0.63d0*mag) ! a: +/-0.55,  b: +/-0.08
        len_surf_rupt(3) = 10.0d0 ** (-2.01d0 + 0.50d0*mag) ! a: +/-0.65,  b: +/-0.10
        len_surf_rupt(4) = 10.0d0 ** (-3.22d0 + 0.69d0*mag) ! a: +/-0.27,  b: +/-0.04

        ! Total rupture area
        area(1) = 10.0d0 ** (-3.42d0 + 0.90d0*mag)          ! a: +/-0.18,  b: +/-0.03
        area(2) = 10.0d0 ** (-3.99d0 + 0.98d0*mag)          ! a: +/-0.36,  b: +/-0.06
        area(3) = 10.0d0 ** (-2.87d0 + 0.82d0*mag)          ! a: +/-0.50,  b: +/-0.08
        area(4) = 10.0d0 ** (-3.49d0 + 0.91d0*mag)          ! a: +/-0.16,  b: +/-0.03

        ! Maximum slip
        slip_max(1) = 10.0d0 ** (-7.03d0 + 1.03d0*mag)      ! a: +/-0.55,  b: +/-0.08
        slip_max(2) = 10.0d0 ** (-1.84d0 + 0.29d0*mag)      ! a: +/-1.14,  b: +/-0.17
        slip_max(3) = 10.0d0 ** (-5.90d0 + 0.89d0*mag)      ! a: +/-1.18,  b: +/-0.18
        slip_max(4) = 10.0d0 ** (-5.46d0 + 0.82d0*mag)      ! a: +/-0.51,  b: +/-0.08

        ! Average slip
        slip_avg(1) = 10.0d0 ** (-6.32d0 + 0.90d0*mag)      ! a: +/-0.61,  b: +/-0.09
        slip_avg(2) = 10.0d0 ** (-0.74d0 + 0.08d0*mag)      ! a: +/-1.40,  b: +/-0.21
        slip_avg(3) = 10.0d0 ** (-4.45d0 + 0.63d0*mag)      ! a: +/-1.59,  b: +/-0.24
        slip_avg(4) = 10.0d0 ** (-4.80d0 + 0.69d0*mag)      ! a: +/-0.57,  b: +/-0.08

        write(stdout,*) 'Scaling Relation: Wells and Coppersmith, BSSA (1994)'
        write(stdout,*)
        write(stdout,*) 'Description:'
        write(stdout,*) '"A worldwide data base of source parameters for 421 historical earthquakes ',&
                   'is compiled for this study.'
        write(stdout,*) ' The data include shallow-focus (hypocentral depth less than 40 km), ',&
                   'continental interplate or intraplate'
        write(stdout,*) ' earthquakes of magnitudes greater than approximately 4.5. Earthquakes ',&
                   'associated with subduction zones,'
        write(stdout,*) ' both plate interface earthquakes and those occurring within oceanic slabs, ',&
                   'are excluded."'
        write(stdout,*)
        write(stdout,5062) 'Magnitude:',mag
        write(stdout,*)
        write(stdout,*) 'Total subsurface length (km)'
        write(stdout,5102) 'strike-slip:',len_total(1)
        write(stdout,5102) 'thrust:     ',len_total(2)
        write(stdout,5102) 'normal:     ',len_total(3)
        write(stdout,5102) 'all:        ',len_total(4)
        write(stdout,*) 'Total down-dip width (km)'
        write(stdout,5102) 'strike-slip:',wid_total(1)
        write(stdout,5102) 'thrust:     ',wid_total(2)
        write(stdout,5102) 'normal:     ',wid_total(3)
        write(stdout,5102) 'all:        ',wid_total(4)
        write(stdout,*) 'Surface rupture length (km)'
        write(stdout,5102) 'strike-slip:',len_surf_rupt(1)
        write(stdout,5102) 'thrust:     ',len_surf_rupt(2)
        write(stdout,5102) 'normal:     ',len_surf_rupt(3)
        write(stdout,5102) 'all:        ',len_surf_rupt(4)
        write(stdout,*) 'Total rupture area (km^2)'
        write(stdout,5101) 'strike-slip:',area(1)
        write(stdout,5101) 'thrust:     ',area(2)
        write(stdout,5101) 'normal:     ',area(3)
        write(stdout,5101) 'all:        ',area(4)
        write(stdout,*) 'Maximum slip (m)'
        write(stdout,5103) 'strike-slip:',slip_max(1)
        write(stdout,5103) 'thrust:     ',slip_max(2)
        write(stdout,5103) 'normal:     ',slip_max(3)
        write(stdout,5103) 'all:        ',slip_max(4)
        write(stdout,*) 'Average slip (m)'
        write(stdout,5103) 'strike-slip:',slip_avg(1)
        write(stdout,5103) 'thrust:     ',slip_avg(2)
        write(stdout,5103) 'normal:     ',slip_avg(3)
        write(stdout,5103) 'all:        ',slip_avg(4)
    endif

    ! Return total fault length and width
    if (fault_type.eq.'strike-slip'.or.fault_type.eq.'ss') then
        len = len_total(1)
        wid = wid_total(1)
    elseif (fault_type.eq.'thrust'.or.fault_type.eq.'th') then
        len = len_total(2)
        wid = wid_total(2)
    elseif (fault_type.eq.'normal'.or.fault_type.eq.'no') then
        len = len_total(3)
        wid = wid_total(3)
    else
        len = len_total(4)
        wid = wid_total(4)
    endif


elseif (relation_name.eq.'mai_beroza_00'.or.relation_name.eq.'MB') then
    ! Mai, P.M., Beroza, G.C. (2000). Source scaling properties from finite- fault-rupture models.
    !   Bulletin of the Seismological Society of America, 90, 604-615.
    call mag2mom(mag,mom)
    log10mom = log10(mom)

    ! 1=ss  2=th  3=no  4=all
    ! Effective subsurface rupture length
    len_total(1) = 10.0d0**(-6.31d0+0.40d0*log10mom)          ! a: +/-1.15,  b: +/-0.06
    len_total(2) = 10.0d0**(-6.39d0+0.40d0*log10mom)          ! a: +/-1.04,  b: +/-0.05
    len_total(3) = 10.0d0**(-6.39d0+0.40d0*log10mom)          ! a: +/-1.04,  b: +/-0.05
    len_total(4) = 10.0d0**(-6.13d0+0.39d0*log10mom)          ! a: +/-0.74,  b: +/-0.04

    ! Effective downdip rupture width
    wid_total(1) = 10.0d0**(-2.18d0+0.17d0*log10mom)          ! a: +/-1.09,  b: +/-0.06
    wid_total(2) = 10.0d0**(-5.51d0+0.35d0*log10mom)          ! a: +/-0.81,  b: +/-0.04
    wid_total(3) = 10.0d0**(-5.51d0+0.35d0*log10mom)          ! a: +/-0.81,  b: +/-0.04
    wid_total(4) = 10.0d0**(-5.05d0+0.32d0*log10mom)          ! a: +/-0.74,  b: +/-0.04

    if (printRelation) then
        ! Effective rupture area
        area(1) = 10.0d0**(-8.49d0+0.57d0*log10mom)            ! a: +/-1.85,  b: +/-0.10
        area(2) = 10.0d0**(-11.90d0+0.75d0*log10mom)           ! a: +/-1.67,  b: +/-0.09
        area(3) = 10.0d0**(-11.90d0+0.75d0*log10mom)           ! a: +/-1.67,  b: +/-0.09
        area(4) = 10.0d0**(-11.18d0+0.72d0*log10mom)           ! a: +/-1.18,  b: +/-0.06

        ! Effective slip
        slip_avg(1) = 10.0d0**(-6.03d0+0.43d0*log10mom)/1.0d2  ! a: +/-0.61,  b: +/-0.09
        slip_avg(2) = 10.0d0**(-2.62d0+0.25d0*log10mom)/1.0d2  ! a: +/-1.40,  b: +/-0.21
        slip_avg(3) = 10.0d0**(-2.62d0+0.25d0*log10mom)/1.0d2  ! a: +/-1.59,  b: +/-0.24
        slip_avg(4) = 10.0d0**(-3.34d0+0.29d0*log10mom)/1.0d2  ! a: +/-0.57,  b: +/-0.08


        write(stdout,*) 'Scaling Relation: Mai and Beroza, BSSA (2000)'
        write(stdout,*)
        write(stdout,*) 'Description:'
        write(stdout,*) '"Finite-source images of earthquake rupture show that fault slip is spatially ',&
                   'variable at all resolvable'
        write(stdout,*) ' scales. In this study we develop scaling laws that account for this variability by ',&
                   'measuring effective'
        write(stdout,*) ' fault dimensions derived from the autocorrelation of the slip function for ',&
                   '31 published slip models '
        write(stdout,*) ' of 18 earthquakes, 8 strike-slip events, and 10 dip-slip (reverse, normal, ',&
                   'or oblique) events."'
        write(stdout,*)
        write(stdout,5062) 'Magnitude:',mag
        write(stdout,*)
        write(stdout,*) 'Effective subsurface length (km)'
        write(stdout,5102) 'strike-slip:',len_total(1)
        write(stdout,5102) 'thrust:     ',len_total(2)
        write(stdout,5102) 'normal:     ',len_total(3)
        write(stdout,5102) 'all:        ',len_total(4)
        write(stdout,*) 'Effective down-dip width (km)'
        write(stdout,5102) 'strike-slip:',wid_total(1)
        write(stdout,5102) 'thrust:     ',wid_total(2)
        write(stdout,5102) 'normal:     ',wid_total(3)
        write(stdout,5102) 'all:        ',wid_total(4)
        write(stdout,*) 'Effective rupture area (km^2)'
        write(stdout,5101) 'strike-slip:',area(1)
        write(stdout,5101) 'thrust:     ',area(2)
        write(stdout,5101) 'normal:     ',area(3)
        write(stdout,5101) 'all:        ',area(4)
        write(stdout,*) 'Effective slip (m)'
        write(stdout,5103) 'strike-slip:',slip_avg(1)
        write(stdout,5103) 'thrust:     ',slip_avg(2)
        write(stdout,5103) 'normal:     ',slip_avg(3)
        write(stdout,5103) 'all:        ',slip_avg(4)
    endif

    ! Return effective fault length and width
    if (fault_type.eq.'strike-slip'.or.fault_type.eq.'ss') then
        len = len_total(1)
        wid = wid_total(1)
    elseif (fault_type.eq.'thrust'.or.fault_type.eq.'th') then
        len = len_total(2)
        wid = wid_total(2)
    elseif (fault_type.eq.'normal'.or.fault_type.eq.'no') then
        len = len_total(3)
        wid = wid_total(3)
    else
        len = len_total(4)
        wid = wid_total(4)
    endif


elseif (relation_name.eq.'blaser_et_al_10'.or.relation_name.eq.'B') then
    ! Blaser, L., Kruger, F., Ohrnberger, M., Scherbaum, F. (2010). Scaling relations of earthquake
    !   source parameter estimates with special focus on subduction environment. Bulletin of the
    !   Seismological Society of America, vol. 100, no. 6, pp. 2914-2926.

    ! 1=ss  2=th  3=no  4=all
    ! Subsurface rupture length
    len_total(1) = 10.0d0**(-2.69d0+0.64d0*mag)     ! a: +/-0.11,  b: +/-0.02
    len_total(2) = 10.0d0**(-2.37d0+0.57d0*mag)     ! a: +/-0.13,  b: +/-0.02
    len_total(3) = 10.0d0**(-1.91d0+0.52d0*mag)     ! a: +/-0.29,  b: +/-0.04
    len_total(4) = 10.0d0**(-2.31d0+0.57d0*mag)     ! a: +/-0.08,  b: +/-0.01

    ! Downdip rupture width
    wid_total(1) = 10.0d0**(-1.12d0+0.33d0*mag)     ! a: +/-0.12,  b: +/-0.02
    wid_total(2) = 10.0d0**(-1.86d0+0.46d0*mag)     ! a: +/-0.12,  b: +/-0.02
    wid_total(3) = 10.0d0**(-1.20d0+0.36d0*mag)     ! a: +/-0.25,  b: +/-0.04
    wid_total(4) = 10.0d0**(-1.56d0+0.41d0*mag)     ! a: +/-0.09,  b: +/-0.01

    if (printRelation) then
        write(stdout,*) 'Scaling Relation: Blaser et al., BSSA (2010)'
        write(stdout,*)
        write(stdout,*) 'Description:'
        write(stdout,*) '"Within this study, we compiled a large database of source parameter ',&
                        'estimates of 283 earthquakes....'
        write(stdout,*) ' We distinguish between reverse faulting subduction zone events, ',&
                        'neighboring outer-rise events of'
        write(stdout,*) ' normal faulting slip type and oceanic strike-slip faulting earthquakes."'
        write(stdout,*)
        write(stdout,5062) 'Magnitude:',mag
        write(stdout,*)
        write(stdout,*) 'Subsurface length (km)'
        write(stdout,5102) 'strike-slip:',len_total(1)
        write(stdout,5102) 'thrust:     ',len_total(2)
        write(stdout,5102) 'normal:     ',len_total(3)
        write(stdout,5102) 'all:        ',len_total(4)
        write(stdout,*) 'Down-dip width (km)'
        write(stdout,5102) 'strike-slip:',wid_total(1)
        write(stdout,5102) 'thrust:     ',wid_total(2)
        write(stdout,5102) 'normal:     ',wid_total(3)
        write(stdout,5102) 'all:        ',wid_total(4)
    endif

    if (fault_type.eq.'strike-slip'.or.fault_type.eq.'ss') then
        len = len_total(1)
        wid = wid_total(1)
    elseif (fault_type.eq.'thrust'.or.fault_type.eq.'th') then
        len = len_total(2)
        wid = wid_total(2)
    elseif (fault_type.eq.'normal'.or.fault_type.eq.'no') then
        len = len_total(3)
        wid = wid_total(3)
    else
        len = len_total(4)
        wid = wid_total(4)
    endif


elseif (relation_name.eq.'yen_ma_11'.or.relation_name.eq.'YM') then
    ! Yen, Y.-T., Ma, K.-F. (2011). Source-scaling relationship for M 4.6-8.9 earthquakes,
    !   specifically for earthquakes in the collision zone of Taiwan. Bulletin of the Seismological
    !   Society of America, 101, 464-481.
    call mag2mom(mag,mom)
    log10mom = log10(mom)

    ! 1=ss  2=th  3=no  4=all
    ! Effective subsurface rupture length
    len_total(1) = 10.0d0**(-8.11d0+0.50d0*log10mom)     ! a: +/-1.20,  b: +/-0.07
    len_total(2) = 10.0d0**(-6.66d0+0.42d0*log10mom)     ! a: +/-1.05,  b: +/-0.06
    len_total(3) = 10.0d0**(-6.66d0+0.42d0*log10mom)     ! a: +/-1.05,  b: +/-0.06
    len_total(4) = 10.0d0**(-7.46d0+0.47d0*log10mom)     ! a: +/-0.77,  b: +/-0.04

    ! Effective downdip rupture width
    wid_total(1) = 10.0d0**(-6.67d0+0.42d0*log10mom)     ! a: +/-1.38,  b: +/-0.08
    wid_total(2) = 10.0d0**(-5.76d0+0.37d0*log10mom)     ! a: +/-1.30,  b: +/-0.07
    wid_total(3) = 10.0d0**(-5.76d0+0.37d0*log10mom)     ! a: +/-1.30,  b: +/-0.07
    wid_total(4) = 10.0d0**(-6.30d0+0.40d0*log10mom)     ! a: +/-0.90,  b: +/-0.05

    if (printRelation) then
        ! Effective rupture area
        area(1) = 10.0d0**(-14.77d0+0.92d0*log10mom)         ! a: +/-2.46,  b: +/-0.14
        area(2) = 10.0d0**(-12.45d0+0.80d0*log10mom)         ! a: +/-2.32,  b: +/-0.13
        area(3) = 10.0d0**(-12.45d0+0.80d0*log10mom)         ! a: +/-2.32,  b: +/-0.13
        area(4) = 10.0d0**(-13.79d0+0.87d0*log10mom)         ! a: +/-1.63,  b: +/-0.09

        ! Effective slip
        slip_avg(1) = 10.0d0**( 0.36d0+0.08d0*log10mom)/1.0d2      ! a: +/-2.47,  b: +/-0.14
        slip_avg(2) = 10.0d0**(-2.01d0+0.20d0*log10mom)/1.0d2      ! a: +/-2.32,  b: +/-0.13
        slip_avg(3) = 10.0d0**(-2.01d0+0.20d0*log10mom)/1.0d2      ! a: +/-2.32,  b: +/-0.13
        slip_avg(4) = 10.0d0**(-0.65d0+0.13d0*log10mom)/1.0d2      ! a: +/-1.64,  b: +/-0.09

        write(stdout,*) 'Scaling Relation: Yen and Ma, BSSA (2011)'
        write(stdout,*)
        write(stdout,*) 'Description:'
        write(stdout,*) '"We investigated the source scaling of earthquakes (Mw 4.6–8.9), mostly ',&
                        'from the Taiwan orogenic belt,'
        write(stdout,*) ' and made a global compilation of source parameters to examine the ',&
                        'scaling self-similarity. Finite-fault'
        write(stdout,*) ' slip models (12 dip-slip and 7 strike-slip) using mainly dense ',&
                        'strong-motion data and teleseismic data '
        write(stdout,*) ' from Taiwan were utilized. Seven additional earthquakes (M >7) were ',&
                        'included for further examination'
        write(stdout,*) ' of scaling of large events."'
        write(stdout,*)
        write(stdout,5062) 'Magnitude:',mag
        write(stdout,*)
        write(stdout,*) 'Effective subsurface length (km)'
        write(stdout,5102) 'strike-slip:',len_total(1)
        write(stdout,5102) 'thrust:     ',len_total(2)
        write(stdout,5102) 'normal:     ',len_total(3)
        write(stdout,5102) 'all:        ',len_total(4)
        write(stdout,*) 'Effective down-dip width (km)'
        write(stdout,5102) 'strike-slip:',wid_total(1)
        write(stdout,5102) 'thrust:     ',wid_total(2)
        write(stdout,5102) 'normal:     ',wid_total(3)
        write(stdout,5102) 'all:        ',wid_total(4)
        write(stdout,*) 'Effective rupture area (km^2)'
        write(stdout,5101) 'strike-slip:',area(1)
        write(stdout,5101) 'thrust:     ',area(2)
        write(stdout,5101) 'normal:     ',area(3)
        write(stdout,5101) 'all:        ',area(4)
        write(stdout,*) 'Effective slip (m)'
        write(stdout,5103) 'strike-slip:',slip_avg(1)
        write(stdout,5103) 'thrust:     ',slip_avg(2)
        write(stdout,5103) 'normal:     ',slip_avg(3)
        write(stdout,5103) 'all:        ',slip_avg(4)
    endif

    ! Return effective fault length and width
    if (fault_type.eq.'strike-slip'.or.fault_type.eq.'ss') then
        len = len_total(1)
        wid = wid_total(1)
    elseif (fault_type.eq.'thrust'.or.fault_type.eq.'th') then
        len = len_total(2)
        wid = wid_total(2)
    elseif (fault_type.eq.'normal'.or.fault_type.eq.'no') then
        len = len_total(3)
        wid = wid_total(3)
    else
        len = len_total(4)
        wid = wid_total(4)
    endif


elseif (relation_name.eq.'allen_hayes_17'.or.relation_name.eq.'AH') then
    ! Allen, T.I., Hayes, G.P. (2017). Alternative rupture‐scaling relationships for subduction
    !   interface and other offshore environments. Bulletin of the Seismological Society of America
    !   107, 1240–1253.

    ! 1=ss  2=th  3=no  4=all
    ! Effective subsurface rupture length
    len_total(1) = 10.0d0**(-2.81d0+0.63d0*mag)
    len_total(2) = 10.0d0**(-2.90d0+0.63d0*mag)
    len_total(3) = 10.0d0**(-2.87d0+0.63d0*mag)
    len_total(4) = 10.0d0**(-3.03d0+0.63d0*mag)

    ! Effective downdip rupture width
    wid_total(1) = 10.0d0**(-1.39d0+0.35d0*mag)
    if (mag.le.8.67d0) then
        wid_total(2) = 10.0d0**(-1.91d0+0.48d0*mag)
    else
        wid_total(2) = 10.0d0**( 2.29d0+0.00d0*mag)
    endif
    wid_total(3) = 10.0d0**(-1.18d0+0.35d0*mag)
    wid_total(4) = 10.0d0**(-1.01d0+0.35d0*mag)

    if (printRelation) then
        ! Effective rupture area
        area(1) = 10.0d0**(-4.04d0+0.96d0*mag)
        if (mag.le.8.67d0) then
            area(2) = 10.0d0**(-5.62d0+1.22d0*mag)
        else
            area(2) = 10.0d0**( 2.23d0+0.31d0*mag)
        endif
        area(3) = 10.0d0**(-3.89d0+0.96d0*mag)
        area(4) = 10.0d0**(-3.89d0+0.87d0*mag)

        ! Maximum slip
        slip_max(1) = 10.0d0**(-4.39d0+0.71d0*mag)
        slip_max(2) = 10.0d0**(-4.94d0+0.71d0*mag)
        slip_max(3) = 10.0d0**(-4.58d0+0.71d0*mag)
        slip_max(4) = 10.0d0**(-4.73d0+0.71d0*mag)

        ! Average slip
        slip_avg(1) = 10.0d0**(-4.52d0+0.66d0*mag)
        slip_avg(2) = 10.0d0**(-5.05d0+0.66d0*mag)
        slip_avg(3) = 10.0d0**(-4.70d0+0.66d0*mag)
        slip_avg(4) = 10.0d0**(-4.81d0+0.66d0*mag)

        write(stdout,*) 'Scaling Relation: Allen and Hayes, BSSA (2017)'
        write(stdout,*)
        write(stdout,*) 'Description:'
        write(stdout,*) '"Alternative fault-rupture-scaling relationships are developed for Mw 7.1-',&
                        '9.5 subduction interface '
        write(stdout,*) ' earthquakes using a new database of consistently derived finite-fault ',&
                        'rupture models from '
        write(stdout,*) ' teleseismic inversion....'
        write(stdout,*) ' Fault-rupture-scaling relationships are also developed for [Mw 7.2-8.7] ',&
                        'intraslab (within the'
        write(stdout,*) ' subducting slab), extensional outer-rise and offshore strike-slip',&
                        'environments."'
        write(stdout,*)
        write(stdout,5062) 'Magnitude:',mag
        write(stdout,*)
        write(stdout,*) 'Effective subsurface length (km)'
        write(stdout,5102) 'strike-slip:',len_total(1)
        write(stdout,5102) 'thrust:     ',len_total(2)
        write(stdout,5102) 'normal:     ',len_total(3)
        write(stdout,5102) 'all:        ',len_total(4)
        write(stdout,*) 'Effective down-dip width (km)'
        write(stdout,5102) 'strike-slip:',wid_total(1)
        write(stdout,5102) 'thrust:     ',wid_total(2)
        write(stdout,5102) 'normal:     ',wid_total(3)
        write(stdout,5102) 'all:        ',wid_total(4)
        write(stdout,*) 'Effective rupture area (km^2)'
        write(stdout,5101) 'strike-slip:',area(1)
        write(stdout,5101) 'thrust:     ',area(2)
        write(stdout,5101) 'normal:     ',area(3)
        write(stdout,5101) 'all:        ',area(4)
        write(stdout,*) 'Maximum slip (m)'
        write(stdout,5103) 'strike-slip:',slip_max(1)
        write(stdout,5103) 'thrust:     ',slip_max(2)
        write(stdout,5103) 'normal:     ',slip_max(3)
        write(stdout,5103) 'all:        ',slip_max(4)
        write(stdout,*) 'Average slip (m)'
        write(stdout,5103) 'strike-slip:',slip_avg(1)
        write(stdout,5103) 'thrust:     ',slip_avg(2)
        write(stdout,5103) 'normal:     ',slip_avg(3)
        write(stdout,5103) 'all:        ',slip_avg(4)
    endif

    ! Return effective fault length and width
    if (fault_type.eq.'strike-slip'.or.fault_type.eq.'ss') then
        len = len_total(1)
        wid = wid_total(1)
    elseif (fault_type.eq.'thrust'.or.fault_type.eq.'th') then
        len = len_total(2)
        wid = wid_total(2)
    elseif (fault_type.eq.'normal'.or.fault_type.eq.'no') then
        len = len_total(3)
        wid = wid_total(3)
    else
        len = len_total(4)
        wid = wid_total(4)
    endif

else
    write(stderr,*) 'empirical: no relation named "',trim(relation_name),'"; returning zeros'
    len = 0.0d0
    wid = 0.0d0
endif

5062 format(X,A,F6.2)
5101 format(3X,A,F10.1)
5102 format(3X,A,F10.2)
5103 format(3X,A,F10.3)

return
end subroutine

end module
