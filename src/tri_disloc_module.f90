module tri_disloc_module
!----
! Translation of Brendan Meade's Matlab triangular dislocation code into Fortran module
! By Matthew Herman, Utrecht University
! February 2019
! Original Matlab code lines start with "!**"
!
! Please cite: Meade, B.J. (2007). Algorithms for the calculation of exact displacements, strains,
! and stresses for Triangular Dislocation Elements in a uniform elastic half space. Computers &
! Geosciences. doi:10.1016/j.cageo.2006.12.003.
!----

! Angular dislocation inputs
double precision :: y1, y2, y3, a, b, nu, B1, B2, B3

! Common variables in dislocation expressions
double precision :: sinb, cosb, cotb
double precision :: R, R2, R3, R5, Rbar, R2bar, R3bar, R4bar, R5bar
double precision :: z1, z3, y3bar, z1bar, z3bar
double precision :: F, Fbar

! Parameters
double precision, parameter :: pi = datan(1.0d0)*4.0d0
double precision, parameter :: d2r = pi/180.0d0
double precision, parameter :: r2d = 180.0d0/pi

!--------------------------------------------------------------------------------------------------!

contains

!--------------------------------------------------------------------------------------------------!

!**function [U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds)
subroutine tri_disloc_disp(disp, sta_coord, tri_coord, poisson, slip)

!----
! Arguments
! sta_coord(3): station coordinates (formerly sx, sy, sz)
! tri_coord(3,3): triangle coordinates (formerly x, y, z)
! poisson: Poisson's ratio (formerly pr)
! slip(3): strike-slip, dip-slip, tensile-slip (formerly ss, ds, ts)
! disp(3): displacements (formerly U)
! IMPORTANT NOTES:
!   - Input strike-slip is positive RIGHT-LATERAL (Aki & Richards convention is positive left-lateral)
!   - Input dip-slip is positive NORMAL (Aki & Richards convention is positive THRUST)
!   - Output vertical displacement is positive down
!----

!**% CalcTriDisps.m
!**%
!**% Calculates displacements due to slip on a triangular dislocation in an
!**% elastic half space utilizing the Comninou and Dunders (1975) expressions
!**% for the displacements due to an angular dislocation in an elastic half
!**% space.
!**%
!**% Arguments
!**%  sx : x-coordinates of observation points
!**%  sy : y-coordinates of observation points
!**%  sz : z-coordinates of observation points
!**%  x  : x-coordinates of triangle vertices.
!**%  y  : y-coordinates of triangle vertices.
!**%  z  : z-coordinates of triangle vertices.
!**%  pr : Poisson's ratio
!**%  ss : strike slip displacement
!**%  ts : tensile slip displacement
!**%  ds : dip slip displacement
!**%
!**% Returns
!**%  U  : structure containing the displacements (U.x, U.y, U.z)
!**% This paper should and related code should be cited as:
!**% Brendan J. Meade, Algorithms for the calculation of exact
!**% displacements, strains, and stresses for Triangular Dislocation
!**% Elements in a uniform elastic half space, Computers &
!**% Geosciences (2007), doi:10.1016/j.cageo.2006.12.003.
!**%
!**% Use at your own risk and please let me know of any bugs/errors!
!**%
!**% Copyright (c) 2006 Brendan Meade
!**%
!**% Permission is hereby granted, free of charge, to any person obtaining a
!**% copy of this software and associated documentation files (the
!**% "Software"), to deal in the Software without restriction, including
!**% without limitation the rights to use, copy, modify, merge, publish,
!**% distribute, sublicense, and/or sell copies of the Software, and to permit
!**% persons to whom the Software is furnished to do so, subject to the
!**% following conditions:
!**%
!**% The above copyright notice and this permission notice shall be included
!**% in all copies or substantial portions of the Software.
!**%
!**% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
!**% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!**% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
!**% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
!**% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
!**% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
!**% USE OR OTHER DEALINGS IN THE SOFTWARE.

implicit none

double precision :: disp(3), sta_coord(3), tri_coord(3,4), poisson, slip(3)

! Local variables
integer :: i, iLeg
double precision :: u12(3), u13(3), normal(3), strike(3), updip(3), magnitude, slip_vector(3)
double precision :: tmp3x1(3)
double precision :: dx, dy, dz, dh, leg_angle, leg_plunge, beta
double precision :: ds_vector(3), ts_vector(3), ss_vector(3), slip_ss, slip_ts, slip_ds
double precision :: sx1, sy1, sz1, disp1(3), sx2, sy2, sz2, disp2(3), z1, z2
integer :: inOrOut
double precision :: planeDepth

!**% Calculate the slip vector in XYZ coordinates
!**normVec                      = cross([x(2);y(2);z(2)]-[x(1);y(1);z(1)], [x(3);y(3);z(3)]-[x(1);y(1);z(1)]);
!**normVec                      = normVec./norm(normVec);

! Vectors pointing from triangle vertices 1->2 and vertices 1->3
do i = 1,3
    u12(i) = tri_coord(i,2) - tri_coord(i,1)
    u13(i) = tri_coord(i,3) - tri_coord(i,1)
enddo

normal(1) = u12(2)*u13(3) - u12(3)*u13(2)
normal(2) = u12(3)*u13(1) - u12(1)*u13(3)
normal(3) = u12(1)*u13(2) - u12(2)*u13(1)

! Normalize the normal vector
magnitude = dsqrt(normal(1)*normal(1)+normal(2)*normal(2)+normal(3)*normal(3))
do i = 1,3
    normal(i) = normal(i)/magnitude
enddo

!**if (normVec(3) < 0) % Enforce clockwise circulation
!**   normVec                   = -normVec;
!**   [x(2) x(3)]               = swap(x(2), x(3));
!**   [y(2) y(3)]               = swap(y(2), y(3));
!**   [z(2) z(3)]               = swap(z(2), z(3));
!**end

! Enforce normal vector points up by swapping points 2 and 3 if necessary
if (normal(3).lt.0.0d0) then
    tmp3x1 = tri_coord(:,2)
    tri_coord(:,2) = tri_coord(:,3)
    tri_coord(:,3) = tmp3x1
    do i = 1,3
        normal(i) = -normal(i)
    enddo
endif

!**strikeVec                    = [-sin(atan2(normVec(2),normVec(1))) cos(atan2(normVec(2),normVec(1))) 0];
!**dipVec                       = cross(normVec, strikeVec);

! Calculate strike, up-dip vectors of triangle
strike(1) = -dsin(datan2(normal(2),normal(1)))
strike(2) =  dcos(datan2(normal(2),normal(1)))
strike(3) = 0.0d0
updip(1) = normal(2)*strike(3) - normal(3)*strike(2)
updip(2) = normal(3)*strike(1) - normal(1)*strike(3)
updip(3) = normal(1)*strike(2) - normal(2)*strike(1)

!**slipComp                     = [ss ds ts];
!**slipVec                      = [strikeVec(:) dipVec(:) normVec(:)] * slipComp(:);

! Calculate slip vector (1=ss, 2=ds, 3=ts)
slip_vector(1) = strike(1)*slip(1) + updip(1)*slip(2) + normal(1)*slip(3)
slip_vector(2) = strike(2)*slip(1) + updip(2)*slip(2) + normal(2)*slip(3)
slip_vector(3) = strike(3)*slip(1) + updip(3)*slip(2) + normal(3)*slip(3)

!**% Solution vectors
!**U.x                          = zeros(size(sx));
!**U.y                          = zeros(size(sx));
!**U.z                          = zeros(size(sx));

! Initialize displacement solution
disp = 0.0d0

!**% Add a copy of the first vertex to the vertex list for indexing
!**x(4)                         = x(1);
!**y(4)                         = y(1);
!**z(4)                         = z(1);

tri_coord(:,4) = tri_coord(:,1)

!**for iTri = 1:3
do iLeg = 1,3

!**    % Calculate strike and dip of current leg
!**    strike                   = 180/pi*(atan2(y(iTri+1)-y(iTri), x(iTri+1)-x(iTri)));
!**    segMapLength             = sqrt((x(iTri)-x(iTri+1))^2 + (y(iTri)-y(iTri+1))^2);
!**    [rx ry]                  = RotateXyVec(x(iTri+1)-x(iTri), y(iTri+1)-y(iTri), -strike);
!**    dip                      = 180/pi*(atan2(z(iTri+1)-z(iTri), rx));

    ! Azimuth (CCW from X) and plunge of leg
    dx = tri_coord(1,iLeg+1)-tri_coord(1,iLeg)
    dy = tri_coord(2,iLeg+1)-tri_coord(2,iLeg)
    dz = tri_coord(3,iLeg+1)-tri_coord(3,iLeg)
    dh = dsqrt(dx*dx+dy*dy)
    leg_angle = r2d*datan2(dy,dx)
    leg_plunge = r2d*datan2(dz,dh)

!**    if dip >= 0
!**       beta                  = pi/180*(90-dip);
!**       if beta > pi/2
!**          beta               = pi/2-beta;
!**       end
!**    else
!**       beta                  = -pi/180*(90+dip);
!**       if beta < -pi/2
!**          beta = pi/2-abs(beta);
!**       end
!**    end

    ! Adjust plunge angle
    if (leg_plunge.ge.0) then
        beta = d2r*(90.0d0-leg_plunge)
        if (beta.gt.pi/2.0d0) then
            beta = pi/2.0d0 - beta
        endif
    else
        beta = -d2r*(90.0d0+leg_plunge)
        if (beta.lt.-pi/2.0d0) then
            beta = pi/2.0d0 - dabs(beta)
        endif
    endif

!**    ssVec                    = [cos(strike/180*pi) sin(strike/180*pi) 0];
!**    tsVec                    = [-sin(strike/180*pi) cos(strike/180*pi) 0];
!**    dsVec                    = cross(ssVec, tsVec);
!**    lss                      = dot(slipVec, ssVec);
!**    lts                      = dot(slipVec, tsVec);
!**    lds                      = dot(slipVec, dsVec);

    ! Slip vector projected onto leg
    ss_vector(1) = dcos(leg_angle*d2r)
    ss_vector(2) = dsin(leg_angle*d2r)
    ss_vector(3) = 0.0d0
    ts_vector(1) = -dsin(leg_angle*d2r)
    ts_vector(2) = dcos(leg_angle*d2r)
    ts_vector(3) = 0.0d0
    ds_vector(1) = ss_vector(2)*ts_vector(3) - ss_vector(3)*ts_vector(2)
    ds_vector(2) = ss_vector(3)*ts_vector(1) - ss_vector(1)*ts_vector(3)
    ds_vector(3) = ss_vector(1)*ts_vector(2) - ss_vector(2)*ts_vector(1)
    slip_ss = slip_vector(1)*ss_vector(1)+slip_vector(2)*ss_vector(2)+slip_vector(3)*ss_vector(3)
    slip_ts = slip_vector(1)*ts_vector(1)+slip_vector(2)*ts_vector(2)+slip_vector(3)*ts_vector(3)
    slip_ds = slip_vector(1)*ds_vector(1)+slip_vector(2)*ds_vector(2)+slip_vector(3)*ds_vector(3)

!**    if (abs(beta) > 0.000001) && (abs(beta-pi) > 0.000001)
!**       % First angular dislocation
!**       [sx1 sy1]                 = RotateXyVec(sx-x(iTri), sy-y(iTri), -strike);
!**       [ux1 uy1 uz1]             = adv(sx1, sy1, sz-z(iTri), z(iTri), beta, pr, lss, lts, lds);
!**
!**       % Second angular dislocation
!**       [sx2 sy2]                 = RotateXyVec(sx-x(iTri+1), sy-y(iTri+1), -strike);
!**       [ux2 uy2 uz2]             = adv(sx2, sy2, sz-z(iTri+1), z(iTri+1), beta, pr, lss, lts, lds);
!**
!**       % Rotate vectors to correct for strike
!**       [uxn uyn]                 = RotateXyVec(ux1-ux2, uy1-uy2, strike);
!**       uzn                       = uz1-uz2;
!**
!**       % Add the displacements from current leg
!**       U.x                       = U.x + uxn;
!**       U.y                       = U.y + uyn;
!**       U.z                       = U.z + uzn;
!**    end

    ! Calculate angular dislocations for each point in leg
    if (dabs(beta).gt.1.0d-6 .and. dabs(beta-pi).gt.1.0d-6) then

        ! First angular dislocation
        sx1 = dcos(-leg_angle*d2r)*(sta_coord(1)-tri_coord(1,iLeg)) - &
              dsin(-leg_angle*d2r)*(sta_coord(2)-tri_coord(2,iLeg))
        sy1 = dsin(-leg_angle*d2r)*(sta_coord(1)-tri_coord(1,iLeg)) + &
              dcos(-leg_angle*d2r)*(sta_coord(2)-tri_coord(2,iLeg))
        sz1 = sta_coord(3)-tri_coord(3,iLeg)
        z1 = tri_coord(3,iLeg)
        call ang_disloc_disp(disp1,sx1,sy1,sz1,z1,beta,poisson,slip_ss,slip_ts,slip_ds)

        ! Second angular dislocation
        sx2 = dcos(-leg_angle*d2r)*(sta_coord(1)-tri_coord(1,iLeg+1)) - &
              dsin(-leg_angle*d2r)*(sta_coord(2)-tri_coord(2,iLeg+1))
        sy2 = dsin(-leg_angle*d2r)*(sta_coord(1)-tri_coord(1,iLeg+1)) + &
              dcos(-leg_angle*d2r)*(sta_coord(2)-tri_coord(2,iLeg+1))
        sz2 = sta_coord(3)-tri_coord(3,iLeg+1)
        z2 = tri_coord(3,iLeg+1)
        call ang_disloc_disp(disp2,sx2,sy2,sz2,z2,beta,poisson,slip_ss,slip_ts,slip_ds)

        ! Rotate back
        disp(1) = disp(1) + dcos(leg_angle*d2r)*(disp1(1)-disp2(1)) - &
                            dsin(leg_angle*d2r)*(disp1(2)-disp2(2))
        disp(2) = disp(2) + dsin(leg_angle*d2r)*(disp1(1)-disp2(1)) + &
                            dcos(leg_angle*d2r)*(disp1(2)-disp2(2))
        disp(3) = disp(3) + disp1(3) - disp2(3)
    endif

!**end
enddo

!**% Identify indices for stations under current triangle
!**inPolyIdx                       = find(inpolygon(sx, sy, x, y) == 1);
!**underIdx = [];
!**for iIdx = 1 : numel(inPolyIdx)
!**   d                            = LinePlaneIntersect(x, y, z, sx(inPolyIdx(iIdx)), sy(inPolyIdx(iIdx)), sz(inPolyIdx(iIdx)));
!**   if d(3)-sz(inPolyIdx(iIdx)) < 0
!**      underIdx = [underIdx ; inPolyIdx(iIdx)];
!**   end
!**end

!**% Apply static offset to the points that lie underneath the current triangle
!**U.x(underIdx)                = U.x(underIdx) - slipVec(1);
!**U.y(underIdx)                = U.y(underIdx) - slipVec(2);
!**U.z(underIdx)                = U.z(underIdx) - slipVec(3);

! Check for points under triangle
call pnpoly(sta_coord(1),sta_coord(2),tri_coord(1,:),tri_coord(2,:),4,inOrOut)

! If point is under, add slip to displacement
if (inOrOut.ge.0) then
    call verticalLinePlaneIntersect(planeDepth,sta_coord,tri_coord(:,1),normal)
    if (planeDepth.lt.sta_coord(3)) then
        disp(1) = disp(1) - slip_vector(1)
        disp(2) = disp(2) - slip_vector(2)
        disp(3) = disp(3) - slip_vector(3)
    endif
endif


return
end subroutine

!--------------------------------------------------------------------------------------------------!

!**function [S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds)
subroutine tri_disloc_strain(strain_output, sta_coord, tri_coord, poisson, slip)

!----
! Arguments
! sta_coord(3): station coordinates (formerly sx, sy, sz)
! tri_coord(3,3): triangle coordinates (formerly x, y, z)
! poisson: Poisson's ratio (formerly pr)
! slip(3): strike-slip, dip-slip, tensile-slip (formerly ss, ds, ts)
! strain(3,3): strains (formerly S)
! IMPORTANT NOTES:
!   - Input strike-slip is positive RIGHT-LATERAL (Aki & Richards convention is positive left-lateral)
!   - Input dip-slip is positive NORMAL (Aki & Richards convention is positive THRUST)
!   - Output vertical displacement is positive down, so xz and yz are also opposite signs
!----

!**% CalcTriStrains.m
!**%
!**% Calculates strains due to slip on a triangular dislocation in an
!**% elastic half space utilizing the symbolically differentiated
!**% displacement gradient tensor derived from the expressions for
!**% the displacements due to an angular dislocation in an elastic half
!**% space (Comninou and Dunders, 1975).
!**%
!**% Arguments
!**%  sx : x-coordinates of observation points
!**%  sy : y-coordinates of observation points
!**%  sz : z-coordinates of observation points
!**%  x  : x-coordinates of triangle vertices.
!**%  y  : y-coordinates of triangle vertices.
!**%  z  : z-coordinates of triangle vertices.
!**%  pr : Poisson's ratio
!**%  ss : strike slip displacement
!**%  ts : tensile slip displacement
!**%  ds : dip slip displacement
!**%
!**% Returns
!**%  S  : structure containing the strains (S.xx, S.yy, S.zz, S.xy, S.xz, S.yz)
!**%
!**% This paper should and related code should be cited as:
!**% Brendan J. Meade, Algorithms for the calculation of exact
!**% displacements, strains, and stresses for Triangular Dislocation
!**% Elements in a uniform elastic half space, Computers &
!**% Geosciences (2007), doi:10.1016/j.cageo.2006.12.003.
!**%
!**% Use at your own risk and please let me know of any bugs/errors.
!**%
!**% Copyright (c) 2006 Brendan Meade
!**%
!**% Permission is hereby granted, free of charge, to any person obtaining a
!**% copy of this software and associated documentation files (the
!**% "Software"), to deal in the Software without restriction, including
!**% without limitation the rights to use, copy, modify, merge, publish,
!**% distribute, sublicense, and/or sell copies of the Software, and to permit
!**% persons to whom the Software is furnished to do so, subject to the
!**% following conditions:
!**%
!**% The above copyright notice and this permission notice shall be included
!**% in all copies or substantial portions of the Software.
!**%
!**% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
!**% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!**% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
!**% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
!**% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
!**% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
!**% USE OR OTHER DEALINGS IN THE SOFTWARE.

implicit none

double precision :: strain_output(3,3), sta_coord(3), tri_coord(3,4), poisson, slip(3)

! Local variables
integer :: i
double precision :: u12(3), u13(3), magnitude, normal(3), strike(3), updip(3), slip_vector(3)
double precision :: tmp3x1(3)
integer :: iLeg
double precision :: dx, dy, dz, dh, leg_angle, sina, cosa, leg_plunge, beta
double precision :: ds_vector(3), ts_vector(3), ss_vector(3), slip_ss, slip_ts, slip_ds
double precision :: sx1, sy1, sz1, z1, sx2, sy2, sz2, z2
double precision :: stn1(3,3), stn2(3,3), stn_diff(3,3), strain_leg(3,3)

!**% Calculate the slip vector in XYZ coordinates
!**normVec                      = cross([x(2);y(2);z(2)]-[x(1);y(1);z(1)], [x(3);y(3);z(3)]-[x(1);y(1);z(1)]);
!**normVec                      = normVec./norm(normVec);

! Vectors pointing from triangle vertices 1->2 and vertices 1->3
do i = 1,3
    u12(i) = tri_coord(i,2) - tri_coord(i,1)
    u13(i) = tri_coord(i,3) - tri_coord(i,1)
enddo

normal(1) = u12(2)*u13(3) - u12(3)*u13(2)
normal(2) = u12(3)*u13(1) - u12(1)*u13(3)
normal(3) = u12(1)*u13(2) - u12(2)*u13(1)

! Normalize the normal vector
magnitude = dsqrt(normal(1)*normal(1)+normal(2)*normal(2)+normal(3)*normal(3))
do i = 1,3
    normal(i) = normal(i)/magnitude
enddo

!**if (normVec(3) < 0) % Enforce clockwise circulation
!**   normVec                   = -normVec;
!**   [x(2) x(3)]               = swap(x(2), x(3));
!**   [y(2) y(3)]               = swap(y(2), y(3));
!**   [z(2) z(3)]               = swap(z(2), z(3));
!**end

! Enforce normal vector points up by swapping points 2 and 3 if necessary
if (normal(3).lt.0.0d0) then
    tmp3x1 = tri_coord(:,2)
    tri_coord(:,2) = tri_coord(:,3)
    tri_coord(:,3) = tmp3x1
    do i = 1,3
        normal(i) = -normal(i)
    enddo
endif

!**strikeVec                    = [-sin(atan2(normVec(2),normVec(1))) cos(atan2(normVec(2),normVec(1))) 0];
!**dipVec                       = cross(normVec, strikeVec);

! Calculate strike, up-dip vectors of triangle
strike(1) = -dsin(datan2(normal(2),normal(1)))
strike(2) =  dcos(datan2(normal(2),normal(1)))
strike(3) = 0.0d0
updip(1) = normal(2)*strike(3) - normal(3)*strike(2)
updip(2) = normal(3)*strike(1) - normal(1)*strike(3)
updip(3) = normal(1)*strike(2) - normal(2)*strike(1)

!**slipComp                     = [ss ds ts];
!**slipVec                      = [strikeVec(:) dipVec(:) normVec(:)] * slipComp(:);

! Calculate slip vector (1=ss, 2=ds, 3=ts)
slip_vector(1) = strike(1)*slip(1) + updip(1)*slip(2) + normal(1)*slip(3)
slip_vector(2) = strike(2)*slip(1) + updip(2)*slip(2) + normal(2)*slip(3)
slip_vector(3) = strike(3)*slip(1) + updip(3)*slip(2) + normal(3)*slip(3)

!**% Solution vectors
!**S.xx                         = zeros(size(sx));
!**S.yy                         = zeros(size(sx));
!**S.zz                         = zeros(size(sx));
!**S.xy                         = zeros(size(sx));
!**S.xz                         = zeros(size(sx));
!**S.yz                         = zeros(size(sx));

! Initialize output strain tensor
strain_output = 0.0d0

!**% Add a copy of the first vertex to the vertex list for indexing
!**x(4)                         = x(1);
!**y(4)                         = y(1);
!**z(4)                         = z(1);

tri_coord(:,4) = tri_coord(:,1)

!**for iTri = 1:3
do iLeg = 1,3

!**   % Calculate strike and dip of current leg
!**   strike                   = 180/pi*(atan2(y(iTri+1)-y(iTri), x(iTri+1)-x(iTri)));
!**   segMapLength             = sqrt((x(iTri)-x(iTri+1))^2 + (y(iTri)-y(iTri+1))^2);
!**   [rx ry]                  = RotateXyVec(x(iTri+1)-x(iTri), y(iTri+1)-y(iTri), -strike);
!**   dip                      = 180/pi*(atan2(z(iTri+1)-z(iTri), rx));

    ! Azimuth (CCW from X) and plunge of leg
    dx = tri_coord(1,iLeg+1)-tri_coord(1,iLeg)
    dy = tri_coord(2,iLeg+1)-tri_coord(2,iLeg)
    dz = tri_coord(3,iLeg+1)-tri_coord(3,iLeg)
    dh = dsqrt(dx*dx+dy*dy)
    leg_angle = r2d*datan2(dy,dx)
    leg_plunge = r2d*datan2(dz,dh)

    cosa = dcos(leg_angle*d2r)
    sina = dsin(leg_angle*d2r)


!**   if dip >= 0
!**      beta                  = pi/180*(90-dip);
!**      if beta > pi/2
!**         beta               = pi/2-beta;
!**      end
!**   else
!**      beta                  = -pi/180*(90+dip);
!**      if beta < -pi/2
!**         beta = pi/2-abs(beta);
!**      end
!**   end

    ! Adjust plunge angle
    if (leg_plunge.ge.0) then
        beta = d2r*(90.0d0-leg_plunge)
        if (beta.gt.pi/2.0d0) then
            beta = pi/2.0d0 - beta
        endif
    else
        beta = -d2r*(90.0d0+leg_plunge)
        if (beta.lt.-pi/2.0d0) then
            beta = pi/2.0d0 - dabs(beta)
        endif
    endif

!**   ssVec                    = [cos(strike/180*pi) sin(strike/180*pi) 0];
!**   tsVec                    = [-sin(strike/180*pi) cos(strike/180*pi) 0];
!**   dsVec                    = cross(ssVec, tsVec);
!**   lss                      = dot(slipVec, ssVec);
!**   lts                      = dot(slipVec, tsVec);
!**   lds                      = dot(slipVec, dsVec);

    ! Slip vector projected onto leg
    ss_vector(1) = cosa
    ss_vector(2) = sina
    ss_vector(3) = 0.0d0
    ts_vector(1) = -sina
    ts_vector(2) = cosa
    ts_vector(3) = 0.0d0
    ds_vector(1) = ss_vector(2)*ts_vector(3) - ss_vector(3)*ts_vector(2)
    ds_vector(2) = ss_vector(3)*ts_vector(1) - ss_vector(1)*ts_vector(3)
    ds_vector(3) = ss_vector(1)*ts_vector(2) - ss_vector(2)*ts_vector(1)
    slip_ss = slip_vector(1)*ss_vector(1)+slip_vector(2)*ss_vector(2)+slip_vector(3)*ss_vector(3)
    slip_ts = slip_vector(1)*ts_vector(1)+slip_vector(2)*ts_vector(2)+slip_vector(3)*ts_vector(3)
    slip_ds = slip_vector(1)*ds_vector(1)+slip_vector(2)*ds_vector(2)+slip_vector(3)*ds_vector(3)

!**   if (abs(beta) > 0.000001) && (abs(beta-pi) > 0.000001)
    if (dabs(beta).gt.1.0d-6 .and. dabs(beta-pi).gt.1.0d-6) then

!**      % First angular dislocation
!**      [sx1 sy1]                 = RotateXyVec(sx-x(iTri), sy-y(iTri), -strike);
!**      [a11 a22 a33 a12 a13 a23] = advs(sx1, sy1, sz-z(iTri), z(iTri), beta, pr, lss, lts, lds);
!**
!**      % Second angular dislocation
!**      [sx2 sy2]                 = RotateXyVec(sx-x(iTri+1), sy-y(iTri+1), -strike);
!**      [b11 b22 b33 b12 b13 b23] = advs(sx2, sy2, sz-z(iTri+1), z(iTri+1), beta, pr, lss, lts, lds);

        ! First angular dislocation
        sx1 = dcos(-leg_angle*d2r)*(sta_coord(1)-tri_coord(1,iLeg)) - &
              dsin(-leg_angle*d2r)*(sta_coord(2)-tri_coord(2,iLeg))
        sy1 = dsin(-leg_angle*d2r)*(sta_coord(1)-tri_coord(1,iLeg)) + &
              dcos(-leg_angle*d2r)*(sta_coord(2)-tri_coord(2,iLeg))
        sz1 = sta_coord(3)-tri_coord(3,iLeg)
        z1 = tri_coord(3,iLeg)
        call ang_disloc_strain(stn1,sx1,sy1,sz1,z1,beta,poisson,slip_ss,slip_ts,slip_ds)

        ! Second angular dislocation
        sx2 = dcos(-leg_angle*d2r)*(sta_coord(1)-tri_coord(1,iLeg+1)) - &
              dsin(-leg_angle*d2r)*(sta_coord(2)-tri_coord(2,iLeg+1))
        sy2 = dsin(-leg_angle*d2r)*(sta_coord(1)-tri_coord(1,iLeg+1)) + &
              dcos(-leg_angle*d2r)*(sta_coord(2)-tri_coord(2,iLeg+1))
        sz2 = sta_coord(3)-tri_coord(3,iLeg+1)
        z2 = tri_coord(3,iLeg+1)
        call ang_disloc_strain(stn2,sx2,sy2,sz2,z2,beta,poisson,slip_ss,slip_ts,slip_ds)

!**      % Rotate tensors to correct for strike
!**      bxx                       = a11-b11;
!**      byy                       = a22-b22;
!**      bzz                       = a33-b33;
!**      bxy                       = a12-b12;
!**      bxz                       = a13-b13;
!**      byz                       = a23-b23;

        stn_diff = stn1-stn2

!**      g                         = pi/180*strike;
        ! Already defined rotation matrix components as cosa, sina

!**      e11n                      = (cos(g)*bxx-sin(g)*bxy)*cos(g)-(cos(g)*bxy-sin(g)*byy)*sin(g);
!**      e12n                      = (cos(g)*bxx-sin(g)*bxy)*sin(g)+(cos(g)*bxy-sin(g)*byy)*cos(g);
!**      e13n                      = cos(g)*bxz-sin(g)*byz;
!**      e22n                      = (sin(g)*bxx+cos(g)*bxy)*sin(g)+(sin(g)*bxy+cos(g)*byy)*cos(g);
!**      e23n                      = sin(g)*bxz+cos(g)*byz;
!**      e33n                      = bzz;

        strain_leg = 0.0d0
        strain_leg(1,1) = (cosa*stn_diff(1,1)-sina*stn_diff(1,2))*cosa - &
                          (cosa*stn_diff(1,2)-sina*stn_diff(2,2))*sina
        strain_leg(1,2) = (cosa*stn_diff(1,1)-sina*stn_diff(1,2))*sina + &
                          (cosa*stn_diff(1,2)-sina*stn_diff(2,2))*cosa
        strain_leg(1,3) = cosa*stn_diff(1,3)-sina*stn_diff(2,3)
        strain_leg(2,1) = strain_leg(1,2)
        strain_leg(2,2) = (sina*stn_diff(1,1)+cosa*stn_diff(1,2))*sina + &
                          (sina*stn_diff(1,2)+cosa*stn_diff(2,2))*cosa
        strain_leg(2,3) = sina*stn_diff(1,3)+cosa*stn_diff(2,3)
        strain_leg(3,1) = strain_leg(1,3)
        strain_leg(3,2) = strain_leg(2,3)
        strain_leg(3,3) = stn_diff(3,3)

!**      % Add the strains from current leg
!**      S.xx                      = S.xx + e11n;
!**      S.yy                      = S.yy + e22n;
!**      S.zz                      = S.zz + e33n;
!**      S.xy                      = S.xy + e12n;
!**      S.xz                      = S.xz + e13n;
!**      S.yz                      = S.yz + e23n;

        strain_output = strain_output + strain_leg

    endif
!**end
enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

!**function [v1 v2 v3] = adv(y1, y2, y3, a, beta, nu, B1, B2, B3)
subroutine ang_disloc_disp(disp,y1in,y2in,y3in,ain,bin,nuin,B1in,B2in,B3in)

!**% These are the displacements in a uniform elastic half space due to slip
!**% on an angular dislocation (Comninou and Dunders, 1975).  Some of the
!**% equations for the B2 and B3 cases have been corrected following Thomas
!**% 1993.  The equations are coded in way such that they roughly correspond
!**% to each line in original text.  Exceptions have been made where it made
!**% more sense because of grouping symbols.

implicit none
! Arguments
double precision :: disp(3), y1in, y2in, y3in, ain, bin, nuin, B1in, B2in, B3in

! THESE VALUES GET CALCULATED FOR STRAINS, TOO. MOVED TO THEIR OWN SUBROUTINE....
!**sinbeta           = sin(beta);
!**cosbeta           = cos(beta);
!**cotbeta           = cot(beta);
!**z1                = y1.*cosbeta - y3.*sinbeta;
!**z3                = y1.*sinbeta + y3.*cosbeta;
!**R2                = y1.*y1 + y2.*y2 + y3.*y3;
!**R                 = sqrt(R2);
!**y3bar             = y3 + 2.*a;
!**z1bar             = y1.*cosbeta + y3bar.*sinbeta;
!**z3bar             = -y1.*sinbeta + y3bar.*cosbeta;
!**R2bar             = y1.*y1 + y2.*y2 + y3bar.*y3bar;
!**Rbar              = sqrt(R2bar);
!**F                 = -atan2(y2, y1) + atan2(y2, z1) + atan2(y2.*R.*sinbeta, y1.*z1+(y2.*y2).*cosbeta);
!**Fbar              = -atan2(y2, y1) + atan2(y2, z1bar) + atan2(y2.*Rbar.*sinbeta, y1.*z1bar+(y2.*y2).*cosbeta);

y1 = y1in
y2 = y2in
y3 = y3in
a = ain
b = bin
nu = nuin
B1 = B1in
B2 = B2in
B3 = B3in
call angular_disloc_vars()

!% Sum the for each slip component
disp(1) = B1*v1B1() + B2*v1B2() + B3*v1B3()
disp(2) = B1*v2B1() + B2*v2B2() + B3*v2B3()
disp(3) = B1*v3B1() + B2*v3B2() + B3*v3B3()

return
end subroutine

!--------------------------------------------------------------------------------------------------!

!**function [e11 e22 e33 e12 e13 e23] = advs(y1, y2, y3, a, b, nu, B1, B2, B3)
subroutine ang_disloc_strain(ang_disloc_stn, y1in, y2in, y3in, ain, bin, nuin, B1in, B2in, B3in)

!----
! Compute the strain tensor corresponding to an angular dislocation in an elastic half-space
! In the off chance that someone wants to use this subroutine directly to compute an
! angular dislocation (e.g., to build non-triangular elements), arguments are passed
! rather than using hidden module variables. They are then explicitly saved as the module variables.
!----

!**% These are the strains in a uniform elastic half space due to slip
!**% on an angular dislocation.  They were calculated by symbolically
!**% differentiating the expressions for the displacements (Comninou and
!**% Dunders, 1975, with typos noted by Thomas 1993) then combining the
!**% elements of the displacement gradient tensor to form the strain tensor.

implicit none
! Arguments
double precision :: ang_disloc_stn(3,3), y1in, y2in, y3in, ain, bin, nuin, B1in, B2in, B3in

! Save arguments to module variables
y1 = y1in
y2 = y2in
y3 = y3in
a = ain
b = bin
nu = nuin
B1 = B1in
B2 = B2in
B3 = B3in

! Pre-compute common values in angular dislocation expressions
call angular_disloc_vars()

ang_disloc_stn(1,1) = B1*(e11_termB1a()+e11_termB1b()) + &
                      B2*(e11_termB2a()+e11_termB2b()) + &
                      B3*(e11_termB3a()+e11_termB3b())

ang_disloc_stn(2,2) = B1*(e22_termB1a()+e22_termB1b()) + &
                      B2*(e22_termB2a()+e22_termB2b()) + &
                      B3*(e22_termB3a()+e22_termB3b())

ang_disloc_stn(3,3) = B1*(e33_termB1a()+e33_termB1b()) + &
                      B2*(e33_termB2a()+e33_termB2b()) + &
                      B3*(e33_termB3a()+e33_termB3b())

ang_disloc_stn(1,2) = B1*(e12_termB1a()+e12_termB1b()) + &
                      B2*(e12_termB2a()+e12_termB2b()) + &
                      B3*(e12_termB3a()+e12_termB3b())
ang_disloc_stn(2,1) = ang_disloc_stn(1,2)

ang_disloc_stn(1,3) = B1*(e13_termB1a()+e13_termB1b()) + &
                      B2*(e13_termB2a()+e13_termB2b()) + &
                      B3*(e13_termB3a()+e13_termB3b())
ang_disloc_stn(3,1) = ang_disloc_stn(1,3)

ang_disloc_stn(2,3) = B1*(e23_termB1a()+e23_termB1b()) + &
                      B2*(e23_termB2a()+e23_termB2b()) + &
                      B3*(e23_termB3a()+e23_termB3b())
ang_disloc_stn(3,2) = ang_disloc_stn(2,3)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine angular_disloc_vars()
implicit none
sinb = dsin(b)
cosb = dcos(b)
cotb = cosb/sinb
R2 = y1*y1 + y2*y2 + y3*y3
R  = dsqrt(R2)
R3 = R*R*R
R5 = R*R*R*R*R
z1 = y1*cosb - y3*sinb
z3 = y1*sinb + y3*cosb
y3bar = y3 + 2.0d0*a
z1bar = y1*cosb + y3bar*sinb
z3bar = -y1*sinb + y3bar*cosb
R2bar = y1*y1 + y2*y2 + y3bar*y3bar
Rbar = dsqrt(R2bar)
R3bar = R2bar*Rbar
R4bar = R2bar*R2bar
R5bar = R3bar*R2bar
F = -datan2(y2,y1) + datan2(y2,z1) + datan2(y2*R*sinb,y1*z1+(y2*y2)*cosb)
Fbar = -datan2(y2,y1) + datan2(y2,z1bar) + datan2(y2*Rbar*sinb,y1*z1bar+(y2*y2)*cosb)
return
end subroutine

!--------------------------------------------------------------------------------------------------!

!**function d = LinePlaneIntersect(x, y, z, sx, sy, sz)
subroutine verticalLinePlaneIntersect(depthOfPlaneAtPoint,pointOfInterest,origin,normal)

!**% Calculate the intersection of a line and a plane using a parametric
!**% representation of the plane.  This is hardcoded for a vertical line.
!**numerator                       = [1 1 1 1 ; x(1) x(2) x(3) sx ; y(1) y(2) y(3) sy ; z(1) z(2) z(3) sz];
!**numerator                       = det(numerator);
!**denominator                     = [1 1 1 0 ; x(1) x(2) x(3) 0 ; y(1) y(2) y(3) 0 ; z(1) z(2) z(3) -sz];
!**denominator                     = det(denominator);
!**if denominator == 0;
!**   denominator                  = eps;
!**end
!**t                               = numerator/denominator; % parametric curve parameter
!**d                               = [sx sy sz]-([sx sy 0]-[sx sy sz])*t;

implicit none
! Arguments
double precision :: depthOfPlaneAtPoint, pointOfInterest(3), origin(3), normal(3)
! Local
double precision :: d

! Compute parameterization of plane (nx*x+ny*y+nz*z+d=0)
d = normal(1)*origin(1) + normal(2)*origin(2) + normal(3)*origin(3)
d = -d
!write(0,*) 'd:',d

! Compute depth of plane at point of interest
depthOfPlaneAtPoint = (-d-normal(1)*pointOfInterest(1)-normal(2)*pointOfInterest(2))/normal(3)
!write(0,*) 'z:',depthOfPlaneAtPoint

return
end subroutine

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
! COMPONENTS OF THE DISPLACEMENT VECTOR FOR AN ANGULAR DISLOCATION....EDIT AT YOUR OWN RISK!       !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

function v1B1()
implicit none
double precision :: v1B1, v1InfB1, v1CB1
v1InfB1 = 2.0d0*(1.0d0-nu)*(F+Fbar) - &
              y1*y2*(1.0d0/(R*(R-y3)) + &
              1.0d0/(Rbar*(Rbar+y3bar))) - &
              y2*cosb*((R*sinb-y1)/(R*(R-z3)) + &
              (Rbar*sinb-y1)/(Rbar*(Rbar+z3bar)))

v1InfB1 = v1InfB1/(8.0d0*pi*(1.0d0-nu))

v1CB1 = -2.0d0*(1.0d0-nu)*(1.0d0-2.0d0*nu)*Fbar*(cotb*cotb) + &
            (1.0d0-2.0d0*nu)*y2/(Rbar+y3bar)*((1.0d0-2.0d0*nu-a/Rbar)*cotb - &
            y1/(Rbar+y3bar)*(nu+a/Rbar)) + &
            (1.0d0-2.0d0*nu)*y2*cosb*cotb/(Rbar+z3bar)*(cosb+a/Rbar) + &
            a*y2*(y3bar-a)*cotb/(Rbar*Rbar*Rbar) + &
            y2*(y3bar-a)/(Rbar*(Rbar+y3bar))*(-(1.0d0-2.0d0*nu)*cotb + &
            y1/(Rbar+y3bar)*(2.0d0*nu+a/Rbar) + &
            a*y1/(Rbar*Rbar)) + &
            y2*(y3bar-a)/(Rbar*(Rbar+z3bar))*( &
                cosb/(Rbar+z3bar)*((Rbar*cosb+y3bar)*&
                          ((1.0d0-2.0d0*nu)*cosb-a/Rbar)*cotb + &
                2.0d0*(1.0d0-nu)*(Rbar*sinb-y1)*cosb) - &
                a*y3bar*cosb*cotb/(Rbar*Rbar))

v1CB1 = v1CB1/(4.0d0*pi*(1.0d0-nu))

v1B1 = v1InfB1 + v1CB1
return
end function

function v2B1()
implicit none
double precision :: v2B1, v2InfB1, v2CB1
v2InfB1 = (1.0d0-2.0d0*nu)*(dlog(R-y3)+dlog(Rbar+y3bar) - &
              cosb*(dlog(R-z3)+dlog(Rbar+z3bar))) - &
              y2*y2*(1.0d0/(R*(R-y3))+1.0d0/(Rbar*(Rbar+y3bar)) - &
              cosb*(1.0d0/(R*(R-z3))+1.0d0/(Rbar*(Rbar+z3bar))))

v2InfB1 = v2InfB1/(8.0d0*pi*(1.0d0-nu))

v2CB1 = (1.0d0-2.0d0*nu)*((2.0d0*(1.0d0-nu)*(cotb*cotb)-nu)*dlog(Rbar+y3bar) - &
            (2.0d0*(1.0d0-nu)*(cotb*cotb)+1.0d0-2.0d0*nu)*cosb*dlog(Rbar+z3bar)) - &
            (1.0d0-2.0d0*nu)/(Rbar+y3bar)*(y1*cotb*(1.0d0-2.0d0*nu-a/Rbar) + &
            nu*y3bar - a + &
            (y2*y2)/(Rbar+y3bar)*(nu+a/Rbar)) - &
            (1.0d0-2.0d0*nu)*z1bar*cotb/(Rbar+z3bar)*(cosb+a/Rbar) - &
            a*y1*(y3bar-a)*cotb/(Rbar*Rbar*Rbar) + &
            (y3bar-a)/(Rbar+y3bar)*(-2.0d0*nu + 1.0d0/Rbar*((1.0d0-2.0d0*nu)*y1*cotb-a) + &
            (y2*y2)/(Rbar*(Rbar+y3bar))*(2.0d0*nu+a/Rbar)+a*(y2*y2)/(Rbar*Rbar*Rbar)) + &
            (y3bar-a)/(Rbar+z3bar)*((cosb*cosb) - &
            1.0d0/Rbar*((1.0d0-2.0d0*nu)*z1bar*cotb+a*cosb) + &
            a*y3bar*z1bar*cotb/(Rbar*Rbar*Rbar) - &
            1.0d0/(Rbar*(Rbar+z3bar))*((y2*y2)*(cosb*cosb)-&
                   a*z1bar*cotb/Rbar*(Rbar*cosb+y3bar)))

v2CB1 = v2CB1/(4.0d0*pi*(1.0d0-nu))

v2B1 = v2InfB1 + v2CB1
return
end function

function v3B1()
implicit none
double precision :: v3B1, v3InfB1, v3CB1
v3InfB1 = y2*(1.0d0/R - &
              1.0d0/Rbar - &
              cosb*((R*cosb-y3)/(R*(R-z3)) - (Rbar*cosb+y3bar)/(Rbar*(Rbar+z3bar))))
v3InfB1 = v3InfB1/(8.0d0*pi*(1.0d0-nu))

v3CB1 = 2.0d0*(1.0d0-nu)*(((1.0d0-2.0d0*nu)*Fbar*cotb) + &
            (y2/(Rbar+y3bar)*(2.0d0*nu+a/Rbar)) - &
            (y2*cosb/(Rbar+z3bar)*(cosb+a/Rbar))) + &
            y2*(y3bar-a)/Rbar*(2.0d0*nu/(Rbar+y3bar)+a/(Rbar*Rbar)) + &
            y2*(y3bar-a)*cosb/(Rbar*(Rbar+z3bar))*(1.0d0-2.0d0*nu-&
                (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar)-a*y3bar/(Rbar*Rbar))
v3CB1 = v3CB1/(4.0d0*pi*(1.0d0-nu))

v3B1 = v3InfB1 + v3CB1
return
end function

function v1B2()
implicit none
double precision :: v1B2, v1InfB2, v1CB2
v1InfB2 = -(1.0d0-2.0d0*nu)*(dlog(R-y3)+dlog(Rbar+y3bar)-cosb*(dlog(R-z3)+dlog(Rbar+z3bar))) + &
              y1*y1*(1.0d0/(R*(R-y3))+1.0d0/(Rbar*(Rbar+y3bar))) + &
              z1*(R*sinb-y1)/(R*(R-z3)) + &
              z1bar*(Rbar*sinb-y1)/(Rbar*(Rbar+z3bar))

v1InfB2 = v1InfB2/(8.0d0*pi*(1.0d0-nu))

v1CB2 = (1.0d0-2.0d0*nu)*((2.0d0*(1.0d0-nu)*(cotb*cotb)+nu)*dlog(Rbar+y3bar) - &
            (2.0d0*(1.0d0-nu)*(cotb*cotb)+1.0d0)*cosb*dlog(Rbar+z3bar)) + &
            (1.0d0-2.0d0*nu)/(Rbar+y3bar)*(-(1.0d0-2.0d0*nu)*y1*cotb+nu*y3bar-&
                                                 a+a*y1*cotb/Rbar + &
            (y1*y1)/(Rbar+y3bar)*(nu+a/Rbar)) - &
            (1.0d0-2.0d0*nu)*cotb/(Rbar+z3bar)*(z1bar*cosb - &
            a*(Rbar*sinb-y1)/(Rbar*cosb)) - &
            a*y1*(y3bar-a)*cotb/(Rbar*Rbar*Rbar) + &
            (y3bar-a)/(Rbar+y3bar)*(2.0d0*nu+1.0d0/Rbar*((1.0d0-2.0d0*nu)*y1*cotb+a) - &
            (y1*y1)/(Rbar*(Rbar+y3bar))*(2.0d0*nu+a/Rbar) - a*(y1*y1)/(Rbar*Rbar*Rbar)) + &
            (y3bar-a)*cotb/(Rbar+z3bar)*(-cosb*sinb+a*y1*y3bar/(Rbar*Rbar*Rbar*cosb) + &
            (Rbar*sinb-y1)/Rbar*(2.0d0*(1.0d0-nu)*cosb-(Rbar*cosb+y3bar)/(Rbar+z3bar)*&
                (1.0d0+a/(Rbar*cosb))))

v1CB2 = v1CB2/(4.0d0*pi*(1.0d0-nu))

v1B2 = v1InfB2 + v1CB2
return
end function

function v2B2()
implicit none
double precision :: v2B2, v2InfB2, v2CB2
v2InfB2 = 2.0d0*(1.0d0-nu)*(F+Fbar) + &
              y1*y2*(1.0d0/(R*(R-y3))+1.0d0/(Rbar*(Rbar+y3bar))) - &
              y2*(z1/(R*(R-z3))+z1bar/(Rbar*(Rbar+z3bar)))

v2InfB2 = v2InfB2/(8.0d0*pi*(1.0d0-nu))

v2CB2 = 2.0d0*(1.0d0-nu)*(1.0d0-2.0d0*nu)*Fbar*cotb*cotb + &
            (1.0d0-2.0d0*nu)*y2/(Rbar+y3bar)*(-(1.0d0-2.0d0*nu-a/Rbar)*cotb + &
            y1/(Rbar+y3bar)*(nu+a/Rbar)) - &
            (1.0d0-2.0d0*nu)*y2*cotb/(Rbar+z3bar)*(1.0d0+a/(Rbar*cosb)) - &
            a*y2*(y3bar-a)*cotb/(Rbar*Rbar*Rbar) + &
            y2*(y3bar-a)/(Rbar*(Rbar+y3bar))*((1.0d0-2.0d0*nu)*cotb - &
            2.0d0*nu*y1/(Rbar+y3bar) - &
            a*y1/Rbar*(1.0d0/Rbar+1.0d0/(Rbar+y3bar))) + &
            y2*(y3bar-a)*cotb/(Rbar*(Rbar+z3bar))*(-2.0d0*(1.0d0-nu)*cosb + &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/(Rbar*cosb)) + &
            a*y3bar/((Rbar*Rbar)*cosb))

v2CB2 = v2CB2/(4.0d0*pi*(1.0d0-nu))

v2B2 = v2InfB2 + v2CB2
return
end function

function v3B2()
implicit none
double precision :: v3B2, v3InfB2, v3CB2
v3InfB2 = -(1.0d0-2.0d0*nu)*sinb*(dlog(R-z3)-dlog(Rbar+z3bar)) - &
              y1*(1.0d0/R-1.0d0/Rbar) + &
              z1*(R*cosb-y3)/(R*(R-z3)) - z1bar*(Rbar*cosb+y3bar)/(Rbar*(Rbar+z3bar))
v3InfB2 = v3InfB2/(8.0d0*pi*(1.0d0-nu))

v3CB2 = -2.0d0*(1.0d0-nu)*(1.0d0-2.0d0*nu)*cotb*&
            (dlog(Rbar+y3bar)-cosb*dlog(Rbar+z3bar)) - &
            2.0d0*(1.0d0-nu)*y1/(Rbar+y3bar)*(2.0d0*nu+a/Rbar) + &
            2.0d0*(1.0d0-nu)*z1bar/(Rbar+z3bar)*(cosb+a/Rbar) + &
            (y3bar-a)/Rbar*((1.0d0-2.0d0*nu)*cotb-2.0d0*nu*y1/(Rbar+y3bar)-a*y1/(Rbar*Rbar)) - &
            (y3bar-a)/(Rbar+z3bar)*(cosb*sinb + &
            (Rbar*cosb+y3bar)*cotb/Rbar*(2.0d0*(1.0d0-nu)*cosb - &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)) + &
            a/Rbar*(sinb - y3bar*z1bar/(Rbar*Rbar) - &
            z1bar*(Rbar*cosb+y3bar)/(Rbar*(Rbar+z3bar))))
v3CB2 = v3CB2/(4.0d0*pi*(1.0d0-nu))

v3B2 = v3InfB2 + v3CB2
return
end function

function v1B3()
implicit none
double precision :: v1B3, v1InfB3, v1CB3
v1InfB3 = y2*sinb*((R*sinb-y1)/(R*(R-z3))+(Rbar*sinb-y1)/(Rbar*(Rbar+z3bar)))

v1InfB3 = v1InfB3/(8.0d0*pi*(1.0d0-nu))

v1CB3 = (1.0d0-2.0d0*nu)*(y2/(Rbar+y3bar)*(1.0d0+a/Rbar) - &
            y2*cosb/(Rbar+z3bar)*(cosb+a/Rbar)) - &
            y2*(y3bar-a)/Rbar*(a/(Rbar*Rbar) + &
            1.0d0/(Rbar+y3bar)) + &
            y2*(y3bar-a)*cosb/(Rbar*(Rbar+z3bar))*((Rbar*cosb+y3bar)/&
                                    (Rbar+z3bar)*(cosb+a/Rbar) + a*y3bar/(Rbar*Rbar))

v1CB3 = v1CB3/(4.0d0*pi*(1.0d0-nu))

v1B3 = v1InfB3 + v1CB3
return
end function

function v2B3()
implicit none
double precision :: v2B3, v2InfB3, v2CB3
v2InfB3 = (1.0d0-2.0d0*nu)*sinb*(dlog(R-z3)+dlog(Rbar+z3bar)) - &
              (y2*y2)*sinb*(1.0d0/(R*(R-z3))+1.0d0/(Rbar*(Rbar+z3bar)))

v2InfB3 = v2InfB3/(8.0d0*pi*(1.0d0-nu))

v2CB3 = (1.0d0-2.0d0*nu)*(-sinb*dlog(Rbar+z3bar) - &
            y1/(Rbar+y3bar)*(1.0d0+a/Rbar) + &
            z1bar/(Rbar+z3bar)*(cosb+a/Rbar)) + &
            y1*(y3bar-a)/Rbar*(a/(Rbar*Rbar) + &
            1.0d0/(Rbar+y3bar)) - &
            (y3bar-a)/(Rbar+z3bar)*(sinb*(cosb-a/Rbar) + &
            z1bar/Rbar*(1.0d0+a*y3bar/(Rbar*Rbar)) - &
            1.0d0/(Rbar*(Rbar+z3bar))*((y2*y2)*cosb*sinb - &
            a*z1bar/Rbar*(Rbar*cosb+y3bar)))

v2CB3 = v2CB3/(4.0d0*pi*(1.0d0-nu))

v2B3 = v2InfB3 + v2CB3
return
end function

function v3B3()
implicit none
double precision :: v3B3, v3InfB3, v3CB3
v3InfB3 = 2.0d0*(1.0d0-nu)*(F-Fbar) + &
              y2*sinb*((R*cosb-y3)/(R*(R-z3))-(Rbar*cosb+y3bar)/(Rbar*(Rbar+z3bar)))
v3InfB3 = v3InfB3/(8.0d0*pi*(1.0d0-nu))

v3CB3 = 2.0d0*(1.0d0-nu)*Fbar + &
            2.0d0*(1.0d0-nu)*(y2*sinb/(Rbar+z3bar)*(cosb + a/Rbar)) + &
            y2*(y3bar-a)*sinb/(Rbar*(Rbar+z3bar))*(1.0d0+(Rbar*cosb+y3bar)/&
                (Rbar+z3bar)*(cosb+a/Rbar) + a*y3bar/(Rbar*Rbar))
v3CB3 = v3CB3/(4.0d0*pi*(1.0d0-nu))

v3B3 = v3InfB3 + v3CB3
return
end function


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
! COMPONENTS OF THE STRAIN TENSOR FOR AN ANGULAR DISLOCATION....EDIT AT YOUR OWN RISK!             !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

function e11_termB1a()
implicit none
double precision :: e11_termB1a
e11_termB1a = 1.0d0/8.0d0*( &
    (2.0d0-2.0d0*nu)*( &
        2.0d0*y2/y1**2/(1.0d0+y2**2/y1**2) - &
        y2/z1**2*cosb/(1.0d0+y2**2/z1**2) + &
        ( &
            y2/R*sinb/(y1*z1+y2**2*cosb)*y1 - &
            y2*R*sinb/(y1*z1+y2**2*cosb)**2*(2*y1*cosb-y3*sinb) &
        )/(1.0d0+y2**2*R2*sinb**2/(y1*z1+y2**2*cosb)**2) - &
        y2/z1bar**2*cosb/(1.0d0+y2**2/z1bar**2) + &
        ( &
            y2/Rbar*sinb/(y1*z1bar+y2**2*cosb)*y1 - &
            y2*Rbar*sinb/(y1*z1bar+y2**2*cosb)**2*(2.0d0*y1*cosb+y3bar*sinb) &
        )/(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
    ) - &
    y2*(1.0d0/R/(R-y3)+1.0d0/Rbar/(Rbar+y3bar)) - &
    y1*y2*( &
        -1.0d0/R3/(R-y3)*y1 - &
        1.0d0/R2/(R-y3)**2*y1 - &
        1.0d0/R3bar/(Rbar+y3bar)*y1 - &
        1.0d0/R2bar/(Rbar+y3bar)**2*y1 &
    ) - &
    y2*cosb*( &
        (1.0d0/R*sinb*y1-1.0d0)/R/(R-z3) - &
        (R*sinb-y1)/R3/(R-z3)*y1 - &
        (R*sinb-y1)/R/(R-z3)**2*(1.0d0/R*y1-sinb) + &
        (1.0d0/Rbar*sinb*y1-1.0d0)/Rbar/(Rbar+z3bar) - &
        (Rbar*sinb-y1)/R3bar/(Rbar+z3bar)*y1 - &
        (Rbar*sinb-y1)/Rbar/(Rbar+z3bar)**2*(1.0d0/Rbar*y1-sinb) &
    ) &
)/pi/(1.0d0-nu)
return
end function

function e11_termB1b()
implicit none
double precision :: e11_termB1b
e11_termB1b = 1.0d0/4.0d0*( &
    (-2.0d0+2.0d0*nu)*(1.0d0-2.0d0*nu)*( &
        y2/y1**2/(1.0d0+y2**2/y1**2) - &
        y2/z1bar**2*cosb/(1.0d0+y2**2/z1bar**2) + &
        (y2/Rbar*sinb/(y1*z1bar+y2**2*cosb)*y1 - &
        y2*Rbar*sinb/(y1*z1bar+y2**2*cosb)**2*(2.0d0*y1*cosb+y3bar*sinb)) &
         /(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
    )*cotb**2 - &
    (1.0d0-2.0d0*nu)*y2/(Rbar+y3bar)**2*( &
        (1.0d0-2.0d0*nu-a/Rbar)*cotb-y1/(Rbar+y3bar)*(nu+a/Rbar) &
    )/Rbar*y1 + &
    (1.0d0-2.0d0*nu)*y2/(Rbar+y3bar)*( &
        a/R3bar*y1*cotb - &
        1.0d0/(Rbar+y3bar)*(nu+a/Rbar) + &
        y1**2/(Rbar+y3bar)**2*(nu+a/Rbar)/Rbar + &
        y1**2/(Rbar+y3bar)*a/R3bar &
    ) - &
    (1.0d0-2.0d0*nu)*y2*cosb*cotb/(Rbar+z3bar)**2*(cosb+a/Rbar)*(1.0d0/Rbar*y1-sinb) - &
    (1.0d0-2.0d0*nu)*y2*cosb*cotb/(Rbar+z3bar)*a/R3bar*y1 - &
    3.0d0*a*y2*(y3+a)*cotb/R5bar*y1 - &
    y2*(y3+a)/R3bar/(Rbar+y3bar)*( &
        (-1.0d0+2.0d0*nu)*cotb+y1/(Rbar+y3bar)*(2.0d0*nu+a/Rbar)+a*y1/R2bar &
    )*y1 - &
    y2*(y3+a)/R2bar/(Rbar+y3bar)**2*( &
        (-1.0d0+2.0d0*nu)*cotb+y1/(Rbar+y3bar)*(2.0d0*nu+a/Rbar)+a*y1/R2bar &
    )*y1 + &
    y2*(y3+a)/Rbar/(Rbar+y3bar)*( &
        1.0d0/(Rbar+y3bar)*(2.0d0*nu+a/Rbar) - &
        y1**2/(Rbar+y3bar)**2*(2.0d0*nu+a/Rbar)/Rbar - &
        y1**2/(Rbar+y3bar)*a/R3bar+a/R2bar - &
        2.0d0*a*y1**2/R4bar &
    ) - &
    y2*(y3+a)/R3bar/(Rbar+z3bar)*( &
        cosb/(Rbar+z3bar)*( &
            (Rbar*cosb+y3bar)*((1.0d0-2.0d0*nu)*cosb-a/Rbar)*cotb + &
            (2.0d0-2.0d0*nu)*(Rbar*sinb-y1)*cosb &
        ) - &
        a*y3bar*cosb*cotb/R2bar &
    )*y1 - &
    y2*(y3+a)/Rbar/(Rbar+z3bar)**2*( &
        cosb/(Rbar+z3bar)*( &
            (Rbar*cosb+y3bar)*((1.0d0-2.0d0*nu)*cosb-a/Rbar)*cotb + &
            (2.0d0-2.0d0*nu)*(Rbar*sinb-y1)*cosb &
        ) - &
        a*y3bar*cosb*cotb/R2bar &
    )*(1.0d0/Rbar*y1-sinb) + &
    y2*(y3+a)/Rbar/(Rbar+z3bar)*( &
        -cosb/(Rbar+z3bar)**2*( &
            (Rbar*cosb+y3bar)*((1.0d0-2.0d0*nu)*cosb-a/Rbar)*cotb + &
            (2.0d0-2.0d0*nu)*(Rbar*sinb-y1)*cosb &
        )*(1.0d0/Rbar*y1-sinb) + &
        cosb/(Rbar+z3bar)*( &
            1.0d0/Rbar*cosb*y1*((1.0d0-2.0d0*nu)*cosb-a/Rbar)*cotb + &
            (Rbar*cosb+y3bar)*a/R3bar*y1*cotb + &
            (2.0d0-2.0d0*nu)*(1.0d0/Rbar*sinb*y1-1.0d0)*cosb &
        )+2.0d0*a*y3bar*cosb*cotb/R4bar*y1 &
    ) &
)/pi/(1.0d0-nu)
return
end function

function e11_termB2a()
implicit none
double precision :: e11_termB2a
e11_termB2a = 1.0d0/8.0d0*( &
    (-1.0d0+2.0d0*nu)*( &
        1.0d0/R*y1/(R-y3) + &
        1.0d0/Rbar*y1/(Rbar+y3bar) - &
        cosb*((1.0d0/R*y1-sinb)/(R-z3)+(1.0d0/Rbar*y1-sinb)/(Rbar+z3bar)) &
    ) + &
    2.0d0*y1*(1.0d0/R/(R-y3)+1.0d0/Rbar/(Rbar+y3bar)) + &
    y1**2*( &
        -1.0d0/R3/(R-y3)*y1 - &
        1.0d0/R2/(R-y3)**2*y1 - &
        1.0d0/R3bar/(Rbar+y3bar)*y1 - &
        1.0d0/(y1**2+y2**2+y3bar**2)/(Rbar+y3bar)**2*y1 &
    ) + &
    cosb*(R*sinb-y1)/R/(R-z3) + &
    z1*(1.0d0/R*sinb*y1-1.0d0)/R/(R-z3) - &
    z1*(R*sinb-y1)/R3/(R-z3)*y1 - &
    z1*(R*sinb-y1)/R/(R-z3)**2*(1.0d0/R*y1-sinb) + &
    cosb*(Rbar*sinb-y1)/Rbar/(Rbar+z3bar) + &
    z1bar*(1.0d0/Rbar*sinb*y1-1.0d0)/Rbar/(Rbar+z3bar) - &
    z1bar*(Rbar*sinb-y1)/R3bar/(Rbar+z3bar)*y1 - &
    z1bar*(Rbar*sinb-y1)/Rbar/(Rbar+z3bar)**2*(1.0d0/Rbar*y1-sinb) &
)/pi/(1.0d0-nu)
return
end function

function e11_termB2b()
implicit none
double precision :: e11_termB2b
e11_termB2b = 1.0d0/4.0d0*( &
    (1.0d0-2.0d0*nu)*( &
        ((2.0d0-2.0d0*nu)*cotb**2+nu)/Rbar*y1/(Rbar+y3bar) - &
        ((2.0d0-2.0d0*nu)*cotb**2+1.0d0)*cosb*(1.0d0/Rbar*y1-sinb)/(Rbar+z3bar) &
    ) - &
    (1.0d0-2.0d0*nu)/(Rbar+y3bar)**2*( &
        (-1.0d0+2.0d0*nu)*y1*cotb + &
        nu*y3bar - &
        a + &
        a*y1*cotb/Rbar + &
        y1**2/(Rbar+y3bar)*(nu+a/Rbar) &
    )/Rbar*y1 + &
    (1.0d0-2.0d0*nu)/(Rbar+y3bar)*( &
        (-1.0d0+2.0d0*nu)*cotb + &
        a*cotb/Rbar - &
        a*y1**2*cotb/R3bar + &
        2.0d0*y1/(Rbar+y3bar)*(nu+a/Rbar) - &
        y1**3/(Rbar+y3bar)**2*(nu+a/Rbar)/Rbar - &
        y1**3/(Rbar+y3bar)*a/R3bar &
     ) + &
    (1.0d0-2.0d0*nu)*cotb/(Rbar+z3bar)**2*( &
        z1bar*cosb - &
        a*(Rbar*sinb-y1)/Rbar/cosb &
    )*(1.0d0/Rbar*y1-sinb) - &
    (1.0d0-2.0d0*nu)*cotb/(Rbar+z3bar)*( &
        cosb**2 - &
        a*(1.0d0/Rbar*sinb*y1-1.0d0)/Rbar/cosb + &
        a*(Rbar*sinb-y1)/R3bar/cosb*y1 &
    ) - &
    a*(y3+a)*cotb/R3bar + &
    3.0d0*a*y1**2*(y3+a)*cotb/R5bar - &
    (y3+a)/(Rbar+y3bar)**2*( &
        2.0d0*nu + &
        1.0d0/Rbar*((1.0d0-2.0d0*nu)*y1*cotb+a) - &
        y1**2/Rbar/(Rbar+y3bar)*(2.0d0*nu+a/Rbar) - &
        a*y1**2/R3bar &
    )/Rbar*y1 + &
    (y3+a)/(Rbar+y3bar)*( &
        -1.0d0/R3bar*((1.0d0-2.0d0*nu)*y1*cotb+a)*y1 + &
        1.0d0/Rbar*(1.0d0-2.0d0*nu)*cotb - &
        2.0d0*y1/Rbar/(Rbar+y3bar)*(2.0d0*nu+a/Rbar) + &
        y1**3/R3bar/(Rbar+y3bar)*(2.0d0*nu+a/Rbar) + &
        y1**3/R2bar/(Rbar+y3bar)**2*(2.0d0*nu+a/Rbar) + &
        y1**3/R4bar/(Rbar+y3bar)*a - &
        2.0d0*a/R3bar*y1 + &
        3.0d0*a*y1**3/R5bar &
    ) - &
    (y3+a)*cotb/(Rbar+z3bar)**2*( &
        -cosb*sinb + &
        a*y1*y3bar/R3bar/cosb + &
        (Rbar*sinb-y1)/Rbar*( &
            (2.0d0-2.0d0*nu)*cosb - &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) &
        ) &
    )*(1.0d0/Rbar*y1-sinb) + &
    (y3+a)*cotb/(Rbar+z3bar)*( &
        a*y3bar/R3bar/cosb - &
        3.0d0*a*y1**2*y3bar/R5bar/cosb + &
        (1.0d0/Rbar*sinb*y1-1.0d0)/Rbar*( &
            (2.0d0-2.0d0*nu)*cosb - &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) &
        ) - &
        (Rbar*sinb-y1)/R3bar*( &
            (2.0d0-2.0d0*nu)*cosb - &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) &
        )*y1+(Rbar*sinb-y1)/Rbar*( &
            -1.0d0/Rbar*cosb*y1/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) + &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)**2*(1.0d0+a/Rbar/cosb)*(1.0d0/Rbar*y1-sinb) + &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)*a/R3bar/cosb*y1 &
        ) &
    ) &
)/pi/(1.0d0-nu)
return
end function

function e11_termB3a()
implicit none
double precision :: e11_termB3a
e11_termB3a = 1.0d0/8.0d0*y2*sinb*( &
    (1.0d0/R*sinb*y1-1.0d0)/R/(R-z3) - &
    (R*sinb-y1)/R3/(R-z3)*y1 - &
    (R*sinb-y1)/R/(R-z3)**2*(1.0d0/R*y1-sinb) + &
    (1.0d0/Rbar*sinb*y1-1.0d0)/Rbar/(Rbar+z3bar) - &
    (Rbar*sinb-y1)/R3bar/(Rbar+z3bar)*y1 - &
    (Rbar*sinb-y1)/Rbar/(Rbar+z3bar)**2*(1.0d0/Rbar*y1-sinb) &
)/pi/(1.0d0-nu)
return
end function

function e11_termB3b()
implicit none
double precision :: e11_termB3b
e11_termB3b = 1.0d0/4.0d0*( &
    (1.0d0-2.0d0*nu)*( &
        -y2/(Rbar+y3bar)**2*(1.0d0+a/Rbar)/Rbar*y1 - &
        y2/(Rbar+y3bar)*a/R3bar*y1 + &
        y2*cosb/(Rbar+z3bar)**2*(cosb+a/Rbar)*(1.0d0/Rbar*y1-sinb) + &
        y2*cosb/(Rbar+z3bar)*a/R3bar*y1 &
    ) + &
    y2*(y3+a)/R3bar*(a/(y1**2+y2**2+y3bar**2)+1/(Rbar+y3bar))*y1 - &
    y2*(y3+a)/Rbar*(-2.0d0*a/R4bar*y1-1.0d0/(Rbar+y3bar)**2/Rbar*y1) - &
    y2*(y3+a)*cosb/R3bar/(Rbar+z3bar)*( &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) + &
        a*y3bar/(y1**2+y2**2+y3bar**2) &
    )*y1 - &
    y2*(y3+a)*cosb/Rbar/(Rbar+z3bar)**2*( &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) + &
        a*y3bar/(y1**2+y2**2+y3bar**2) &
    )*(1.0d0/Rbar*y1-sinb) + &
    y2*(y3+a)*cosb/Rbar/(Rbar+z3bar)*( &
        1.0d0/Rbar*cosb*y1/(Rbar+z3bar)*(cosb+a/Rbar) - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)**2*(cosb+a/Rbar)*(1.0d0/Rbar*y1-sinb) - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*a/R3bar*y1 - &
        2.0d0*a*y3bar/R4bar*y1 &
    ) &
)/pi/(1.0d0-nu)
return
end function

function e22_termB1a()
implicit none
double precision :: e22_termB1a
e22_termB1a = 1.0d0/8.0d0*( &
    (1.0d0-2.0d0*nu)*( &
        1.0d0/R*y2/(R-y3) + &
        1.0d0/Rbar*y2/(Rbar+y3bar) - &
        cosb*(1.0d0/R*y2/(R-z3)+1.0d0/Rbar*y2/(Rbar+z3bar)) &
    ) - &
    2.0d0*y2*( &
        1.0d0/R/(R-y3) + &
        1.0d0/Rbar/(Rbar+y3bar) - &
        cosb*(1.0d0/R/(R-z3)+1.0d0/Rbar/(Rbar+z3bar)) &
    ) - &
    y2**2*( &
        -1.0d0/R3/(R-y3)*y2 - &
        1.0d0/R2/(R-y3)**2*y2 - &
        1.0d0/R3bar/(Rbar+y3bar)*y2 - &
        1.0d0/R2bar/(Rbar+y3bar)**2*y2 - &
        cosb*( &
            -1.0d0/R3/(R-z3)*y2 - &
            1.0d0/R2/(R-z3)**2*y2 - &
            1.0d0/R3bar/(Rbar+z3bar)*y2 - &
            1.0d0/R2bar/(Rbar+z3bar)**2*y2 &
        ) &
    ) &
)/pi/(1.0d0-nu)
return
end function

function e22_termB1b()
implicit none
double precision :: e22_termB1b
e22_termB1b = 1.0d0/4.0d0*( &
    (1.0d0-2.0d0*nu)*( &
        ((2.0d0-2.0d0*nu)*cotb**2-nu)/Rbar*y2/(Rbar+y3bar) - &
        ((2.0d0-2.0d0*nu)*cotb**2+1-2.0d0*nu)*cosb/Rbar*y2/(Rbar+z3bar) &
    ) + &
    (1.0d0-2.0d0*nu)/(Rbar+y3bar)**2*( &
        y1*cotb*(1.0d0-2.0d0*nu-a/Rbar)+nu*y3bar - &
        a+y2**2/(Rbar+y3bar)*(nu+a/Rbar) &
    )/Rbar*y2 - &
    (1.0d0-2.0d0*nu)/(Rbar+y3bar)*( &
        a*y1*cotb/R3bar*y2 + &
        2.0d0*y2/(Rbar+y3bar)*(nu+a/Rbar) - &
        y2**3/(Rbar+y3bar)**2*(nu+a/Rbar)/Rbar - &
        y2**3/(Rbar+y3bar)*a/R3bar &
    ) + &
    (1.0d0-2.0d0*nu)*z1bar*cotb/(Rbar+z3bar)**2*(cosb+a/Rbar)/Rbar*y2 + &
    (1.0d0-2.0d0*nu)*z1bar*cotb/(Rbar+z3bar)*a/R3bar*y2 + &
    3.0d0*a*y2*(y3+a)*cotb/R5bar*y1 - &
    (y3+a)/(Rbar+y3bar)**2*( &
        -2.0d0*nu+1/Rbar*((1.0d0-2.0d0*nu)*y1*cotb-a) + &
        y2**2/Rbar/(Rbar+y3bar)*(2*nu+a/Rbar) + &
        a*y2**2/R3bar &
    )/Rbar*y2 + &
    (y3+a)/(Rbar+y3bar)*( &
        -1/R3bar*((1.0d0-2.0d0*nu)*y1*cotb-a)*y2 + &
        2.0d0*y2/Rbar/(Rbar+y3bar)*(2*nu+a/Rbar) - &
        y2**3/R3bar/(Rbar+y3bar)*(2*nu+a/Rbar) - &
        y2**3/R2bar/(Rbar+y3bar)**2*(2*nu+a/Rbar) - &
        y2**3/R4bar/(Rbar+y3bar)*a + &
        2.0d0*a/R3bar*y2 - &
        3.0d0*a*y2**3/R5bar &
    ) - &
    (y3+a)/(Rbar+z3bar)**2*( &
        cosb**2 - &
        1.0d0/Rbar*((1.0d0-2.0d0*nu)*z1bar*cotb+a*cosb) + &
        a*y3bar*z1bar*cotb/R3bar - &
        1.0d0/Rbar/(Rbar+z3bar)*( &
            y2**2*cosb**2 - &
            a*z1bar*cotb/Rbar*(Rbar*cosb+y3bar) &
        ) &
    )/Rbar*y2 + &
    (y3+a)/(Rbar+z3bar)*( &
        1.0d0/R3bar*((1.0d0-2.0d0*nu)*z1bar*cotb+a*cosb)*y2 - &
        3.0d0*a*y3bar*z1bar*cotb/R5bar*y2 + &
        1.0d0/R3bar/(Rbar+z3bar)*( &
            y2**2*cosb**2 - &
            a*z1bar*cotb/Rbar*(Rbar*cosb+y3bar) &
        )*y2 + &
        1.0d0/R2bar/(Rbar+z3bar)**2*( &
            y2**2*cosb**2 - &
            a*z1bar*cotb/Rbar*(Rbar*cosb+y3bar) &
        )*y2 - &
        1.0d0/Rbar/(Rbar+z3bar)*( &
            2.0d0*y2*cosb**2 + &
            a*z1bar*cotb/R3bar*(Rbar*cosb+y3bar)*y2 - &
            a*z1bar*cotb/R2bar*cosb*y2 &
        ) &
    ) &
)/pi/(1.0d0-nu)
return
end function

function e22_termB2a()
implicit none
double precision :: e22_termB2a
e22_termB2a = 1.0d0/8.0d0*( &
    (2.0d0-2.0d0*nu)*( &
        -2/y1/(1.0d0+y2**2/y1**2) + &
        1.0d0/z1/(1.0d0+y2**2/z1**2) + &
        (R*sinb/(y1*z1+y2**2*cosb) + &
        y2**2/R*sinb/(y1*z1+y2**2*cosb) - &
        2.0d0*y2**2*R*sinb/(y1*z1+y2**2*cosb)**2*cosb)/( &
            1 + &
            y2**2*R2*sinb**2/(y1*z1+y2**2*cosb)**2 &
        ) + &
        1.0d0/z1bar/(1.0d0+y2**2/z1bar**2) + &
        (Rbar*sinb/(y1*z1bar+y2**2*cosb) + &
        y2**2/Rbar*sinb/(y1*z1bar+y2**2*cosb) - &
        2.0d0*y2**2*Rbar*sinb/(y1*z1bar+y2**2*cosb)**2*cosb)/( &
            1 + &
            y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2 &
        ) &
    ) + &
    y1*(1.0d0/R/(R-y3)+1/Rbar/(Rbar+y3bar)) + &
    y1*y2*( &
        -1/R3/(R-y3)*y2 - &
        1.0d0/R2/(R-y3)**2*y2 - &
        1.0d0/R3bar/(Rbar+y3bar)*y2 - &
        1.0d0/R2bar/(Rbar+y3bar)**2*y2 &
    ) - &
    z1/R/(R-z3) - &
    z1bar/Rbar/(Rbar+z3bar) - &
    y2*( &
        -z1/R3/(R-z3)*y2 - &
        z1/R2/(R-z3)**2*y2 - &
        z1bar/R3bar/(Rbar+z3bar)*y2 - &
        z1bar/R2bar/(Rbar+z3bar)**2*y2 &
    ) &
)/pi/(1.0d0-nu)
return
end function

function e22_termB2b()
implicit none
double precision :: e22_termB2b
e22_termB2b = 1.0/4.0*( &
    (2.0d0-2.0d0*nu)*(1.0d0-2.0d0*nu)*( &
        -1/y1/(1.0d0+y2**2/y1**2) + &
        1.0d0/z1bar/(1.0d0+y2**2/z1bar**2) + &
        ( &
            Rbar*sinb/(y1*z1bar+y2**2*cosb) + &
            y2**2/Rbar*sinb/(y1*z1bar+y2**2*cosb) - &
            2.0d0*y2**2*Rbar*sinb/(y1*z1bar+y2**2*cosb)**2*cosb &
        )/( &
            1 + &
            y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2 &
        ) &
    )*cotb**2 + &
    (1.0d0-2.0d0*nu)/(Rbar+y3bar)*( &
        (-1+2*nu+a/Rbar)*cotb + &
        y1/(Rbar+y3bar)*(nu+a/Rbar) &
    ) - &
    (1.0d0-2.0d0*nu)*y2**2/(Rbar+y3bar)**2*( &
        (-1+2*nu+a/Rbar)*cotb + &
        y1/(Rbar+y3bar)*(nu+a/Rbar))/Rbar + &
    (1.0d0-2.0d0*nu)*y2/(Rbar+y3bar)*( &
        -a/R3bar*y2*cotb - &
        y1/(Rbar+y3bar)**2*(nu+a/Rbar)/Rbar*y2 - &
        y2/(Rbar+y3bar)*a/R3bar*y1 &
    ) - &
    (1.0d0-2.0d0*nu)*cotb/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) + &
    (1.0d0-2.0d0*nu)*y2**2*cotb/(Rbar+z3bar)**2*(1.0d0+a/Rbar/cosb)/Rbar + &
    (1.0d0-2.0d0*nu)*y2**2*cotb/(Rbar+z3bar)*a/R3bar/cosb - &
    a*(y3+a)*cotb/R3bar + &
    3.0d0*a*y2**2*(y3+a)*cotb/R5bar + &
    (y3+a)/Rbar/(Rbar+y3bar)*( &
        (1.0d0-2.0d0*nu)*cotb - &
        2.0d0*nu*y1/(Rbar+y3bar) - &
        a*y1/Rbar*(1.0d0/Rbar+1/(Rbar+y3bar)) &
    ) - &
    y2**2*(y3+a)/R3bar/(Rbar+y3bar)*( &
        (1.0d0-2.0d0*nu)*cotb - &
        2.0d0*nu*y1/(Rbar+y3bar) - &
        a*y1/Rbar*(1.0d0/Rbar+1/(Rbar+y3bar)) &
    ) - &
    y2**2*(y3+a)/R2bar/(Rbar+y3bar)**2*( &
        (1.0d0-2.0d0*nu)*cotb - &
        2.0d0*nu*y1/(Rbar+y3bar) - &
        a*y1/Rbar*(1.0d0/Rbar+1/(Rbar+y3bar)) &
    ) + &
    y2*(y3+a)/Rbar/(Rbar+y3bar)*( &
        2.0d0*nu*y1/(Rbar+y3bar)**2/Rbar*y2 + &
        a*y1/R3bar*(1.0d0/Rbar+1/(Rbar+y3bar))*y2 - &
        a*y1/Rbar*(-1/R3bar*y2-1/(Rbar+y3bar)**2/Rbar*y2) &
    ) + &
    (y3+a)*cotb/Rbar/(Rbar+z3bar)*( &
        (-2.0d0+2*nu)*cosb + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) + &
        a*y3bar/R2bar/cosb &
    ) - &
    y2**2*(y3+a)*cotb/R3bar/(Rbar+z3bar)*( &
        (-2.0d0+2*nu)*cosb + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) + &
        a*y3bar/R2bar/cosb &
    ) - &
    y2**2*(y3+a)*cotb/R2bar/(Rbar+z3bar)**2*( &
        (-2.0d0+2*nu)*cosb + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) + &
        a*y3bar/R2bar/cosb &
    ) + &
    y2*(y3+a)*cotb/Rbar/(Rbar+z3bar)*( &
        1.0d0/Rbar*cosb*y2/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)**2*(1.0d0+a/Rbar/cosb)/Rbar*y2 - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*a/R3bar/cosb*y2 - &
        2.0d0*a*y3bar/R4bar/cosb*y2 &
    ) &
)/pi/(1.0d0-nu)
return
end function

function e22_termB3a()
implicit none
double precision :: e22_termB3a
e22_termB3a = 1.0d0/8.0d0*( &
    (1.0d0-2.0d0*nu)*sinb*(1.0d0/R*y2/(R-z3)+1/Rbar*y2/(Rbar+z3bar)) - &
    2.0d0*y2*sinb*(1.0d0/R/(R-z3)+1/Rbar/(Rbar+z3bar)) - &
    y2**2*sinb*( &
        -1/R3/(R-z3)*y2 - &
        1.0d0/R2/(R-z3)**2*y2 - &
        1.0d0/R3bar/(Rbar+z3bar)*y2 - &
    1.0d0/R2bar/(Rbar+z3bar)**2*y2 &
    ) &
)/pi/(1.0d0-nu)
return
end function

function e22_termB3b()
implicit none
double precision :: e22_termB3b
e22_termB3b = 1.0d0/4.0d0*( &
    (1.0d0-2.0d0*nu)*( &
        -sinb/Rbar*y2/(Rbar+z3bar) + &
        y2/(Rbar+y3bar)**2*(1.0d0+a/Rbar)/Rbar*y1 + &
        y2/(Rbar+y3bar)*a/R3bar*y1 - &
        z1bar/(Rbar+z3bar)**2*(cosb+a/Rbar)/Rbar*y2 - &
        z1bar/(Rbar+z3bar)*a/R3bar*y2 &
    ) - &
    y2*(y3+a)/R3bar*(a/R2bar+1/(Rbar+y3bar))*y1 + &
    y1*(y3+a)/Rbar*(-2.0d0*a/R4bar*y2-1/(Rbar+y3bar)**2/Rbar*y2) + &
    (y3+a)/(Rbar+z3bar)**2*( &
        sinb*(cosb-a/Rbar) + &
        z1bar/Rbar*(1.0d0+a*y3bar/R2bar) - &
        1.0d0/Rbar/(Rbar+z3bar)*(y2**2*cosb*sinb-a*z1bar/Rbar*(Rbar*cosb+y3bar)) &
    )/Rbar*y2 - &
    (y3+a)/(Rbar+z3bar)*( &
        sinb*a/R3bar*y2 - &
        z1bar/R3bar*(1.0d0+a*y3bar/R2bar)*y2 - &
        2.0d0*z1bar/R5bar*a*y3bar*y2 + &
        1.0d0/R3bar/(Rbar+z3bar)*(y2**2*cosb*sinb-a*z1bar/Rbar*(Rbar*cosb+y3bar))*y2 + &
        1.0d0/R2bar/(Rbar+z3bar)**2*( &
            y2**2*cosb*sinb - &
            a*z1bar/Rbar*(Rbar*cosb+y3bar) &
        )*y2 - &
        1.0d0/Rbar/(Rbar+z3bar)*( &
            2.0d0*y2*cosb*sinb + &
            a*z1bar/R3bar*(Rbar*cosb+y3bar)*y2 - &
            a*z1bar/R2bar*cosb*y2 &
        ) &
    ) &
)/pi/(1.0d0-nu)
return
end function

function e33_termB1a()
implicit none
double precision :: e33_termB1a
e33_termB1a = 1.0d0/8.0d0*y2*( &
    -1.0d0/R3*y3+1.0d0/2.0d0/R3bar*(2.0d0*y3+4.0d0*a) - &
    cosb*( &
        (1.0d0/R*cosb*y3-1.0d0)/R/(R-z3) - &
        (R*cosb-y3)/R3/(R-z3)*y3 - &
        (R*cosb-y3)/R/(R-z3)**2*(1.0d0/R*y3-cosb) - &
        (1.0d0/2.0d0/Rbar*cosb*(2.0d0*y3+4.0d0*a)+1.0d0)/Rbar/(Rbar+z3bar) + &
        1.0d0/2.0d0*(Rbar*cosb+y3bar)/R3bar/(Rbar+z3bar)*(2.0d0*y3+4.0d0*a) + &
        (Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) &
    ) &
)/pi/(1.0d0-nu)
return
end function

function e33_termB1b()
implicit none
double precision :: e33_termB1b
e33_termB1b = 1.0d0/4.0d0*( &
    (2.0d0-2.0d0*nu)*( &
        (1.0d0-2.0d0*nu)*( &
            -y2/z1bar**2*sinb/(1.0d0+y2**2/z1bar**2) + &
            ( &
                1.0d0/2.0d0*y2/Rbar*sinb/(y1*z1bar+y2**2*cosb)*(2.0d0*y3+4.0d0*a) - &
                y2*Rbar*sinb**2/(y1*z1bar+y2**2*cosb)**2*y1 &
             )/(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
        )*cotb - &
        y2/(Rbar+y3bar)**2*(2.0d0*nu+a/Rbar)*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) - &
        1.0d0/2.0d0*y2/(Rbar+y3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) + &
        y2*cosb/(Rbar+z3bar)**2*(cosb+a/Rbar)*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) + &
        1.0d0/2.0d0*y2*cosb/(Rbar+z3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) &
    ) + &
    y2/Rbar*(2.0d0*nu/(Rbar+y3bar)+a/R2bar) - &
    1.0d0/2.0d0*y2*(y3+a)/R3bar*(2.0d0*nu/(Rbar+y3bar)+a/R2bar)*(2.0d0*y3+4.0d0*a) + &
    y2*(y3+a)/Rbar*( &
        -2.0d0*nu/(Rbar+y3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) - &
        a/R4bar*(2.0d0*y3+4.0d0*a) &
    ) + &
    y2*cosb/Rbar/(Rbar+z3bar)*( &
        1.0d0-2.0d0*nu-(Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar)-a*y3bar/R2bar &
    ) - &
    1.0d0/2.0d0*y2*(y3+a)*cosb/R3bar/(Rbar+z3bar)*( &
        1.0d0-2.0d0*nu-(Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar)-a*y3bar/R2bar &
    )*(2.0d0*y3+4.0d0*a) - &
    y2*(y3+a)*cosb/Rbar/(Rbar+z3bar)**2*( &
        1.0d0-2.0d0*nu-(Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar)-a*y3bar/R2bar &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) + &
    y2*(y3+a)*cosb/Rbar/(Rbar+z3bar)*( &
        -(1.0d0/2.0d0/Rbar*cosb*(2.0d0*y3+4.0d0*a)+1.0d0)/(Rbar+z3bar)*(cosb+a/Rbar) + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)**2*(cosb+a/Rbar)*( &
            1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb &
        ) + &
        1.0d0/2.0d0*(Rbar*cosb+y3bar)/(Rbar+z3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) - &
        a/R2bar + &
        a*y3bar/R4bar*(2.0d0*y3+4.0d0*a) &
    ) &
)/pi/(1.0d0-nu)
return
end function

function e33_termB2a()
implicit none
double precision :: e33_termB2a
e33_termB2a = 1.0d0/8.0d0*( &
    (-1.0d0+2.0d0*nu)*sinb*( &
        (1.0d0/R*y3-cosb)/(R-z3) - &
        (1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb)/(Rbar+z3bar) &
    ) - &
    y1*(-1.0d0/R3*y3+1.0d0/2.0d0/R3bar*(2.0d0*y3+4.0d0*a)) - &
    sinb*(R*cosb-y3)/R/(R-z3) + &
    z1*(1.0d0/R*cosb*y3-1.0d0)/R/(R-z3) - &
    z1*(R*cosb-y3)/R3/(R-z3)*y3 - &
    z1*(R*cosb-y3)/R/(R-z3)**2*(1.0d0/R*y3-cosb) - &
    sinb*(Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar) - &
    z1bar*(1.0d0/2.0d0/Rbar*cosb*(2.0d0*y3+4.0d0*a)+1.0d0)/Rbar/(Rbar+z3bar) + &
    1.0d0/2.0d0*z1bar*(Rbar*cosb+y3bar)/R3bar/(Rbar+z3bar)*(2.0d0*y3+4.0d0*a) + &
    z1bar*(Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) &
)/pi/(1.0d0-nu)
return
end function

function e33_termB2b()
implicit none
double precision :: e33_termB2b
e33_termB2b = 1.0d0/4.0d0*( &
    (-2.0d0+2.0d0*nu)*(1.0d0-2.0d0*nu)*cotb*( &
        (1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0)/(Rbar+y3bar) - &
        cosb*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb)/(Rbar+z3bar) &
    ) + &
    (2.0d0-2.0d0*nu)*y1/(Rbar+y3bar)**2*(2.0d0*nu+a/Rbar)*( &
        1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a) + &
        1.0d0 &
    ) + &
    1.0d0/2.0d0*(2.0d0-2.0d0*nu)*y1/(Rbar+y3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) + &
    (2.0d0-2.0d0*nu)*sinb/(Rbar+z3bar)*(cosb+a/Rbar) - &
    (2.0d0-2.0d0*nu)*z1bar/(Rbar+z3bar)**2*(cosb+a/Rbar)*( &
        1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a) + &
        cosb &
    ) - &
    1.0d0/2.0d0*(2.0d0-2.0d0*nu)*z1bar/(Rbar+z3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) + &
    1.0d0/Rbar*((1.0d0-2.0d0*nu)*cotb-2.0d0*nu*y1/(Rbar+y3bar)-a*y1/R2bar) - &
    1.0d0/2.0d0*(y3+a)/R3bar*( &
        (1.0d0-2.0d0*nu)*cotb - &
        2.0d0*nu*y1/(Rbar+y3bar) - &
        a*y1/R2bar &
    )*(2.0d0*y3+4.0d0*a) + &
    (y3+a)/Rbar*( &
        2.0d0*nu*y1/(Rbar+y3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) + &
        a*y1/R4bar*(2.0d0*y3+4.0d0*a) &
    ) - &
    1.0d0/(Rbar+z3bar)*( &
        cosb*sinb+(Rbar*cosb+y3bar)*cotb/Rbar*( &
            (2.0d0-2.0d0*nu)*cosb - &
            (Rbar*cosb+y3bar)/(Rbar+z3bar) &
        ) + &
        a/Rbar*( &
            sinb - &
            y3bar*z1bar/R2bar - &
            z1bar*(Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar) &
        ) &
    ) + &
    (y3+a)/(Rbar+z3bar)**2*( &
        cosb*sinb + &
        (Rbar*cosb+y3bar)*cotb/Rbar*( &
            (2.0d0-2.0d0*nu)*cosb - &
            (Rbar*cosb+y3bar)/(Rbar+z3bar) &
        ) + &
        a/Rbar*( &
            sinb - &
            y3bar*z1bar/R2bar - &
            z1bar*(Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar) &
        ) &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) - &
    (y3+a)/(Rbar+z3bar)*( &
        (1.0d0/2.0d0/Rbar*cosb*(2.0d0*y3+4.0d0*a)+1.0d0)*cotb/Rbar*( &
            (2.0d0-2.0d0*nu)*cosb - &
            (Rbar*cosb+y3bar)/(Rbar+z3bar) &
        ) - &
        1.0d0/2.0d0*(Rbar*cosb+y3bar)*cotb/R3bar*( &
            (2.0d0-2.0d0*nu)*cosb - &
            (Rbar*cosb+y3bar)/(Rbar+z3bar) &
        )*(2.0d0*y3+4.0d0*a) + &
        (Rbar*cosb+y3bar)*cotb/Rbar*( &
            -(1.0d0/2.0d0/Rbar*cosb*(2.0d0*y3+4.0d0*a)+1.0d0)/(Rbar+z3bar) + &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) &
        ) - &
        1.0d0/2.0d0*a/R3bar*( &
            sinb - &
            y3bar*z1bar/R2bar - &
            z1bar*(Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar) &
        )*(2.0d0*y3+4.0d0*a) + &
        a/Rbar*( &
            -z1bar/R2bar - &
            y3bar*sinb/R2bar + &
            y3bar*z1bar/R4bar*(2.0d0*y3+4.0d0*a) - &
            sinb*(Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar) - &
            z1bar*(1.0d0/2.0d0/Rbar*cosb*(2.0d0*y3+4.0d0*a)+1.0d0)/Rbar/(Rbar+z3bar) + &
            1.0d0/2.0d0*z1bar*(Rbar*cosb+y3bar)/R3bar/(Rbar+z3bar)*(2.0d0*y3+4.0d0*a) + &
            z1bar*(Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar)**2*( &
                1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a) + &
                cosb &
            ) &
        ) &
    ) &
)/pi/(1.0d0-nu)
return
end function

function e33_termB3a()
implicit none
double precision :: e33_termB3a
e33_termB3a = 1.0d0/8.d0*( &
    (2.0d0-2.0d0*nu)*( &
        y2/z1**2*sinb/(1.0d0+y2**2/z1**2) + &
        ( &
            y2/R*sinb/(y1*z1+y2**2*cosb)*y3 + &
            y2*R*sinb**2/(y1*z1+y2**2*cosb)**2*y1 &
        )/(1.0d0+y2**2*R2*sinb**2/(y1*z1+y2**2*cosb)**2) + &
        y2/z1bar**2*sinb/(1.0d0+y2**2/z1bar**2) - &
        ( &
            1.0d0/2.0d0*y2/Rbar*sinb/(y1*z1bar+y2**2*cosb)*(2.0d0*y3+4.0d0*a) - &
            y2*Rbar*sinb**2/(y1*z1bar+y2**2*cosb)**2*y1 &
        )/(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
    ) + &
    y2*sinb*( &
        (1.0d0/R*cosb*y3-1.0d0)/R/(R-z3) - &
        (R*cosb-y3)/R3/(R-z3)*y3 - &
        (R*cosb-y3)/R/(R-z3)**2*(1.0d0/R*y3-cosb) - &
        (1.0d0/2.0d0/Rbar*cosb*(2.0d0*y3+4.0d0*a)+1.0d0)/Rbar/(Rbar+z3bar) + &
        1.0d0/2.0d0*(Rbar*cosb+y3bar)/R3bar/(Rbar+z3bar)*(2.0d0*y3+4.0d0*a) + &
        (Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) &
    ) &
)/pi/(1.0d0-nu)
return
end function

function e33_termB3b()
implicit none
double precision :: e33_termB3b
e33_termB3b = 1.0d0/4.0d0*( &
    (2.0d0-2.0d0*nu)*( &
        -y2/z1bar**2*sinb/(1.0d0+y2**2/z1bar**2) + &
        ( &
            1.0d0/2.0d0*y2/Rbar*sinb/(y1*z1bar+y2**2*cosb)*(2.0d0*y3+4.0d0*a) - &
            y2*Rbar*sinb**2/(y1*z1bar+y2**2*cosb)**2*y1 &
        )/(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
    ) - &
    (2.0d0-2.0d0*nu)*y2*sinb/(Rbar+z3bar)**2*(cosb+a/Rbar)*( &
        1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a) + &
        cosb &
    ) - &
    1.0d0/2.0d0*(2.0d0-2.0d0*nu)*y2*sinb/(Rbar+z3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) + &
    y2*sinb/Rbar/(Rbar+z3bar)*( &
        1.0d0 + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) + &
        a*y3bar/R2bar &
    ) - &
    1.0d0/2.0d0*y2*(y3+a)*sinb/R3bar/(Rbar+z3bar)*( &
        1.0d0 + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) + &
        a*y3bar/R2bar &
    )*(2.0d0*y3+4.0d0*a) - &
    y2*(y3+a)*sinb/Rbar/(Rbar+z3bar)**2*( &
        1.0d0 + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) + &
        a*y3bar/R2bar &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) + &
    y2*(y3+a)*sinb/Rbar/(Rbar+z3bar)*( &
        (1.0d0/2.0d0/Rbar*cosb*(2.0d0*y3+4.0d0*a)+1.0d0)/(Rbar+z3bar)*(cosb+a/Rbar) - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)**2*(cosb+a/Rbar)*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) - &
        1.0d0/2.0d0*(Rbar*cosb+y3bar)/(Rbar+z3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) + &
        a/R2bar-a*y3bar/R4bar*(2.0d0*y3+4.0d0*a) &
    ) &
)/pi/(1.0d0-nu)
return
end function

function e12_termB1a()
implicit none
double precision :: e12_termB1a, e12_termB1a_1, e12_termB1a_2
e12_termB1a_1 = 1.0d0/8.0d0*( &
    (2.0d0-2.0d0*nu)*( &
        -2.0d0/y1/(1.0d0+y2**2/y1**2) + &
        1.0d0/z1/(1.0d0+y2**2/z1**2) + &
        ( &
            R*sinb/(y1*z1+y2**2*cosb) + &
            y2**2/R*sinb/(y1*z1+y2**2*cosb) - &
            2.0d0*y2**2*R*sinb/(y1*z1+y2**2*cosb)**2*cosb &
        )/(1.0d0+y2**2*R2*sinb**2/(y1*z1+y2**2*cosb)**2) + &
        1.0d0/z1bar/(1.0d0+y2**2/z1bar**2)+( &
            Rbar*sinb/(y1*z1bar+y2**2*cosb) + &
            y2**2/Rbar*sinb/(y1*z1bar+y2**2*cosb) - &
            2.0d0*y2**2*Rbar*sinb/(y1*z1bar+y2**2*cosb)**2*cosb &
        )/(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
    ) - &
    y1*(1.0d0/R/(R-y3)+1.0d0/Rbar/(Rbar+y3bar)) - &
    y1*y2*( &
        -1.0d0/R3/(R-y3)*y2 - &
        1.0d0/R2/(R-y3)**2*y2 - &
        1.0d0/R3bar/(Rbar+y3bar)*y2 - &
        1.0d0/R2bar/(Rbar+y3bar)**2*y2 &
    ) - &
    cosb*((R*sinb-y1)/R/(R-z3)+(Rbar*sinb-y1)/Rbar/(Rbar+z3bar)) - &
    y2*cosb*( &
        1.0d0/R2*sinb*y2/(R-z3) - &
        (R*sinb-y1)/R3/(R-z3)*y2 - &
        (R*sinb-y1)/R2/(R-z3)**2*y2 + &
        1.0d0/R2bar*sinb*y2/(Rbar+z3bar) - &
        (Rbar*sinb-y1)/R3bar/(Rbar+z3bar)*y2 - &
        (Rbar*sinb-y1)/R2bar/(Rbar+z3bar)**2*y2 &
    ) &
)/pi/(1.0d0-nu)

e12_termB1a_2 = 1.0d0/4.0d0*( &
    (-2.0d0+2.0d0*nu)*(1.0d0-2.0d0*nu)*( &
        -1.0d0/y1/(1.0d0+y2**2/y1**2) + &
        1.0d0/z1bar/(1.0d0+y2**2/z1bar**2) + &
        ( &
            Rbar*sinb/(y1*z1bar+y2**2*cosb) + &
                y2**2/Rbar*sinb/(y1*z1bar+y2**2*cosb) - &
                2.0d0*y2**2*Rbar*sinb/(y1*z1bar+y2**2*cosb)**2*cosb &
         )/(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
    )*cotb**2 + &
    (1.0d0-2.0d0*nu)/(Rbar+y3bar)*((1.0d0-2.0d0*nu-a/Rbar)*cotb-y1/(Rbar+y3bar)*(nu+a/Rbar)) - &
    (1.0d0-2.0d0*nu)*y2**2/(Rbar+y3bar)**2*( &
        (1.0d0-2.0d0*nu-a/Rbar)*cotb - &
        y1/(Rbar+y3bar)*(nu+a/Rbar) &
    )/Rbar + &
    (1.0d0-2.0d0*nu)*y2/(Rbar+y3bar)*( &
        a/R3bar*y2*cotb + &
        y1/(Rbar+y3bar)**2*(nu+a/Rbar)/Rbar*y2 + &
        y2/(Rbar+y3bar)*a/R3bar*y1 &
    ) + &
    (1.0d0-2.0d0*nu)*cosb*cotb/(Rbar+z3bar)*(cosb+a/Rbar) - &
    (1.0d0-2.0d0*nu)*y2**2*cosb*cotb/(Rbar+z3bar)**2*(cosb+a/Rbar)/Rbar - &
    (1.0d0-2.0d0*nu)*y2**2*cosb*cotb/(Rbar+z3bar)*a/R3bar + &
    a*(y3+a)*cotb/R3bar - &
    3.0d0*a*y2**2*(y3+a)*cotb/R5bar + &
    (y3+a)/Rbar/(Rbar+y3bar)*( &
        (-1.0d0+2.0d0*nu)*cotb + &
        y1/(Rbar+y3bar)*(2.0d0*nu+a/Rbar) + &
        a*y1/R2bar &
    ) - &
    y2**2*(y3+a)/R3bar/(Rbar+y3bar)*( &
        (-1.0d0+2.0d0*nu)*cotb + &
        y1/(Rbar+y3bar)*(2.0d0*nu+a/Rbar) + &
        a*y1/R2bar &
    ) - &
    y2**2*(y3+a)/R2bar/(Rbar+y3bar)**2*( &
        (-1.0d0+2.0d0*nu)*cotb + &
        y1/(Rbar+y3bar)*(2.0d0*nu+a/Rbar) + &
        a*y1/R2bar &
    ) + &
    y2*(y3+a)/Rbar/(Rbar+y3bar)*( &
        -y1/(Rbar+y3bar)**2*(2.0d0*nu+a/Rbar)/Rbar*y2 - &
        y2/(Rbar+y3bar)*a/R3bar*y1 - &
        2.0d0*a*y1/R4bar*y2 &
    ) + &
    (y3+a)/Rbar/(Rbar+z3bar)*( &
        cosb/(Rbar+z3bar)*( &
            (Rbar*cosb+y3bar)*((1.0d0-2.0d0*nu)*cosb-a/Rbar)*cotb + &
            (2.0d0-2.0d0*nu)*(Rbar*sinb-y1)*cosb &
        ) - &
        a*y3bar*cosb*cotb/R2bar &
    ) - &
    y2**2*(y3+a)/R3bar/(Rbar+z3bar)*( &
        cosb/(Rbar+z3bar)*( &
            (Rbar*cosb+y3bar)*((1.0d0-2.0d0*nu)*cosb-a/Rbar)*cotb + &
            (2.0d0-2.0d0*nu)*(Rbar*sinb-y1)*cosb &
        ) - &
        a*y3bar*cosb*cotb/R2bar &
    ) - &
    y2**2*(y3+a)/R2bar/(Rbar+z3bar)**2*( &
        cosb/(Rbar+z3bar)*( &
            (Rbar*cosb+y3bar)*((1.0d0-2.0d0*nu)*cosb-a/Rbar)*cotb + &
            (2.0d0-2.0d0*nu)*(Rbar*sinb-y1)*cosb &
        ) - &
        a*y3bar*cosb*cotb/R2bar &
    ) + &
    y2*(y3+a)/Rbar/(Rbar+z3bar)*( &
        -cosb/(Rbar+z3bar)**2*( &
            (Rbar*cosb+y3bar)*((1.0d0-2.0d0*nu)*cosb-a/Rbar)*cotb + &
            (2.0d0-2.0d0*nu)*(Rbar*sinb-y1)*cosb &
        )/Rbar*y2 + &
        cosb/(Rbar+z3bar)*( &
            1.0d0/Rbar*cosb*y2*((1.0d0-2.0d0*nu)*cosb-a/Rbar)*cotb + &
            (Rbar*cosb+y3bar)*a/R3bar*y2*cotb + &
            (2.0d0-2.0d0*nu)/Rbar*sinb*y2*cosb &
        ) + &
        2.0d0*a*y3bar*cosb*cotb/R4bar*y2 &
    ) &
)/pi/(1.0d0-nu)

e12_termB1a = 1.0d0/2.0d0*(e12_termB1a_1+e12_termB1a_2)
return
end function

function e12_termB1b()
implicit none
double precision :: e12_termB1b, e12_termB1b_1, e12_termB1b_2
e12_termB1b_1 = 1.0d0/8.0d0*( &
    (1.0d0-2.0d0*nu)*( &
        1.0d0/R*y1/(R-y3) + &
        1.0d0/Rbar*y1/(Rbar+y3bar) - &
        cosb*((1.0d0/R*y1-sinb)/(R-z3)+(1.0d0/Rbar*y1-sinb)/(Rbar+z3bar)) &
    ) - &
    y2**2*( &
        -1.0d0/R3/(R-y3)*y1 - &
        1.0d0/R2/(R-y3)**2*y1 - &
        1.0d0/R3bar/(Rbar+y3bar)*y1 - &
        1.0d0/R2bar/(Rbar+y3bar)**2*y1 - &
        cosb*( &
            -1.0d0/R3/(R-z3)*y1 - &
            1.0d0/R/(R-z3)**2*(1.0d0/R*y1-sinb) - &
            1.0d0/R3bar/(Rbar+z3bar)*y1 - &
            1.0d0/Rbar/(Rbar+z3bar)**2*(1.0d0/Rbar*y1-sinb) &
        ) &
    ) &
)/pi/(1.0d0-nu)

e12_termB1b_2 = 1.0d0/4.0d0*( &
    (1.0d0-2.0d0*nu)*( &
        ((2.0d0-2.0d0*nu)*cotb**2-nu)/Rbar*y1/(Rbar+y3bar) - &
        ((2.0d0-2.0d0*nu)*cotb**2+1.0d0-2.0d0*nu)*cosb*(1.0d0/Rbar*y1-sinb)/(Rbar+z3bar) &
    ) + &
    (1.0d0-2.0d0*nu)/(Rbar+y3bar)**2*( &
        y1*cotb*(1.0d0-2.0d0*nu-a/Rbar) + &
        nu*y3bar - &
        a + &
        y2**2/(Rbar+y3bar)*(nu+a/Rbar) &
    )/Rbar*y1 - &
    (1.0d0-2.0d0*nu)/(Rbar+y3bar)*( &
        (1.0d0-2.0d0*nu-a/Rbar)*cotb + &
        a*y1**2*cotb/R3bar - &
        y2**2/(Rbar+y3bar)**2*(nu+a/Rbar)/Rbar*y1 - &
        y2**2/(Rbar+y3bar)*a/R3bar*y1 &
    ) - &
    (1.0d0-2.0d0*nu)*cosb*cotb/(Rbar+z3bar)*(cosb+a/Rbar) + &
    (1.0d0-2.0d0*nu)*z1bar*cotb/(Rbar+z3bar)**2*(cosb+a/Rbar)*(1.0d0/Rbar*y1-sinb) + &
    (1.0d0-2.0d0*nu)*z1bar*cotb/(Rbar+z3bar)*a/R3bar*y1 - &
    a*(y3+a)*cotb/R3bar + &
    3.0d0*a*y1**2*(y3+a)*cotb/R5bar - &
    (y3+a)/(Rbar+y3bar)**2*( &
        -2.0d0*nu + &
        1.0d0/Rbar*((1.0d0-2.0d0*nu)*y1*cotb-a) + &
        y2**2/Rbar/(Rbar+y3bar)*(2*nu+a/Rbar) + &
        a*y2**2/R3bar &
    )/Rbar*y1 + &
    (y3+a)/(Rbar+y3bar)*( &
        -1.0d0/R3bar*((1.0d0-2.0d0*nu)*y1*cotb-a)*y1 + &
        1.0d0/Rbar*(1.0d0-2.0d0*nu)*cotb - &
        y2**2/R3bar/(Rbar+y3bar)*(2.0d0*nu+a/Rbar)*y1 - &
        y2**2/R2bar/(Rbar+y3bar)**2*(2.0d0*nu+a/Rbar)*y1 - &
        y2**2/R4bar/(Rbar+y3bar)*a*y1 - &
        3.0d0*a*y2**2/R5bar*y1 &
    ) - &
    (y3+a)/(Rbar+z3bar)**2*( &
        cosb**2 - &
        1.0d0/Rbar*((1.0d0-2.0d0*nu)*z1bar*cotb+a*cosb) + &
        a*y3bar*z1bar*cotb/R3bar - &
        1.0d0/Rbar/(Rbar+z3bar)*(y2**2*cosb**2 - &
        a*z1bar*cotb/Rbar*(Rbar*cosb+y3bar)) &
    )*(1.0d0/Rbar*y1-sinb) + &
    (y3+a)/(Rbar+z3bar)*( &
        1.0d0/R3bar*((1.0d0-2.0d0*nu)*z1bar*cotb+a*cosb)*y1 - &
        1.0d0/Rbar*(1.0d0-2.0d0*nu)*cosb*cotb + &
        a*y3bar*cosb*cotb/R3bar - &
        3.0d0*a*y3bar*z1bar*cotb/R5bar*y1 + &
        1.0d0/R3bar/(Rbar+z3bar)*(y2**2*cosb**2-a*z1bar*cotb/Rbar*(Rbar*cosb+y3bar))*y1 + &
        1.0d0/Rbar/(Rbar+z3bar)**2*( &
            y2**2*cosb**2 - &
            a*z1bar*cotb/Rbar*(Rbar*cosb+y3bar) &
        )*(1.0d0/Rbar*y1-sinb) - &
        1.0d0/Rbar/(Rbar+z3bar)*( &
            -a*cosb*cotb/Rbar*(Rbar*cosb+y3bar) + &
            a*z1bar*cotb/R3bar*(Rbar*cosb+y3bar)*y1 - &
            a*z1bar*cotb/R2bar*cosb*y1 &
        ) &
    ) &
)/pi/(1.0d0-nu)

e12_termB1b = 1.0d0/2.0d0*(e12_termB1b_1+e12_termB1b_2)
return
end function

function e12_termB2a()
implicit none
double precision :: e12_termB2a, e12_termB2a_1, e12_termB2a_2
e12_termB2a_1 = 1.0d0/8.0d0*( &
    (-1+2*nu)*( &
        1.0d0/R*y2/(R-y3) + &
        1.0d0/Rbar*y2/(Rbar+y3bar) - &
        cosb*(1.0d0/R*y2/(R-z3)+1/Rbar*y2/(Rbar+z3bar)) &
    ) + &
    y1**2*( &
        -1/R3/(R-y3)*y2 - &
        1.0d0/R2/(R-y3)**2*y2 - &
        1.0d0/R3bar/(Rbar+y3bar)*y2 - &
        1.0d0/R2bar/(Rbar+y3bar)**2*y2 &
    ) + &
    z1/R2*sinb*y2/(R-z3) - &
    z1*(R*sinb-y1)/R3/(R-z3)*y2 - &
    z1*(R*sinb-y1)/R2/(R-z3)**2*y2 + &
    z1bar/R2bar*sinb*y2/(Rbar+z3bar) - &
    z1bar*(Rbar*sinb-y1)/R3bar/(Rbar+z3bar)*y2 - &
    z1bar*(Rbar*sinb-y1)/R2bar/(Rbar+z3bar)**2*y2 &
)/pi/(1.0d0-nu)

e12_termB2a_2 = 1.0d0/4.0d0*( &
    (1.0d0-2.0d0*nu)*( &
        ((2.0d0-2.0d0*nu)*cotb**2+nu)/Rbar*y2/(Rbar+y3bar) - &
        ((2.0d0-2.0d0*nu)*cotb**2+1.0d0)*cosb/Rbar*y2/(Rbar+z3bar) &
    ) - &
    (1.0d0-2.0d0*nu)/(Rbar+y3bar)**2*( &
        (-1+2*nu)*y1*cotb + &
        nu*y3bar - &
        a + &
        a*y1*cotb/Rbar + &
        y1**2/(Rbar+y3bar)*(nu+a/Rbar) &
    )/Rbar*y2 + &
    (1.0d0-2.0d0*nu)/(Rbar+y3bar)*( &
        -a*y1*cotb/R3bar*y2 - &
        y1**2/(Rbar+y3bar)**2*(nu+a/Rbar)/Rbar*y2 - &
        y1**2/(Rbar+y3bar)*a/R3bar*y2 &
    ) + &
    (1.0d0-2.0d0*nu)*cotb/(Rbar+z3bar)**2*( &
        z1bar*cosb - &
        a*(Rbar*sinb-y1)/Rbar/cosb &
    )/Rbar*y2 - &
    (1.0d0-2.0d0*nu)*cotb/(Rbar+z3bar)*( &
        -a/R2bar*sinb*y2/cosb + &
        a*(Rbar*sinb-y1)/R3bar/cosb*y2 &
    ) + &
    3.0d0*a*y2*(y3+a)*cotb/R5bar*y1 - &
    (y3+a)/(Rbar+y3bar)**2*( &
        2.0d0*nu + &
        1.0d0/Rbar*((1.0d0-2.0d0*nu)*y1*cotb+a) - &
        y1**2/Rbar/(Rbar+y3bar)*(2*nu+a/Rbar) - &
        a*y1**2/R3bar &
    )/Rbar*y2 + &
    (y3+a)/(Rbar+y3bar)*( &
        -1/R3bar*((1.0d0-2.0d0*nu)*y1*cotb+a)*y2 + &
        y1**2/R3bar/(Rbar+y3bar)*(2*nu+a/Rbar)*y2 + &
        y1**2/R2bar/(Rbar+y3bar)**2*(2*nu+a/Rbar)*y2 + &
        y1**2/R4bar/(Rbar+y3bar)*a*y2 + &
        3.0d0*a*y1**2/R5bar*y2 &
    ) - &
    (y3+a)*cotb/(Rbar+z3bar)**2*( &
        -cosb*sinb + &
        a*y1*y3bar/R3bar/cosb + &
        (Rbar*sinb-y1)/Rbar*( &
            (2.0d0-2.0d0*nu)*cosb - &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) &
        ) &
    )/Rbar*y2 + &
    (y3+a)*cotb/(Rbar+z3bar)*( &
        -3*a*y1*y3bar/R5bar/cosb*y2 + &
        1.0d0/R2bar*sinb*y2*( &
            (2.0d0-2.0d0*nu)*cosb - &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) &
        ) - &
        (Rbar*sinb-y1)/R3bar*( &
            (2.0d0-2.0d0*nu)*cosb - &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) &
        )*y2 + &
        (Rbar*sinb-y1)/Rbar*( &
            -1/Rbar*cosb*y2/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) + &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)**2*(1.0d0+a/Rbar/cosb)/Rbar*y2 + &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)*a/R3bar/cosb*y2 &
        ) &
    ) &
)/pi/(1.0d0-nu)

e12_termB2a = 1.0d0/2.0d0*(e12_termB2a_1+e12_termB2a_2)
return
end function

function e12_termB2b()
implicit none
double precision :: e12_termB2b, e12_termB2b_1, e12_termB2b_2
e12_termB2b_1 = 1.0d0/8.0d0*( &
    (2.0d0-2.0d0*nu)*( &
        2.0d0*y2/y1**2/(1.0d0+y2**2/y1**2) - &
        y2/z1**2*cosb/(1.0d0+y2**2/z1**2) + &
        ( &
            y2/R*sinb/(y1*z1+y2**2*cosb)*y1 - &
            y2*R*sinb/(y1*z1+y2**2*cosb)**2*(2.0d0*y1*cosb-y3*sinb) &
        )/(1.0d0+y2**2*R2*sinb**2/(y1*z1+y2**2*cosb)**2) - &
        y2/z1bar**2*cosb/(1.0d0+y2**2/z1bar**2) + &
        ( &
            y2/Rbar*sinb/(y1*z1bar+y2**2*cosb)*y1 - &
            y2*Rbar*sinb/(y1*z1bar+y2**2*cosb)**2*(2.0d0*y1*cosb+y3bar*sinb) &
        )/(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
    ) + &
    y2*(1.0d0/R/(R-y3)+1.0d0/Rbar/(Rbar+y3bar)) + &
    y1*y2*( &
        -1.0d0/R3/(R-y3)*y1 - &
        1.0d0/R2/(R-y3)**2*y1 - &
        1.0d0/R3bar/(Rbar+y3bar)*y1 - &
        1.0d0/R2bar/(Rbar+y3bar)**2*y1 &
    ) - &
    y2*( &
        cosb/R/(R-z3) - &
        z1/R3/(R-z3)*y1 - &
        z1/R/(R-z3)**2*(1.0d0/R*y1-sinb) + &
        cosb/Rbar/(Rbar+z3bar) - &
        z1bar/R3bar/(Rbar+z3bar)*y1 - &
        z1bar/Rbar/(Rbar+z3bar)**2*(1.0d0/Rbar*y1-sinb) &
    ) &
)/pi/(1.0d0-nu)

e12_termB2b_2 = 1.0d0/4.0d0*( &
    (2.0d0-2.0d0*nu)*(1.0d0-2.0d0*nu)*( &
        y2/y1**2/(1.0d0+y2**2/y1**2) - &
        y2/z1bar**2*cosb/(1.0d0+y2**2/z1bar**2) + &
        ( &
            y2/Rbar*sinb/(y1*z1bar+y2**2*cosb)*y1 - &
            y2*Rbar*sinb/(y1*z1bar+y2**2*cosb)**2*(2.0d0*y1*cosb+y3bar*sinb) &
        )/(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
    )*cotb**2 - &
    (1.0d0-2.0d0*nu)*y2/(Rbar+y3bar)**2*( &
        (-1.0d0+2.0d0*nu+a/Rbar)*cotb + &
        y1/(Rbar+y3bar)*(nu+a/Rbar) &
    )/Rbar*y1 + &
    (1.0d0-2.0d0*nu)*y2/(Rbar+y3bar)*( &
        -a/R3bar*y1*cotb + &
        1.0d0/(Rbar+y3bar)*(nu+a/Rbar) - &
        y1**2/(Rbar+y3bar)**2*(nu+a/Rbar)/Rbar - &
        y1**2/(Rbar+y3bar)*a/R3bar &
    ) + &
    (1.0d0-2.0d0*nu)*y2*cotb/(Rbar+z3bar)**2*(1.0d0+a/Rbar/cosb)*(1.0d0/Rbar*y1-sinb) + &
    (1.0d0-2.0d0*nu)*y2*cotb/(Rbar+z3bar)*a/R3bar/cosb*y1 + &
    3.0d0*a*y2*(y3+a)*cotb/R5bar*y1 - &
    y2*(y3+a)/R3bar/(Rbar+y3bar)*( &
        (1.0d0-2.0d0*nu)*cotb - &
        2.0d0*nu*y1/(Rbar+y3bar) - &
        a*y1/Rbar*(1.0d0/Rbar+1.0d0/(Rbar+y3bar)) &
    )*y1 - &
    y2*(y3+a)/R2bar/(Rbar+y3bar)**2*( &
        (1.0d0-2.0d0*nu)*cotb - &
        2.0d0*nu*y1/(Rbar+y3bar) - &
        a*y1/Rbar*(1.0d0/Rbar+1.0d0/(Rbar+y3bar)) &
    )*y1 + &
    y2*(y3+a)/Rbar/(Rbar+y3bar)*( &
        -2.0d0*nu/(Rbar+y3bar) + &
        2.0d0*nu*y1**2/(Rbar+y3bar)**2/Rbar - &
        a/Rbar*(1.0d0/Rbar+1.0d0/(Rbar+y3bar)) + &
        a*y1**2/R3bar*(1.0d0/Rbar+1.0d0/(Rbar+y3bar)) - &
        a*y1/Rbar*(-1.0d0/R3bar*y1-1.0d0/(Rbar+y3bar)**2/Rbar*y1) &
    ) - &
    y2*(y3+a)*cotb/R3bar/(Rbar+z3bar)*( &
        (-2.0d0+2.0d0*nu)*cosb + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) + &
        a*y3bar/R2bar/cosb &
    )*y1 - &
    y2*(y3+a)*cotb/Rbar/(Rbar+z3bar)**2*( &
        (-2.0d0+2.0d0*nu)*cosb + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) + &
        a*y3bar/R2bar/cosb &
    )*(1.0d0/Rbar*y1-sinb) + &
    y2*(y3+a)*cotb/Rbar/(Rbar+z3bar)*( &
        1.0d0/Rbar*cosb*y1/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)**2*(1.0d0+a/Rbar/cosb)*(1.0d0/Rbar*y1-sinb) - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*a/R3bar/cosb*y1 - &
        2.0d0*a*y3bar/R4bar/cosb*y1 &
    ) &
)/pi/(1.0d0-nu)

e12_termB2b = 1.0d0/2.0d0*(e12_termB2b_1+e12_termB2b_2)
return
end function

function e12_termB3a()
implicit none
double precision :: e12_termB3a, e12_termB3a_1, e12_termB3a_2
e12_termB3a_1 = 1.0d0/8.0d0*( &
    sinb*((R*sinb-y1)/R/(R-z3)+(Rbar*sinb-y1)/Rbar/(Rbar+z3bar)) + &
    y2*sinb*( &
        1.0d0/R2*sinb*y2/(R-z3) - &
        (R*sinb-y1)/R3/(R-z3)*y2 - &
        (R*sinb-y1)/R2/(R-z3)**2*y2 + &
        1.0d0/R2bar*sinb*y2/(Rbar+z3bar) - &
        (Rbar*sinb-y1)/R3bar/(Rbar+z3bar)*y2 - &
        (Rbar*sinb-y1)/R2bar/(Rbar+z3bar)**2*y2 &
    ) &
)/pi/(1.0d0-nu)

e12_termB3a_2 = 1.0d0/4.0d0*( &
    (1.0d0-2.0d0*nu)*( &
        1.0d0/(Rbar+y3bar)*(1.0d0+a/Rbar) - &
        y2**2/(Rbar+y3bar)**2*(1.0d0+a/Rbar)/Rbar - &
        y2**2/(Rbar+y3bar)*a/R3bar - &
        cosb/(Rbar+z3bar)*(cosb+a/Rbar) + &
        y2**2*cosb/(Rbar+z3bar)**2*(cosb+a/Rbar)/Rbar + &
        y2**2*cosb/(Rbar+z3bar)*a/R3bar &
    ) - &
    (y3+a)/Rbar*(a/R2bar+1.0d0/(Rbar+y3bar)) + &
    y2**2*(y3+a)/R3bar*(a/R2bar+1.0d0/(Rbar+y3bar)) - &
    y2*(y3+a)/Rbar*(-2.0d0*a/R4bar*y2-1.0d0/(Rbar+y3bar)**2/Rbar*y2) + &
    (y3+a)*cosb/Rbar/(Rbar+z3bar)*((Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar)+a*y3bar/R2bar) - &
    y2**2*(y3+a)*cosb/R3bar/(Rbar+z3bar)*( &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) + &
        a*y3bar/R2bar &
    ) - &
    y2**2*(y3+a)*cosb/R2bar/(Rbar+z3bar)**2*( &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar)+a*y3bar/R2bar &
    ) + &
    y2*(y3+a)*cosb/Rbar/(Rbar+z3bar)*( &
        1.0d0/Rbar*cosb*y2/(Rbar+z3bar)*(cosb+a/Rbar) - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)**2*(cosb+a/Rbar)/Rbar*y2 - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*a/R3bar*y2 - &
        2.0d0*a*y3bar/R4bar*y2 &
    ) &
)/pi/(1.0d0-nu)

e12_termB3a = 1.0d0/2.0d0*(e12_termB3a_1+e12_termB3a_2)
return
end function

function e12_termB3b()
implicit none
double precision :: e12_termB3b, e12_termB3b_1, e12_termB3b_2
e12_termB3b_1 = 1.0d0/8.0d0*( &
    (1.0d0-2.0d0*nu)*sinb*((1.0d0/R*y1-sinb)/(R-z3)+(1.0d0/Rbar*y1-sinb)/(Rbar+z3bar)) - &
    y2**2*sinb*( &
        -1.0d0/R3/(R-z3)*y1 - &
        1.0d0/R/(R-z3)**2*(1.0d0/R*y1-sinb) - &
        1.0d0/R3bar/(Rbar+z3bar)*y1 - &
        1.0d0/Rbar/(Rbar+z3bar)**2*(1.0d0/Rbar*y1-sinb) &
    ) &
)/pi/(1.0d0-nu)

e12_termB3b_2 = 1.0d0/4.0d0*( &
    (1.0d0-2.0d0*nu)*( &
        -sinb*(1.0d0/Rbar*y1-sinb)/(Rbar+z3bar) - &
        1.0d0/(Rbar+y3bar)*(1.0d0+a/Rbar) + &
        y1**2/(Rbar+y3bar)**2*(1.0d0+a/Rbar)/Rbar + &
        y1**2/(Rbar+y3bar)*a/R3bar+cosb/(Rbar+z3bar)*(cosb+a/Rbar) - &
        z1bar/(Rbar+z3bar)**2*(cosb+a/Rbar)*(1.0d0/Rbar*y1-sinb) - &
        z1bar/(Rbar+z3bar)*a/R3bar*y1 &
    ) + &
    (y3+a)/Rbar*(a/R2bar+1.0d0/(Rbar+y3bar)) - &
    y1**2*(y3+a)/R3bar*(a/R2bar+1.0d0/(Rbar+y3bar)) + &
    y1*(y3+a)/Rbar*(-2.0d0*a/R4bar*y1-1.0d0/(Rbar+y3bar)**2/Rbar*y1) + &
    (y3+a)/(Rbar+z3bar)**2*( &
        sinb*(cosb-a/Rbar) + &
        z1bar/Rbar*(1.0d0+a*y3bar/R2bar) - &
        1.0d0/Rbar/(Rbar+z3bar)*(y2**2*cosb*sinb-a*z1bar/Rbar*(Rbar*cosb+y3bar)) &
    )*(1.0d0/Rbar*y1-sinb) - &
    (y3+a)/(Rbar+z3bar)*( &
        sinb*a/R3bar*y1 + &
        cosb/Rbar*(1.0d0+a*y3bar/R2bar) - &
        z1bar/R3bar*(1.0d0+a*y3bar/R2bar)*y1 - &
        2.0d0*z1bar/R5bar*a*y3bar*y1 + &
        1.0d0/R3bar/(Rbar+z3bar)*(y2**2*cosb*sinb-a*z1bar/Rbar*(Rbar*cosb+y3bar))*y1 + &
        1.0d0/Rbar/(Rbar+z3bar)**2*( &
            y2**2*cosb*sinb - &
            a*z1bar/Rbar*(Rbar*cosb+y3bar) &
        )*(1.0d0/Rbar*y1-sinb) - &
        1.0d0/Rbar/(Rbar+z3bar)*( &
            -a*cosb/Rbar*(Rbar*cosb+y3bar) + &
            a*z1bar/R3bar*(Rbar*cosb+y3bar)*y1 - &
            a*z1bar/R2bar*cosb*y1 &
        ) &
    ) &
)/pi/(1.0d0-nu)

e12_termB3b = 1.0d0/2.0d0*(e12_termB3b_1+e12_termB3b_2)
return
end function

function e13_termB1a()
implicit none
double precision :: e13_termB1a, e13_termB1a_1, e13_termB1a_2
e13_termB1a_1 = 1.0d0/8.0d0*( &
    (2.0d0-2.0d0*nu)*( &
        y2/z1**2*sinb/(1.0d0+y2**2/z1**2) + &
        ( &
            y2/R*sinb/(y1*z1+y2**2*cosb)*y3 + &
            y2*R*sinb**2/(y1*z1+y2**2*cosb)**2*y1 &
        )/(1.0d0+y2**2*R2*sinb**2/(y1*z1+y2**2*cosb)**2) - &
        y2/z1bar**2*sinb/(1.0d0+y2**2/z1bar**2) + &
        ( &
            1.0d0/2.0d0*y2/Rbar*sinb/(y1*z1bar+y2**2*cosb)*(2.0d0*y3+4.0d0*a) - &
            y2*Rbar*sinb**2/(y1*z1bar+y2**2*cosb)**2*y1 &
        )/(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
    ) - &
    y1*y2*( &
        -1.0d0/R3/(R-y3)*y3 - &
        1.0d0/R/(R-y3)**2*(1.0d0/R*y3-1.0d0) - &
        1.0d0/2.0d0/R3bar/(Rbar+y3bar)*(2.0d0*y3+4.0d0*a) - &
        1.0d0/Rbar/(Rbar+y3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) &
    ) - &
    y2*cosb*( &
        1.0d0/R2*sinb*y3/(R-z3) - &
        (R*sinb-y1)/R3/(R-z3)*y3 - &
        (R*sinb-y1)/R/(R-z3)**2*(1.0d0/R*y3-cosb) + &
        1.0d0/2.0d0/R2bar*sinb*(2*y3+4*a)/(Rbar+z3bar) - &
        1.0d0/2.0d0*(Rbar*sinb-y1)/R3bar/(Rbar+z3bar)*(2*y3+4*a) - &
        (Rbar*sinb-y1)/Rbar/(Rbar+z3bar)**2*(1.0d0/2.0d0/Rbar*(2*y3+4*a)+cosb) &
    ) &
)/pi/(1.0d0-nu)

e13_termB1a_2 = 1.0d0/4.0d0*( &
    (-2.0d0+2.0d0*nu)*(1.0d0-2.0d0*nu)*( &
        -y2/z1bar**2*sinb/(1.0d0+y2**2/z1bar**2) + &
        ( &
            1.0d0/2.0d0*y2/Rbar*sinb/(y1*z1bar+y2**2*cosb)*(2.0d0*y3+4.0d0*a) - &
            y2*Rbar*sinb**2/(y1*z1bar+y2**2*cosb)**2*y1 &
        )/(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
    )*cotb**2 - &
    (1.0d0-2.0d0*nu)*y2/(Rbar+y3bar)**2*( &
        (1.0d0-2.0d0*nu-a/Rbar)*cotb - &
        y1/(Rbar+y3bar)*(nu+a/Rbar) &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) + &
    (1.0d0-2.0d0*nu)*y2/(Rbar+y3bar)*( &
        1.0d0/2.0d0*a/R3bar*(2.0d0*y3+4.0d0*a)*cotb + &
        y1/(Rbar+y3bar)**2*(nu+a/Rbar)*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) + &
        1.0d0/2.0d0*y1/(Rbar+y3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) &
    ) - &
    (1.0d0-2.0d0*nu)*y2*cosb*cotb/(Rbar+z3bar)**2*(cosb+a/Rbar)*( &
        1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a) + &
        cosb &
    ) - &
    1.0d0/2.0d0*(1.0d0-2.0d0*nu)*y2*cosb*cotb/(Rbar+z3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) + &
    a/R3bar*y2*cotb - &
    3.0d0/2.0d0*a*y2*(y3+a)*cotb/R5bar*(2.0d0*y3+4.0d0*a) + &
    y2/Rbar/(Rbar+y3bar)*((-1.0d0+2.0d0*nu)*cotb+y1/(Rbar+y3bar)*(2.0d0*nu+a/Rbar)+a*y1/R2bar) - &
    1.0d0/2.0d0*y2*(y3+a)/R3bar/(Rbar+y3bar)*( &
        (-1.0d0+2.0d0*nu)*cotb + &
        y1/(Rbar+y3bar)*(2.0d0*nu+a/Rbar) + &
        a*y1/R2bar &
    )*(2.0d0*y3+4.0d0*a) - &
    y2*(y3+a)/Rbar/(Rbar+y3bar)**2*( &
        (-1.0d0+2.0d0*nu)*cotb+y1/(Rbar+y3bar)*(2.0d0*nu+a/Rbar)+a*y1/R2bar &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) + &
    y2*(y3+a)/Rbar/(Rbar+y3bar)*( &
        -y1/(Rbar+y3bar)**2*(2.0d0*nu+a/Rbar)*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) - &
        1.0d0/2.0d0*y1/(Rbar+y3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) - &
        a*y1/R4bar*(2.0d0*y3+4.0d0*a) &
    ) + &
    y2/Rbar/(Rbar+z3bar)*( &
        cosb/(Rbar+z3bar)*( &
            (Rbar*cosb+y3bar)*((1.0d0-2.0d0*nu)*cosb-a/Rbar)*cotb + &
            (2.0d0-2.0d0*nu)*(Rbar*sinb-y1)*cosb &
        )-a*y3bar*cosb*cotb/R2bar &
    ) - &
    1.0d0/2.0d0*y2*(y3+a)/R3bar/(Rbar+z3bar)*( &
        cosb/(Rbar+z3bar)*( &
            (Rbar*cosb+y3bar)*((1.0d0-2.0d0*nu)*cosb-a/Rbar)*cotb + &
            (2.0d0-2.0d0*nu)*(Rbar*sinb-y1)*cosb &
        ) - &
        a*y3bar*cosb*cotb/R2bar &
    )*(2.0d0*y3+4.0d0*a) - &
    y2*(y3+a)/Rbar/(Rbar+z3bar)**2*( &
        cosb/(Rbar+z3bar)*( &
            (Rbar*cosb+y3bar)*((1.0d0-2.0d0*nu)*cosb-a/Rbar)*cotb + &
            (2.0d0-2.0d0*nu)*(Rbar*sinb-y1)*cosb &
        ) - &
        a*y3bar*cosb*cotb/R2bar &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) + &
    y2*(y3+a)/Rbar/(Rbar+z3bar)*( &
        -cosb/(Rbar+z3bar)**2*((Rbar*cosb+y3bar)*((1.0d0-2.0d0*nu)*cosb-a/Rbar)*cotb + &
        (2.0d0-2.0d0*nu)*(Rbar*sinb-y1)*cosb)*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) + &
        cosb/(Rbar+z3bar)*( &
            (1.0d0/2.0d0/Rbar*cosb*(2.0d0*y3+4.0d0*a)+1.0d0)*((1.0d0-2.0d0*nu)*cosb-a/Rbar)*cotb + &
            1.0d0/2.0d0*(Rbar*cosb+y3bar)*a/R3bar*(2.0d0*y3+4.0d0*a)*cotb + &
            1.0d0/2.0d0*(2.0d0-2.0d0*nu)/Rbar*sinb*(2.0d0*y3+4.0d0*a)*cosb &
        ) - &
        a*cosb*cotb/R2bar + &
        a*y3bar*cosb*cotb/R4bar*(2.0d0*y3+4.0d0*a) &
    ) &
)/pi/(1.0d0-nu)

e13_termB1a = 1.0d0/2.0d0*(e13_termB1a_1+e13_termB1a_2)
return
end function

function e13_termB1b()
implicit none
double precision :: e13_termB1b, e13_termB1b_1, e13_termB1b_2
e13_termB1b_1 = 1.0d0/8.0d0*y2*( &
    -1.0d0/R3*y1 + &
    1.0d0/R3bar*y1 - &
    cosb*( &
        1.0d0/R2*cosb*y1/(R-z3) - &
        (R*cosb-y3)/R3/(R-z3)*y1 - &
        (R*cosb-y3)/R/(R-z3)**2*(1.0d0/R*y1-sinb) - &
        1.0d0/R2bar*cosb*y1/(Rbar+z3bar) + &
        (Rbar*cosb+y3bar)/R3bar/(Rbar+z3bar)*y1 + &
        (Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar)**2*(1.0d0/Rbar*y1-sinb) &
    ) &
)/pi/(1.0d0-nu)

e13_termB1b_2 = 1.0d0/4.0d0*( &
    (2.0d0-2.0d0*nu)*( &
        (1.0d0-2.0d0*nu)*( &
            y2/y1**2/(1.0d0+y2**2/y1**2) - &
            y2/z1bar**2*cosb/(1.0d0+y2**2/z1bar**2) + &
            ( &
                y2/Rbar*sinb/(y1*z1bar+y2**2*cosb)*y1 - &
                y2*Rbar*sinb/(y1*z1bar+y2**2*cosb)**2*(2.0d0*y1*cosb+y3bar*sinb) &
            )/(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
        )*cotb - &
        y1/(Rbar+y3bar)**2*(2.0d0*nu+a/Rbar)/Rbar*y2 - &
        y2/(Rbar+y3bar)*a/R3bar*y1 + &
        y2*cosb/(Rbar+z3bar)**2*(cosb+a/Rbar)*(1.0d0/Rbar*y1-sinb) + &
        y2*cosb/(Rbar+z3bar)*a/R3bar*y1 &
    ) - &
    y2*(y3+a)/R3bar*(2.0d0*nu/(Rbar+y3bar)+a/R2bar)*y1 + &
    y2*(y3+a)/Rbar*(-2.0d0*nu/(Rbar+y3bar)**2/Rbar*y1-2.0d0*a/R4bar*y1) - &
    y2*(y3+a)*cosb/R3bar/(Rbar+z3bar)*( &
        1.0d0 - &
        2.0d0*nu - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) - &
        a*y3bar/R2bar &
    )*y1 - &
    y2*(y3+a)*cosb/Rbar/(Rbar+z3bar)**2*( &
        1.0d0 - &
        2.0d0*nu - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) - &
        a*y3bar/R2bar &
    )*(1.0d0/Rbar*y1-sinb) + &
    y2*(y3+a)*cosb/Rbar/(Rbar+z3bar)*( &
        -1.0d0/Rbar*cosb*y1/(Rbar+z3bar)*(cosb+a/Rbar) + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)**2*(cosb+a/Rbar)*(1.0d0/Rbar*y1-sinb) + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*a/R3bar*y1+2.0d0*a*y3bar/R4bar*y1) &
)/pi/(1.0d0-nu)

e13_termB1b = 1.0d0/2.0d0*(e13_termB1b_1+e13_termB1b_2)
return
end function

function e13_termB2a()
implicit none
double precision :: e13_termB2a, e13_termB2a_1, e13_termB2a_2
e13_termB2a_1 = 1.0d0/8.0d0*( &
    (-1.0d0+2.0d0*nu)*( &
        (1.0d0/R*y3-1.0d0)/(R-y3) + &
        (1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0)/(Rbar+y3bar) - &
        cosb*((1.0d0/R*y3-cosb)/(R-z3)+(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb)/(Rbar+z3bar)) &
    ) + &
    y1**2*( &
        -1.0d0/R3/(R-y3)*y3 - &
        1.0d0/R/(R-y3)**2*(1.0d0/R*y3-1.0d0) - &
        1.0d0/2.0d0/R3bar/(Rbar+y3bar)*(2.0d0*y3+4.0d0*a) - &
        1.0d0/Rbar/(Rbar+y3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) &
    ) - &
    sinb*(R*sinb-y1)/R/(R-z3) + &
    z1/R2*sinb*y3/(R-z3) - &
    z1*(R*sinb-y1)/R3/(R-z3)*y3 - &
    z1*(R*sinb-y1)/R/(R-z3)**2*(1.0d0/R*y3-cosb) + &
    sinb*(Rbar*sinb-y1)/Rbar/(Rbar+z3bar) + &
    1.0d0/2.0d0*z1bar/R2bar*sinb*(2.0d0*y3+4.0d0*a)/(Rbar+z3bar) - &
    1.0d0/2.0d0*z1bar*(Rbar*sinb-y1)/R3bar/(Rbar+z3bar)*(2.0d0*y3+4.0d0*a) - &
    z1bar*(Rbar*sinb-y1)/Rbar/(Rbar+z3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) &
)/pi/(1.0d0-nu)

e13_termB2a_2 = 1.0d0/4.0d0*( &
    (1.0d0-2.0d0*nu)*( &
        ((2.0d0-2.0d0*nu)*cotb**2+nu)*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0)/(Rbar+y3bar) - &
        ((2.0d0-2.0d0*nu)*cotb**2+1.0d0)*cosb*( &
            1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a) + &
            cosb &
        )/(Rbar+z3bar) &
    ) - &
    (1.0d0-2.0d0*nu)/(Rbar+y3bar)**2*( &
        (-1.0d0+2.0d0*nu)*y1*cotb + &
        nu*y3bar - &
        a + &
        a*y1*cotb/Rbar + &
        y1**2/(Rbar+y3bar)*(nu+a/Rbar) &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) + &
    (1.0d0-2.0d0*nu)/(Rbar+y3bar)*( &
        nu - &
        1.0d0/2.0d0*a*y1*cotb/R3bar*(2.0d0*y3+4.0d0*a) - &
        y1**2/(Rbar+y3bar)**2*(nu+a/Rbar)*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) - &
        1.0d0/2.0d0*y1**2/(Rbar+y3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) &
    ) + &
    (1.0d0-2.0d0*nu)*cotb/(Rbar+z3bar)**2*( &
        z1bar*cosb - &
        a*(Rbar*sinb-y1)/Rbar/cosb &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) - &
    (1.0d0-2.0d0*nu)*cotb/(Rbar+z3bar)*( &
        cosb*sinb - &
        1.0d0/2.0d0*a/R2bar*sinb*(2.0d0*y3+4.0d0*a)/cosb + &
        1.0d0/2.0d0*a*(Rbar*sinb-y1)/R3bar/cosb*(2.0d0*y3+4.0d0*a) &
    ) - &
    a/R3bar*y1*cotb + &
    3.0d0/2.0d0*a*y1*(y3+a)*cotb/R5bar*(2.0d0*y3+4.0d0*a) + &
    1.0d0/(Rbar+y3bar)*( &
        2.0d0*nu + &
        1.0d0/Rbar*((1.0d0-2.0d0*nu)*y1*cotb+a) - &
        y1**2/Rbar/(Rbar+y3bar)*(2.0d0*nu+a/Rbar) - &
        a*y1**2/R3bar &
    ) - &
    (y3+a)/(Rbar+y3bar)**2*( &
        2.0d0*nu + &
        1.0d0/Rbar*((1.0d0-2.0d0*nu)*y1*cotb+a) - &
        y1**2/Rbar/(Rbar+y3bar)*(2.0d0*nu+a/Rbar) - &
        a*y1**2/R3bar &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) + &
    (y3+a)/(Rbar+y3bar)*( &
        -1.0d0/2.0d0/R3bar*((1.0d0-2.0d0*nu)*y1*cotb+a)*(2.0d0*y3+4.0d0*a) + &
        1.0d0/2.0d0*y1**2/R3bar/(Rbar+y3bar)*(2.0d0*nu+a/Rbar)*(2.0d0*y3+4.0d0*a) + &
        y1**2/Rbar/(Rbar+y3bar)**2*(2.0d0*nu+a/Rbar)*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) + &
        1.0d0/2.0d0*y1**2/R4bar/(Rbar+y3bar)*a*(2.0d0*y3+4.0d0*a) + &
        3.0d0/2.0d0*a*y1**2/R5bar*(2.0d0*y3+4.0d0*a) &
    ) + &
    cotb/(Rbar+z3bar)*( &
        -cosb*sinb + &
        a*y1*y3bar/R3bar/cosb + &
        (Rbar*sinb-y1)/Rbar*( &
            (2.0d0-2.0d0*nu)*cosb - &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) &
        ) &
    ) - &
    (y3+a)*cotb/(Rbar+z3bar)**2*( &
        -cosb*sinb + &
        a*y1*y3bar/R3bar/cosb + &
        (Rbar*sinb-y1)/Rbar*( &
            (2.0d0-2.0d0*nu)*cosb - &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) &
        ) &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) + &
    (y3+a)*cotb/(Rbar+z3bar)*( &
        a/R3bar/cosb*y1 - &
        3.0d0/2.0d0*a*y1*y3bar/R5bar/cosb*(2.0d0*y3+4.0d0*a) + &
        1.0d0/2.0d0/R2bar*sinb*(2.0d0*y3+4.0d0*a)*( &
            (2.0d0-2.0d0*nu)*cosb - &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) &
        ) - &
        1.0d0/2.0d0*(Rbar*sinb-y1)/R3bar*( &
            (2.0d0-2.0d0*nu)*cosb - &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) &
        )*(2.0d0*y3+4.0d0*a) + &
        (Rbar*sinb-y1)/Rbar*( &
            -(1.0d0/2.0d0/Rbar*cosb*(2.0d0*y3+4.0d0*a)+1.0d0)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) + &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)**2*(1.0d0+a/Rbar/cosb)*( &
                1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a) + &
                cosb &
            ) + &
            1.0d0/2.0d0*(Rbar*cosb+y3bar)/(Rbar+z3bar)*a/R3bar/cosb*(2.0d0*y3+4.0d0*a) &
        ) &
    ) &
)/pi/(1.0d0-nu)

e13_termB2a = 1.0d0/2.0d0*(e13_termB2a_1+e13_termB2a_2)
return
end function

function e13_termB2b()
implicit none
double precision :: e13_termB2b, e13_termB2b_1, e13_termB2b_2
e13_termB2b_1 = 1.0d0/8.0d0*( &
    (-1.0d0+2.0d0*nu)*sinb*((1.0d0/R*y1-sinb)/(R-z3)-(1.0d0/Rbar*y1-sinb)/(Rbar+z3bar)) - &
    1.0d0/R + &
    1.0d0/Rbar - &
    y1*(-1.0d0/R3*y1+1.0d0/R3bar*y1) + &
    cosb*(R*cosb-y3)/R/(R-z3) + &
    z1/R2*cosb*y1/(R-z3) - &
    z1*(R*cosb-y3)/R3/(R-z3)*y1 - &
    z1*(R*cosb-y3)/R/(R-z3)**2*(1.0d0/R*y1-sinb) - &
    cosb*(Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar) - &
    z1bar/R2bar*cosb*y1/(Rbar+z3bar) + &
    z1bar*(Rbar*cosb+y3bar)/R3bar/(Rbar+z3bar)*y1 + &
    z1bar*(Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar)**2*(1.0d0/Rbar*y1-sinb) &
)/pi/(1.0d0-nu)

e13_termB2b_2 = 1.0d0/4.0d0*( &
    (-2.0d0+2.0d0*nu)*(1.0d0-2.0d0*nu)*cotb*( &
        1.0d0/Rbar*y1/(Rbar+y3bar) - &
        cosb*(1.0d0/Rbar*y1-sinb)/(Rbar+z3bar) &
    ) - &
    (2.0d0-2.0d0*nu)/(Rbar+y3bar)*(2.0d0*nu+a/Rbar) + &
    (2.0d0-2.0d0*nu)*y1**2/(Rbar+y3bar)**2*(2.0d0*nu+a/Rbar)/Rbar + &
    (2.0d0-2.0d0*nu)*y1**2/(Rbar+y3bar)*a/R3bar + &
    (2.0d0-2.0d0*nu)*cosb/(Rbar+z3bar)*(cosb+a/Rbar) - &
    (2.0d0-2.0d0*nu)*z1bar/(Rbar+z3bar)**2*(cosb+a/Rbar)*(1.0d0/Rbar*y1-sinb) - &
    (2.0d0-2.0d0*nu)*z1bar/(Rbar+z3bar)*a/R3bar*y1 - &
    (y3+a)/R3bar*( &
        (1.0d0-2.0d0*nu)*cotb - &
        2.0d0*nu*y1/(Rbar+y3bar) - &
        a*y1/R2bar &
    )*y1 + &
    (y3+a)/Rbar*( &
        -2.0d0*nu/(Rbar+y3bar) + &
        2.0d0*nu*y1**2/(Rbar+y3bar)**2/Rbar - &
        a/R2bar + &
        2.0d0*a*y1**2/R4bar &
    ) + &
    (y3+a)/(Rbar+z3bar)**2*( &
        cosb*sinb + &
        (Rbar*cosb+y3bar)*cotb/Rbar*((2.0d0-2.0d0*nu)*cosb-(Rbar*cosb+y3bar)/(Rbar+z3bar)) + &
        a/Rbar*(sinb-y3bar*z1bar/R2bar-z1bar*(Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar)) &
    )*(1.0d0/Rbar*y1-sinb) - &
    (y3+a)/(Rbar+z3bar)*( &
        1.0d0/R2bar*cosb*y1*cotb*((2.0d0-2.0d0*nu)*cosb-(Rbar*cosb+y3bar)/(Rbar+z3bar)) - &
        (Rbar*cosb+y3bar)*cotb/R3bar*((2.0d0-2.0d0*nu)*cosb-(Rbar*cosb+y3bar)/(Rbar+z3bar))*y1 + &
        (Rbar*cosb+y3bar)*cotb/Rbar*( &
            -1.0d0/Rbar*cosb*y1/(Rbar+z3bar) + &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)**2*(1.0d0/Rbar*y1-sinb) &
        ) - &
        a/R3bar*(sinb-y3bar*z1bar/R2bar-z1bar*(Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar))*y1 + &
        a/Rbar*( &
            -y3bar*cosb/R2bar + &
            2.0d0*y3bar*z1bar/R4bar*y1 - &
            cosb*(Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar) - &
            z1bar/R2bar*cosb*y1/(Rbar+z3bar) + &
            z1bar*(Rbar*cosb+y3bar)/R3bar/(Rbar+z3bar)*y1 + &
            z1bar*(Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar)**2*(1.0d0/Rbar*y1-sinb) &
        ) &
    ) &
)/pi/(1.0d0-nu)

e13_termB2b = 1.0d0/2.0d0*(e13_termB2b_1+e13_termB2b_2)
return
end function

function e13_termB3a()
implicit none
double precision :: e13_termB3a, e13_termB3a_1, e13_termB3a_2
e13_termB3a_1 = 1.0d0/8.0d0*y2*sinb*( &
    1.0d0/R2*sinb*y3/(R-z3) - &
    (R*sinb-y1)/R3/(R-z3)*y3 - &
    (R*sinb-y1)/R/(R-z3)**2*(1.0d0/R*y3-cosb) + &
    1.0d0/2.0d0/R2bar*sinb*(2.0d0*y3+4.0d0*a)/(Rbar+z3bar) - &
    1.0d0/2.0d0*(Rbar*sinb-y1)/R3bar/(Rbar+z3bar)*(2.0d0*y3+4.0d0*a) - &
    (Rbar*sinb-y1)/Rbar/(Rbar+z3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) &
)/pi/(1.0d0-nu)

e13_termB3a_2 = 1.0d0/4.0d0*( &
    (1.0d0-2.0d0*nu)*( &
        -y2/(Rbar+y3bar)**2*(1.0d0+a/Rbar)*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) - &
        1.0d0/2.0d0*y2/(Rbar+y3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) + &
        y2*cosb/(Rbar+z3bar)**2*(cosb+a/Rbar)*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) + &
        1.0d0/2.0d0*y2*cosb/(Rbar+z3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) &
    ) - &
    y2/Rbar*(a/R2bar+1.0d0/(Rbar+y3bar)) + &
    1.0d0/2.0d0*y2*(y3+a)/R3bar*(a/R2bar+1.0d0/(Rbar+y3bar))*(2.0d0*y3+4.0d0*a) - &
    y2*(y3+a)/Rbar*( &
       -a/R4bar*(2.0d0*y3+4.0d0*a) - &
       1.0d0/(Rbar+y3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) &
    ) + &
    y2*cosb/Rbar/(Rbar+z3bar)*( &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) + &
        a*y3bar/R2bar &
    ) - &
    1.0d0/2.0d0*y2*(y3+a)*cosb/R3bar/(Rbar+z3bar)*( &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) + &
        a*y3bar/R2bar &
    )*(2.0d0*y3+4.0d0*a) - &
    y2*(y3+a)*cosb/Rbar/(Rbar+z3bar)**2*( &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) + &
        a*y3bar/R2bar &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) + &
    y2*(y3+a)*cosb/Rbar/(Rbar+z3bar)*( &
        (1.0d0/2.0d0/Rbar*cosb*(2.0d0*y3+4.0d0*a)+1.0d0)/(Rbar+z3bar)*(cosb+a/Rbar) - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)**2*(cosb+a/Rbar)*( &
            1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a) + &
            cosb &
        ) - &
        1.0d0/2.0d0*(Rbar*cosb+y3bar)/(Rbar+z3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) + &
        a/R2bar - &
        a*y3bar/R4bar*(2.0d0*y3+4.0d0*a) &
    ) &
)/pi/(1.0d0-nu)

e13_termB3a = 1.0d0/2.0d0*(e13_termB3a_1+e13_termB3a_2)
return
end function

function e13_termB3b()
implicit none
double precision :: e13_termB3b, e13_termB3b_1, e13_termB3b_2
e13_termB3b_1 = 1.0d0/8.0d0*( &
    (2.0d0-2.0d0*nu)*( &
        -y2/z1**2*cosb/(1.0d0+y2**2/z1**2) + &
        ( &
            y2/R*sinb/(y1*z1+y2**2*cosb)*y1 - &
            y2*R*sinb/(y1*z1+y2**2*cosb)**2*(2.0d0*y1*cosb-y3*sinb) &
        )/(1.0d0+y2**2*R2*sinb**2/(y1*z1+y2**2*cosb)**2) + &
        y2/z1bar**2*cosb/(1.0d0+y2**2/z1bar**2) - &
        ( &
            y2/Rbar*sinb/(y1*z1bar+y2**2*cosb)*y1 - &
            y2*Rbar*sinb/(y1*z1bar+y2**2*cosb)**2*(2.0d0*y1*cosb+y3bar*sinb) &
        )/(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
    ) + &
    y2*sinb*( &
        1.0d0/R2*cosb*y1/(R-z3) - &
        (R*cosb-y3)/R3/(R-z3)*y1 - &
        (R*cosb-y3)/R/(R-z3)**2*(1.0d0/R*y1-sinb) - &
        1.0d0/R2bar*cosb*y1/(Rbar+z3bar) + &
        (Rbar*cosb+y3bar)/R3bar/(Rbar+z3bar)*y1 + &
        (Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar)**2*(1.0d0/Rbar*y1-sinb) &
    ) &
)/pi/(1.0d0-nu)

e13_termB3b_2 = 1.0d0/4.0d0*( &
    (2.0d0-2.0d0*nu)*( &
        y2/y1**2/(1.0d0+y2**2/y1**2) - &
        y2/z1bar**2*cosb/(1.0d0+y2**2/z1bar**2) + &
        ( &
            y2/Rbar*sinb/(y1*z1bar+y2**2*cosb)*y1 - &
            y2*Rbar*sinb/(y1*z1bar+y2**2*cosb)**2*(2.0d0*y1*cosb+y3bar*sinb) &
        )/(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
    ) - &
    (2.0d0-2.0d0*nu)*y2*sinb/(Rbar+z3bar)**2*(cosb+a/Rbar)*(1.0d0/Rbar*y1-sinb) - &
    (2.0d0-2.0d0*nu)*y2*sinb/(Rbar+z3bar)*a/R3bar*y1 - &
    y2*(y3+a)*sinb/R3bar/(Rbar+z3bar)*( &
        1.0d0 + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) + &
        a*y3bar/R2bar &
    )*y1 - &
    y2*(y3+a)*sinb/Rbar/(Rbar+z3bar)**2*( &
        1.0d0 + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) + &
        a*y3bar/R2bar &
    )*(1.0d0/Rbar*y1-sinb) + &
    y2*(y3+a)*sinb/Rbar/(Rbar+z3bar)*( &
        1.0d0/Rbar*cosb*y1/(Rbar+z3bar)*(cosb+a/Rbar) - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)**2*(cosb+a/Rbar)*(1.0d0/Rbar*y1-sinb) - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*a/R3bar*y1 - &
        2.0d0*a*y3bar/R4bar*y1 &
    ) &
)/pi/(1.0d0-nu)

e13_termB3b = 1.0d0/2.0d0*(e13_termB3b_1+e13_termB3b_2)
return
end function

function e23_termB1a()
implicit none
double precision :: e23_termB1a, e23_termB1a_1, e23_termB1a_2
e23_termB1a_1 = 1.0d0/8.0d0*( &
    (1.0d0-2.0d0*nu)*( &
        (1.0d0/R*y3-1.0d0)/(R-y3) + &
        (1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0)/(Rbar+y3bar) - &
        cosb*((1.0d0/R*y3-cosb)/(R-z3)+(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb)/(Rbar+z3bar)) &
    ) - &
    y2**2*( &
        -1.0d0/R3/(R-y3)*y3 - &
        1.0d0/R/(R-y3)**2*(1.0d0/R*y3-1.0d0) - &
        1.0d0/2.0d0/R3bar/(Rbar+y3bar)*(2.0d0*y3+4.0d0*a) - &
        1.0d0/Rbar/(Rbar+y3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) - &
        cosb*( &
            -1.0d0/R3/(R-z3)*y3 - &
            1.0d0/R/(R-z3)**2*(1.0d0/R*y3-cosb) - &
            1.0d0/2.0d0/R3bar/(Rbar+z3bar)*(2.0d0*y3+4.0d0*a) - &
            1.0d0/Rbar/(Rbar+z3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) &
        ) &
    ) &
)/pi/(1.0d0-nu)

e23_termB1a_2 = 1.0d0/4.0d0*( &
    (1.0d0-2.0d0*nu)*( &
        ((2.0d0-2.0d0*nu)*cotb**2-nu)*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0)/(Rbar+y3bar) - &
        ((2.0d0-2.0d0*nu)*cotb**2+1.0d0-2.0d0*nu)*cosb*( &
            1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a) + &
            cosb &
        )/(Rbar+z3bar) &
    ) + &
    (1.0d0-2.0d0*nu)/(Rbar+y3bar)**2*( &
        y1*cotb*(1.0d0-2.0d0*nu-a/Rbar)+nu*y3bar-a+y2**2/(Rbar+y3bar)*(nu+a/Rbar) &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) - &
    (1.0d0-2.0d0*nu)/(Rbar+y3bar)*( &
        1.0d0/2.0d0*a*y1*cotb/R3bar*(2.0d0*y3+4.0d0*a)+nu - &
        y2**2/(Rbar+y3bar)**2*(nu+a/Rbar)*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) - &
        1.0d0/2.0d0*y2**2/(Rbar+y3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) &
    ) - &
    (1.0d0-2.0d0*nu)*sinb*cotb/(Rbar+z3bar)*(cosb+a/Rbar) + &
    (1.0d0-2.0d0*nu)*z1bar*cotb/(Rbar+z3bar)**2*(cosb+a/Rbar)*( &
        1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a) + &
        cosb &
    ) + &
    1.0d0/2.0d0*(1.0d0-2.0d0*nu)*z1bar*cotb/(Rbar+z3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) - &
    a/R3bar*y1*cotb + &
    3.0d0/2.0d0*a*y1*(y3+a)*cotb/R5bar*(2.0d0*y3+4.0d0*a) + &
    1.0d0/(Rbar+y3bar)*( &
        -2.0d0*nu+1.0d0/Rbar*((1.0d0-2.0d0*nu)*y1*cotb-a) + &
        y2**2/Rbar/(Rbar+y3bar)*(2.0d0*nu+a/Rbar) + &
        a*y2**2/R3bar &
    ) - &
    (y3+a)/(Rbar+y3bar)**2*( &
        -2.0d0*nu + &
        1.0d0/Rbar*((1.0d0-2.0d0*nu)*y1*cotb-a) + &
        y2**2/Rbar/(Rbar+y3bar)*(2*nu+a/Rbar) + &
        a*y2**2/R3bar &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) + &
    (y3+a)/(Rbar+y3bar)*( &
        -1.0d0/2.0d0/R3bar*((1.0d0-2.0d0*nu)*y1*cotb-a)*(2.0d0*y3+4.0d0*a) - &
        1.0d0/2.0d0*y2**2/R3bar/(Rbar+y3bar)*(2.0d0*nu+a/Rbar)*(2.0d0*y3+4.0d0*a) - &
        y2**2/Rbar/(Rbar+y3bar)**2*(2.0d0*nu+a/Rbar)*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) - &
        1.0d0/2.0d0*y2**2/R4bar/(Rbar+y3bar)*a*(2.0d0*y3+4.0d0*a) - &
        3.0d0/2.0d0*a*y2**2/R5bar*(2.0d0*y3+4.0d0*a) &
    ) + &
    1.0d0/(Rbar+z3bar)*( &
        cosb**2 - &
        1.0d0/Rbar*((1.0d0-2.0d0*nu)*z1bar*cotb+a*cosb) + &
        a*y3bar*z1bar*cotb/R3bar - &
        1.0d0/Rbar/(Rbar+z3bar)*(y2**2*cosb**2-a*z1bar*cotb/Rbar*(Rbar*cosb+y3bar)) &
    ) - &
    (y3+a)/(Rbar+z3bar)**2*( &
        cosb**2 - &
        1.0d0/Rbar*((1.0d0-2.0d0*nu)*z1bar*cotb+a*cosb) + &
        a*y3bar*z1bar*cotb/R3bar - &
        1.0d0/Rbar/(Rbar+z3bar)*(y2**2*cosb**2-a*z1bar*cotb/Rbar*(Rbar*cosb+y3bar)) &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) + &
    (y3+a)/(Rbar+z3bar)*( &
        1.0d0/2.0d0/R3bar*((1.0d0-2.0d0*nu)*z1bar*cotb+a*cosb)*(2.0d0*y3+4.0d0*a) - &
        1.0d0/Rbar*(1.0d0-2.0d0*nu)*sinb*cotb + &
        a*z1bar*cotb/R3bar + &
        a*y3bar*sinb*cotb/R3bar - &
        3.0d0/2.0d0*a*y3bar*z1bar*cotb/R5bar*(2*y3+4*a) + &
        1.0d0/2.0d0/R3bar/(Rbar+z3bar)*( &
            y2**2*cosb**2 - &
            a*z1bar*cotb/Rbar*(Rbar*cosb+y3bar) &
        )*(2.0d0*y3+4.0d0*a) + &
        1.0d0/Rbar/(Rbar+z3bar)**2*( &
            y2**2*cosb**2 - &
            a*z1bar*cotb/Rbar*(Rbar*cosb+y3bar) &
        )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) - &
        1.0d0/Rbar/(Rbar+z3bar)*( &
            -a*sinb*cotb/Rbar*(Rbar*cosb+y3bar) + &
            1.0d0/2.0d0*a*z1bar*cotb/R3bar*(Rbar*cosb+y3bar)*(2.0d0*y3+4.0d0*a) - &
            a*z1bar*cotb/Rbar*(1.0d0/2.0d0/Rbar*cosb*(2.0d0*y3+4.0d0*a)+1.0d0) &
        ) &
    ) &
)/pi/(1.0d0-nu)

e23_termB1a = 1.0d0/2.0d0*(e23_termB1a_1+e23_termB1a_2)
return
end function

function e23_termB1b()
implicit none
double precision :: e23_termB1b, e23_termB1b_1, e23_termB1b_2
e23_termB1b_1 = 1.0d0/8.0d0*( &
    1.0d0/R - &
    1.0d0/Rbar - &
    cosb*((R*cosb-y3)/R/(R-z3)-(Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar)) + &
    y2*( &
        -1/R3*y2 + &
        1.0d0/R3bar*y2 - &
        cosb*( &
            1.0d0/R2*cosb*y2/(R-z3) - &
            (R*cosb-y3)/R3/(R-z3)*y2 - &
            (R*cosb-y3)/R2/(R-z3)**2*y2 - &
            1.0d0/R2bar*cosb*y2/(Rbar+z3bar) + &
            (Rbar*cosb+y3bar)/R3bar/(Rbar+z3bar)*y2 + &
            (Rbar*cosb+y3bar)/R2bar/(Rbar+z3bar)**2*y2 &
        ) &
    ) &
)/pi/(1.0d0-nu)

e23_termB1b_2 = 1.0d0/4.0d0*( &
    (2.0d0-2.0d0*nu)*( &
        (1.0d0-2.0d0*nu)*( &
            -1/y1/(1.0d0+y2**2/y1**2) + &
            1.0d0/z1bar/(1.0d0+y2**2/z1bar**2) + &
            (Rbar*sinb/(y1*z1bar+y2**2*cosb) + &
            y2**2/Rbar*sinb/(y1*z1bar+y2**2*cosb) - &
            2.0d0*y2**2*Rbar*sinb/(y1*z1bar+y2**2*cosb)**2*cosb)/( &
                1 + &
                y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2 &
            ) &
        )*cotb + &
        1.0d0/(Rbar+y3bar)*(2*nu+a/Rbar) - &
        y2**2/(Rbar+y3bar)**2*(2*nu+a/Rbar)/Rbar - &
        y2**2/(Rbar+y3bar)*a/R3bar - &
        cosb/(Rbar+z3bar)*(cosb+a/Rbar) + &
        y2**2*cosb/(Rbar+z3bar)**2*(cosb+a/Rbar)/Rbar + &
        y2**2*cosb/(Rbar+z3bar)*a/R3bar &
    ) + &
    (y3+a)/Rbar*(2*nu/(Rbar+y3bar)+a/R2bar) - &
    y2**2*(y3+a)/R3bar*(2*nu/(Rbar+y3bar)+a/R2bar) + &
    y2*(y3+a)/Rbar*(-2.0d0*nu/(Rbar+y3bar)**2/Rbar*y2-2.0d0*a/R4bar*y2) + &
    (y3+a)*cosb/Rbar/(Rbar+z3bar)*( &
        1 - &
        2.0d0*nu - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) - &
        a*y3bar/R2bar &
    ) - &
    y2**2*(y3+a)*cosb/R3bar/(Rbar+z3bar)*( &
        1 - &
        2.0d0*nu - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) - &
        a*y3bar/R2bar &
    ) - &
    y2**2*(y3+a)*cosb/R2bar/(Rbar+z3bar)**2*( &
        1 - &
        2.0d0*nu - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) - &
        a*y3bar/R2bar &
    ) + &
    y2*(y3+a)*cosb/Rbar/(Rbar+z3bar)*( &
        -1/Rbar*cosb*y2/(Rbar+z3bar)*(cosb+a/Rbar) + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)**2*(cosb+a/Rbar)/Rbar*y2 + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*a/R3bar*y2 + &
        2.0d0*a*y3bar/R4bar*y2 &
    ) &
)/pi/(1.0d0-nu)

e23_termB1b = 1.0d0/2.0d0*(e23_termB1b_1+e23_termB1b_2)
return
end function

function e23_termB2a()
implicit none
double precision :: e23_termB2a, e23_termB2a_1, e23_termB2a_2
e23_termB2a_1 = 1.0d0/8.0d0*( &
    (2.0d0-2.0d0*nu)*( &
        y2/z1**2*sinb/(1.0d0+y2**2/z1**2) + &
        ( &
            y2/R*sinb/(y1*z1+y2**2*cosb)*y3 + &
            y2*R*sinb**2/(y1*z1+y2**2*cosb)**2*y1 &
        )/(1.0d0+y2**2*R2*sinb**2/(y1*z1+y2**2*cosb)**2) - &
        y2/z1bar**2*sinb/(1.0d0+y2**2/z1bar**2) + &
        ( &
            1.0d0/2.0d0*y2/Rbar*sinb/(y1*z1bar+y2**2*cosb)*(2.0d0*y3+4.0d0*a) - &
            y2*Rbar*sinb**2/(y1*z1bar+y2**2*cosb)**2*y1 &
        )/(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
    ) + &
    y1*y2*( &
        -1.0d0/R3/(R-y3)*y3 - &
        1.0d0/R/(R-y3)**2*(1.0d0/R*y3-1.0d0) - &
        1.0d0/2.0d0/R3bar/(Rbar+y3bar)*(2.0d0*y3+4.0d0*a) - &
        1.0d0/Rbar/(Rbar+y3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) &
    ) - &
    y2*( &
        -sinb/R/(R-z3) - &
        z1/R3/(R-z3)*y3 - &
        z1/R/(R-z3)**2*(1.0d0/R*y3-cosb) + &
        sinb/Rbar/(Rbar+z3bar) - &
        1.0d0/2.0d0*z1bar/R3bar/(Rbar+z3bar)*(2.0d0*y3+4.0d0*a) - &
        z1bar/Rbar/(Rbar+z3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) &
    ) &
)/pi/(1.0d0-nu)

e23_termB2a_2 = 1.0d0/4.0d0*( &
    (2.0d0-2.0d0*nu)*(1.0d0-2.0d0*nu)*( &
        -y2/z1bar**2*sinb/(1.0d0+y2**2/z1bar**2) + &
        ( &
            1.0d0/2.0d0*y2/Rbar*sinb/(y1*z1bar+y2**2*cosb)*(2.0d0*y3+4.0d0*a) - &
            y2*Rbar*sinb**2/(y1*z1bar+y2**2*cosb)**2*y1 &
        )/(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
    )*cotb**2 - &
    (1.0d0-2.0d0*nu)*y2/(Rbar+y3bar)**2*( &
        (-1.0d0+2.0d0*nu+a/Rbar)*cotb + &
        y1/(Rbar+y3bar)*(nu+a/Rbar) &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) + &
    (1.0d0-2.0d0*nu)*y2/(Rbar+y3bar)*( &
        -1.0d0/2.0d0*a/R3bar*(2.0d0*y3+4.0d0*a)*cotb - &
        y1/(Rbar+y3bar)**2*(nu+a/Rbar)*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) - &
        1.0d0/2.0d0*y1/(Rbar+y3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) &
    ) + &
    (1.0d0-2.0d0*nu)*y2*cotb/(Rbar+z3bar)**2*(1.0d0+a/Rbar/cosb)*( &
        1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a) + &
        cosb &
    ) + &
    1.0d0/2.0d0*(1.0d0-2.0d0*nu)*y2*cotb/(Rbar+z3bar)*a/R3bar/cosb*(2.0d0*y3+4.0d0*a) - &
    a/R3bar*y2*cotb + &
    3.0d0/2.0d0*a*y2*(y3+a)*cotb/R5bar*(2.0d0*y3+4.0d0*a) + &
    y2/Rbar/(Rbar+y3bar)*( &
        (1.0d0-2.0d0*nu)*cotb - &
        2.0d0*nu*y1/(Rbar+y3bar) - &
        a*y1/Rbar*(1.0d0/Rbar+1.0d0/(Rbar+y3bar)) &
    ) - &
    1.0d0/2.0d0*y2*(y3+a)/R3bar/(Rbar+y3bar)*( &
        (1.0d0-2.0d0*nu)*cotb - &
        2.0d0*nu*y1/(Rbar+y3bar) - &
        a*y1/Rbar*(1.0d0/Rbar+1.0d0/(Rbar+y3bar)) &
    )*(2.0d0*y3+4.0d0*a) - &
    y2*(y3+a)/Rbar/(Rbar+y3bar)**2*( &
        (1.0d0-2.0d0*nu)*cotb - &
        2.0d0*nu*y1/(Rbar+y3bar) - &
        a*y1/Rbar*(1.0d0/Rbar+1.0d0/(Rbar+y3bar)) &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) + &
    y2*(y3+a)/Rbar/(Rbar+y3bar)*( &
        2.0d0*nu*y1/(Rbar+y3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) + &
        1.0d0/2.0d0*a*y1/R3bar*(1.0d0/Rbar+1.0d0/(Rbar+y3bar))*(2.0d0*y3+4.0d0*a) - &
        a*y1/Rbar*( &
            -1.0d0/2.0d0/R3bar*(2.0d0*y3+4.0d0*a) - &
            1.0d0/(Rbar+y3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) &
        ) &
    ) + &
    y2*cotb/Rbar/(Rbar+z3bar)*( &
        (-2.0d0+2.0d0*nu)*cosb + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) + &
        a*y3bar/R2bar/cosb &
    ) - &
    1.0d0/2.0d0*y2*(y3+a)*cotb/R3bar/(Rbar+z3bar)*( &
        (-2.0d0+2.0d0*nu)*cosb + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) + &
        a*y3bar/R2bar/cosb &
    )*(2.0d0*y3+4.0d0*a) - &
    y2*(y3+a)*cotb/Rbar/(Rbar+z3bar)**2*( &
        (-2.0d0+2.0d0*nu)*cosb + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) + &
        a*y3bar/R2bar/cosb &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) + &
    y2*(y3+a)*cotb/Rbar/(Rbar+z3bar)*( &
        (1.0d0/2.0d0/Rbar*cosb*(2.0d0*y3+4.0d0*a)+1.0d0)/(Rbar+z3bar)*(1.0d0+a/Rbar/cosb) - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)**2*(1.0d0+a/Rbar/cosb)*( &
            1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a) + &
            cosb &
        ) - &
        1.0d0/2.0d0*(Rbar*cosb+y3bar)/(Rbar+z3bar)*a/R3bar/cosb*(2.0d0*y3+4.0d0*a) + &
        a/R2bar/cosb - &
        a*y3bar/R4bar/cosb*(2.0d0*y3+4.0d0*a) &
    ) &
)/pi/(1.0d0-nu)

e23_termB2a = 1.0d0/2.0d0*(e23_termB2a_1+e23_termB2a_2)
return
end function

function e23_termB2b()
implicit none
double precision :: e23_termB2b, e23_termB2b_1, e23_termB2b_2
e23_termB2b_1 = 1.0d0/8.0d0*( &
    (-1.0d0+2.0d0*nu)*sinb*(1.0d0/R*y2/(R-z3) - &
    1.0d0/Rbar*y2/(Rbar+z3bar)) - &
    y1*(-1.0d0/R3*y2+1.0d0/R3bar*y2) + &
    z1/R2*cosb*y2/(R-z3) - &
    z1*(R*cosb-y3)/R3/(R-z3)*y2 - &
    z1*(R*cosb-y3)/R2/(R-z3)**2*y2 - &
    z1bar/R2bar*cosb*y2/(Rbar+z3bar) + &
    z1bar*(Rbar*cosb+y3bar)/R3bar/(Rbar+z3bar)*y2 + &
    z1bar*(Rbar*cosb+y3bar)/R2bar/(Rbar+z3bar)**2*y2 &
)/pi/(1.0d0-nu)

e23_termB2b_2 = 1.0d0/4.0d0*( &
    (-2.0d0+2.0d0*nu)*(1.0d0-2.0d0*nu)*cotb*( &
        1.0d0/Rbar*y2/(Rbar+y3bar) - &
        cosb/Rbar*y2/(Rbar+z3bar) &
    ) + &
    (2.0d0-2.0d0*nu)*y1/(Rbar+y3bar)**2*(2.0d0*nu+a/Rbar)/Rbar*y2 + &
    (2.0d0-2.0d0*nu)*y1/(Rbar+y3bar)*a/R3bar*y2 - &
    (2.0d0-2.0d0*nu)*z1bar/(Rbar+z3bar)**2*(cosb+a/Rbar)/Rbar*y2 - &
    (2.0d0-2.0d0*nu)*z1bar/(Rbar+z3bar)*a/R3bar*y2 - &
    (y3+a)/R3bar*((1.0d0-2.0d0*nu)*cotb-2.0d0*nu*y1/(Rbar+y3bar)-a*y1/R2bar)*y2 + &
    (y3+a)/Rbar*(2.0d0*nu*y1/(Rbar+y3bar)**2/Rbar*y2+2.0d0*a*y1/R4bar*y2) + &
    (y3+a)/(Rbar+z3bar)**2*( &
        cosb*sinb + &
        (Rbar*cosb+y3bar)*cotb/Rbar*((2.0d0-2.0d0*nu)*cosb-(Rbar*cosb+y3bar)/(Rbar+z3bar)) + &
        a/Rbar*(sinb-y3bar*z1bar/R2bar-z1bar*(Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar)) &
    )/Rbar*y2 - &
    (y3+a)/(Rbar+z3bar)*( &
        1.0d0/R2bar*cosb*y2*cotb*((2.0d0-2.0d0*nu)*cosb-(Rbar*cosb+y3bar)/(Rbar+z3bar)) - &
        (Rbar*cosb+y3bar)*cotb/R3bar*((2.0d0-2.0d0*nu)*cosb-(Rbar*cosb+y3bar)/(Rbar+z3bar))*y2 + &
        (Rbar*cosb+y3bar)*cotb/Rbar*( &
            -cosb/Rbar*y2/(Rbar+z3bar) + &
            (Rbar*cosb+y3bar)/(Rbar+z3bar)**2/Rbar*y2 &
        ) - &
        a/R3bar*(sinb-y3bar*z1bar/R2bar-z1bar*(Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar))*y2 + &
        a/Rbar*( &
            2.0d0*y3bar*z1bar/R4bar*y2 - &
            z1bar/R2bar*cosb*y2/(Rbar+z3bar) + &
            z1bar*(Rbar*cosb+y3bar)/R3bar/(Rbar+z3bar)*y2 + &
            z1bar*(Rbar*cosb+y3bar)/R2bar/(Rbar+z3bar)**2*y2 &
        ) &
    ) &
)/pi/(1.0d0-nu)

e23_termB2b = 1.0d0/2.0d0*(e23_termB2b_1+e23_termB2b_2)
return
end function

function e23_termB3a()
implicit none
double precision :: e23_termB3a, e23_termB3a_1, e23_termB3a_2
e23_termB3a_1 = 1.0d0/8.0d0*( &
    (1.0d0-2.0d0*nu)*sinb*( &
        (1.0d0/R*y3-cosb)/(R-z3) + &
        (1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb)/(Rbar+z3bar) &
    ) - &
    y2**2*sinb*( &
        -1.0d0/R3/(R-z3)*y3 - &
        1.0d0/R/(R-z3)**2*(1.0d0/R*y3-cosb) - &
        1.0d0/2.0d0/R3bar/(Rbar+z3bar)*(2.0d0*y3+4.0d0*a) - &
        1.0d0/Rbar/(Rbar+z3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) &
    ) &
)/pi/(1.0d0-nu)

e23_termB3a_2 = 1.0d0/4.0d0*( &
    (1.0d0-2.0d0*nu)*( &
        -sinb*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb)/(Rbar+z3bar) + &
        y1/(Rbar+y3bar)**2*(1.0d0+a/Rbar)*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) + &
        1.0d0/2.0d0*y1/(Rbar+y3bar)*a/R3bar*(2.0d0*y3+4.0d0*a)+sinb/(Rbar+z3bar)*(cosb+a/Rbar) - &
        z1bar/(Rbar+z3bar)**2*(cosb+a/Rbar)*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) - &
        1.0d0/2.0d0*z1bar/(Rbar+z3bar)*a/R3bar*(2.0d0*y3+4.0d0*a) &
    ) + &
    y1/Rbar*(a/R2bar+1.0d0/(Rbar+y3bar)) - &
    1.0d0/2.0d0*y1*(y3+a)/R3bar*(a/R2bar+1.0d0/(Rbar+y3bar))*(2.0d0*y3+4.0d0*a) + &
    y1*(y3+a)/Rbar*( &
        -a/R4bar*(2.0d0*y3+4.0d0*a) - &
        1.0d0/(Rbar+y3bar)**2*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+1.0d0) &
    ) - &
    1.0d0/(Rbar+z3bar)*( &
        sinb*(cosb-a/Rbar) + &
        z1bar/Rbar*(1.0d0+a*y3bar/R2bar) - &
        1.0d0/Rbar/(Rbar+z3bar)*(y2**2*cosb*sinb-a*z1bar/Rbar*(Rbar*cosb+y3bar)) &
    ) + &
    (y3+a)/(Rbar+z3bar)**2*( &
        sinb*(cosb-a/Rbar)+z1bar/Rbar*(1.0d0+a*y3bar/R2bar) - &
        1.0d0/Rbar/(Rbar+z3bar)*(y2**2*cosb*sinb-a*z1bar/Rbar*(Rbar*cosb+y3bar)) &
    )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) - &
    (y3+a)/(Rbar+z3bar)*( &
        1.0d0/2.0d0*sinb*a/R3bar*(2.0d0*y3+4.0d0*a) + &
        sinb/Rbar*(1.0d0+a*y3bar/R2bar) - &
        1.0d0/2.0d0*z1bar/R3bar*(1.0d0+a*y3bar/R2bar)*(2.0d0*y3+4.0d0*a) + &
        z1bar/Rbar*(a/R2bar-a*y3bar/R4bar*(2.0d0*y3+4.0d0*a)) + &
        1.0d0/2.0d0/R3bar/(Rbar+z3bar)*( &
            y2**2*cosb*sinb-a*z1bar/Rbar*(Rbar*cosb+y3bar) &
        )*(2.0d0*y3+4.0d0*a) + &
        1.0d0/Rbar/(Rbar+z3bar)**2*( &
            y2**2*cosb*sinb - &
            a*z1bar/Rbar*(Rbar*cosb+y3bar) &
        )*(1.0d0/2.0d0/Rbar*(2.0d0*y3+4.0d0*a)+cosb) - &
        1.0d0/Rbar/(Rbar+z3bar)*( &
            -a*sinb/Rbar*(Rbar*cosb+y3bar) + &
            1.0d0/2.0d0*a*z1bar/R3bar*(Rbar*cosb+y3bar)*(2.0d0*y3+4.0d0*a) - &
            a*z1bar/Rbar*(1.0d0/2.0d0/Rbar*cosb*(2.0d0*y3+4.0d0*a)+1.0d0) &
        ) &
    ) &
)/pi/(1.0d0-nu)

e23_termB3a = 1.0d0/2.0d0*(e23_termB3a_1+e23_termB3a_2)
return
end function

function e23_termB3b()
implicit none
double precision :: e23_termB3b, e23_termB3b_1, e23_termB3b_2
e23_termB3b_1 = 1.0d0/8.0d0*( &
    (2.0d0-2.0d0*nu)*( &
        1.0d0/z1/(1.0d0+y2**2/z1**2) + &
        ( &
            R*sinb/(y1*z1+y2**2*cosb) + &
            y2**2/R*sinb/(y1*z1+y2**2*cosb) - &
            2.0d0*y2**2*R*sinb/(y1*z1+y2**2*cosb)**2*cosb)/( &
                1.0d0 + &
                y2**2*R2*sinb**2/(y1*z1+y2**2*cosb)**2 &
            ) - &
            1.0d0/z1bar/(1.0d0+y2**2/z1bar**2) - &
            ( &
                Rbar*sinb/(y1*z1bar+y2**2*cosb) + &
                y2**2/Rbar*sinb/(y1*z1bar+y2**2*cosb) - &
                2.0d0*y2**2*Rbar*sinb/(y1*z1bar+y2**2*cosb)**2*cosb &
            )/(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
    ) + &
    sinb*((R*cosb-y3)/R/(R-z3)-(Rbar*cosb+y3bar)/Rbar/(Rbar+z3bar)) + &
    y2*sinb*( &
        1.0d0/R2*cosb*y2/(R-z3) - &
        (R*cosb-y3)/R3/(R-z3)*y2 - &
        (R*cosb-y3)/R2/(R-z3)**2*y2 - &
        1.0d0/R2bar*cosb*y2/(Rbar+z3bar) + &
        (Rbar*cosb+y3bar)/R3bar/(Rbar+z3bar)*y2 + &
        (Rbar*cosb+y3bar)/R2bar/(Rbar+z3bar)**2*y2 &
    ) &
)/pi/(1.0d0-nu)

e23_termB3b_2 = 1.0d0/4.0d0*( &
    (2.0d0-2.0d0*nu)*( &
        -1.0d0/y1/(1.0d0+y2**2/y1**2) + &
        1.0d0/z1bar/(1.0d0+y2**2/z1bar**2) + &
        ( &
            Rbar*sinb/(y1*z1bar+y2**2*cosb) + &
            y2**2/Rbar*sinb/(y1*z1bar+y2**2*cosb) - &
            2.0d0*y2**2*Rbar*sinb/(y1*z1bar+y2**2*cosb)**2*cosb &
        )/(1.0d0+y2**2*R2bar*sinb**2/(y1*z1bar+y2**2*cosb)**2) &
    ) + &
    (2.0d0-2.0d0*nu)*sinb/(Rbar+z3bar)*(cosb+a/Rbar) - &
    (2.0d0-2.0d0*nu)*y2**2*sinb/(Rbar+z3bar)**2*(cosb+a/Rbar)/Rbar - &
    (2.0d0-2.0d0*nu)*y2**2*sinb/(Rbar+z3bar)*a/R3bar + &
    (y3+a)*sinb/Rbar/(Rbar+z3bar)*( &
        1.0d0 + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) + &
        a*y3bar/R2bar &
    ) - &
    y2**2*(y3+a)*sinb/R3bar/(Rbar+z3bar)*( &
        1.0d0 + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) + &
        a*y3bar/R2bar &
    ) - &
    y2**2*(y3+a)*sinb/R2bar/(Rbar+z3bar)**2*( &
        1.0d0 + &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*(cosb+a/Rbar) + &
        a*y3bar/R2bar &
    ) + &
    y2*(y3+a)*sinb/Rbar/(Rbar+z3bar)*( &
        1.0d0/Rbar*cosb*y2/(Rbar+z3bar)*(cosb+a/Rbar) - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)**2*(cosb+a/Rbar)/Rbar*y2 - &
        (Rbar*cosb+y3bar)/(Rbar+z3bar)*a/R3bar*y2 - &
        2.0d0*a*y3bar/R4bar*y2 &
    ) &
)/pi/(1.0d0-nu)

e23_termB3b = 1.0d0/2.0d0*(e23_termB3b_1+e23_termB3b_2)
return
end function

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine tri_center(center,p1,p2,p3)
implicit none
double precision :: center(3), p1(3), p2(3), p3(3)
center(1) = (p1(1)+p2(1)+p3(1))/3.0d0
center(2) = (p1(2)+p2(2)+p3(2))/3.0d0
center(3) = (p1(3)+p2(3)+p3(3))/3.0d0
return
end

subroutine tri_geometry(normal,strike,updip,p1in,p2in,p3in)
!----
! Calculate the triangle normal vector, strike-parallel vector, and up-dip vector
! Assumes depths of points are positive down, but vectors have positive z up
!----

implicit none
! Arguments
double precision :: normal(3), strike(3), updip(3), p1in(3), p2in(3), p3in(3)
! Local variables
double precision :: p1(3), p2(3), p3(3), magnitude

p1 = p1in
p2 = p2in
p3 = p3in
p1(3) = -p1(3)
p2(3) = -p2(3)
p3(3) = -p3(3)

! Normal vector
normal(1) = (p2(2)-p1(2))*(p3(3)-p1(3)) - (p2(3)-p1(3))*(p3(2)-p1(2))
normal(2) = (p2(3)-p1(3))*(p3(1)-p1(1)) - (p2(1)-p1(1))*(p3(3)-p1(3))
normal(3) = (p2(1)-p1(1))*(p3(2)-p1(2)) - (p2(2)-p1(2))*(p3(1)-p1(1))
magnitude = normal(1)*normal(1)+normal(2)*normal(2)+normal(3)*normal(3)
magnitude = dsqrt(magnitude)
normal = normal/magnitude
! Normal must point upward
if (normal(3).lt.0.0d0) then
    normal = -normal
endif

! Strike-parallel vector
strike(1) = -dsin(datan2(normal(2),normal(1)))
strike(2) = dcos(datan2(normal(2),normal(1)))
strike(3) = 0.0d0

! Vector pointing up-dip
updip(1) = normal(2)*strike(3) - normal(3)*strike(2)
updip(2) = normal(3)*strike(1) - normal(1)*strike(3)
updip(3) = normal(1)*strike(2) - normal(2)*strike(1)

! write(0,*) 'tri_geometry: normal',normal
! write(0,*) 'tri_geometry: strike',strike
! write(0,*) 'tri_geometry: updip',updip

return
end subroutine

subroutine tri_geo2cart(pt1,pt2,pt3,pt1_in,pt2_in,pt3_in,depthUnits)
implicit none
! Arguments
double precision :: pt1(3),pt2(3),pt3(3),pt1_in(3),pt2_in(3),pt3_in(3)
character(len=*) :: depthUnits
! Local variables
double precision :: center(3), dist, az, radius

radius = 0.0d0
if (trim(depthUnits).eq.'m') then
    radius = 6.371d6
elseif (trim(depthUnits).eq.'km') then
    radius = 6.371d3
else
    write(0,*) 'tri_geo2cart: no depth unit named ',trim(depthUnits)
endif

call tri_center(center,pt1_in,pt2_in,pt3_in)

call ddistaz(dist,az,center(1),center(2),pt1_in(1),pt1_in(2))
dist = dist*radius
pt1(1) = dist*dsin(az)
pt1(2) = dist*dcos(az)
pt1(3) = pt1_in(3)

call ddistaz(dist,az,center(1),center(2),pt2_in(1),pt2_in(2))
dist = dist*radius
pt2(1) = dist*dsin(az)
pt2(2) = dist*dcos(az)
pt2(3) = pt2_in(3)

call ddistaz(dist,az,center(1),center(2),pt3_in(1),pt3_in(2))
dist = dist*radius
pt3(1) = dist*dsin(az)
pt3(2) = dist*dcos(az)
pt3(3) = pt3_in(3)

return
end subroutine

end module
