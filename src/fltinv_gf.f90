!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!-------------------------------- GREEN'S FUNCTIONS SUBROUTINES -----------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
!------------------------------------- DISPLACEMENT GFS -------------------------------------------!
!--------------------------------------------------------------------------------------------------!
subroutine calc_disp_gfs()
!----
! Calculate the Green's functions for three-component displacements (gf_disp) or line-of-sight
! displacements (gf_los), given the specified model (gf_model).
!----

use io, only: stdout, stderr, verbosity, debug, progress_indicator
use trig, only: d2r, r2d
use algebra, only: rotate_vector_angle_axis
use earth, only: radius_earth_m
use geom, only: lola2distaz
use okada92, only: o92_pt_disp, o92_rect_disp
use tri_disloc, only: tri_center, tri_disloc_disp

use fltinv, only: inversion_mode, &
                  fault, &
                  rake_constraint, &
                  displacement, &
                  los, &
                  gf_model, &
                  gf_disp, &
                  gf_los, &
                  poisson, &
                  shearmod, &
                  lame, &
                  coord_type

implicit none

! Local variables
integer :: i, j, ierr, iTri, ndsp, nlos, nflt
double precision :: slip_mag, rak1, rak2, slip(3), mom(4), disp1(3), disp2(3)
double precision :: sta(3), sta_new(3), dist, az, dx, dy
double precision :: evlo, evla, evdp, str, dip, wid, len, area, tri(3,4), tri_new(3,4), center(3)


if (verbosity.ge.1) then
    write(stdout,*) 'calc_disp_gfs: starting'
endif

if (displacement%file.eq.'none'.and.los%file.eq.'none') then
    if (verbosity.ge.1) then
        write(stdout,*) 'calc_disp_gfs: displacements are not used in this inversion'
        write(stdout,*)
    endif
    return
endif

if (gf_model.eq.'precomputed') then
    if (verbosity.ge.1) then
        write(stdout,*) 'calc_disp_gfs: displacement Greens functions already computed'
        write(stdout,*)
    endif
    return
endif

! Set unit slip magnitude for the Green's functions calculation
slip_mag = 1.0d0

ndsp = displacement%nrows
nlos = los%nrows
nflt = fault%nrows


! Calculate displacement Green's function for each fault-station pair
do i = 1,ndsp+nlos
    if (debug) then
        write(stderr,*) 'i:',i,' n:',ndsp+nlos
    endif

    ! Station coordinates
    if (i.le.ndsp) then
        sta(1) = displacement%array(i,1)
        sta(2) = displacement%array(i,2)
        sta(3) = displacement%array(i,3)
    else
        sta(1) = los%array(i-ndsp,1)
        sta(2) = los%array(i-ndsp,2)
        sta(3) = los%array(i-ndsp,3)
    endif
    sta_new(3) = sta(3)

    do j = 1,nflt
        if (debug) then
            write(stderr,*) 'j:',j,' n:',nflt
        endif

        ! Rake angle constraints
        if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
            ! In lsqr mode, rake constraints change the rake for computing Green's functions
            rak1 = rake_constraint%array(j,1)
            if (rake_constraint%ncols.eq.2) then
                rak2 = rake_constraint%array(j,2)
            else
                rak2 = 90.0d0
            endif
        else
            ! Otherwise, calculate Green's functions for strike-slip and dip-slip sources
            rak1 = 0.0d0
            rak2 = 90.0d0
        endif
        if (debug) then
            write(stderr,*) 'rak1:',rak1,' rak2:',rak2
        endif

        ! The details of calculating the Green's functions differ depending on the model, but
        ! the general steps are the same:
        !     1. Calculate the location of the station relative to the center of the fault
        !     2. Calculate the displacement for unit slip in both rake directions

        if (gf_model.eq.'okada_rect') then
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str = fault%array(j,4)
            dip = fault%array(j,5)
            wid = fault%array(j,6)
            len = fault%array(j,7)

            ! Calculate location of station relative to fault center
            if (coord_type.eq.'geographic') then
                call lola2distaz(evlo,evla,sta(1),sta(2),dist,az,'radians','degrees',ierr)
                if (ierr.ne.0) then
                    call usage('calc_disp_gfs: error calculating dist and az (usage:none)')
                endif
                dist = dist*radius_earth_m
            elseif (coord_type.eq.'cartesian') then
                dx = sta(1) - evlo
                dy = sta(2) - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)*r2d
            endif

            ! Rotate coordinate axes so that x=strike, y=horizontal updip
            sta_new(1) = dist*( dcos((az-str)*d2r))
            sta_new(2) = dist*(-dsin((az-str)*d2r))

            ! Calculate displacements for both rake angles
            slip(1) = slip_mag*cos(rak1*d2r)
            slip(2) = slip_mag*sin(rak1*d2r)
            slip(3) = 0.0d0
            call o92_rect_disp(disp1,sta_new,evdp,dip,slip,wid,len,lame,shearmod)
            slip(1) = slip_mag*cos(rak2*d2r)
            slip(2) = slip_mag*sin(rak2*d2r)
            slip(3) = 0.0d0
            call o92_rect_disp(disp2,sta_new,evdp,dip,slip,wid,len,lame,shearmod)

            ! Rotate displacements back to ENZ
            call rotate_vector_angle_axis(disp1,90.0d0-str,'z',disp1,ierr)
            call rotate_vector_angle_axis(disp2,90.0d0-str,'z',disp2,ierr)

        elseif (gf_model.eq.'okada_pt') then
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str = fault%array(j,4)
            dip = fault%array(j,5)
            area = fault%array(j,6)

            ! Calculate location of station relative to fault center
            if (coord_type.eq.'geographic') then
                call lola2distaz(evlo,evla,sta(1),sta(2),dist,az,'radians','degrees',ierr)
                if (ierr.ne.0) then
                    call usage('calc_disp_gfs: error calculating dist and az (usage:none)')
                endif
                dist = dist*radius_earth_m
            elseif (coord_type.eq.'cartesian') then
                dx = sta(1) - evlo
                dy = sta(2) - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)*r2d
            endif

            ! Rotate coordinate axes so that x=strike, y=horizontal updip
            sta_new(1) = dist*( dcos((az-str)*d2r))
            sta_new(2) = dist*(-dsin((az-str)*d2r))

            ! Calculate displacements for both rake angles
            mom(1) = slip_mag*cos(rak1*d2r)*area*shearmod
            mom(2) = slip_mag*sin(rak1*d2r)*area*shearmod
            mom(3:4) = 0.0d0
            call o92_pt_disp(disp1,sta_new,evdp,dip,mom,lame,shearmod)
            mom(1) = slip_mag*cos(rak2*d2r)*area*shearmod
            mom(2) = slip_mag*sin(rak2*d2r)*area*shearmod
            mom(3:4) = 0.0d0
            call o92_pt_disp(disp2,sta_new,evdp,dip,mom,lame,shearmod)

            ! Rotate displacements back to ENZ
            call rotate_vector_angle_axis(disp1,90.0d0-str,'z',disp1,ierr)
            call rotate_vector_angle_axis(disp2,90.0d0-str,'z',disp2,ierr)

        elseif (gf_model.eq.'triangle') then
            tri(1,1) = fault%array(j,1)
            tri(2,1) = fault%array(j,2)
            tri(3,1) = fault%array(j,3)
            tri(1,2) = fault%array(j,4)
            tri(2,2) = fault%array(j,5)
            tri(3,2) = fault%array(j,6)
            tri(1,3) = fault%array(j,7)
            tri(2,3) = fault%array(j,8)
            tri(3,3) = fault%array(j,9)
            tri(:,4) = 0.0d0

            ! Calculate location of station relative to fault center
            if (coord_type.eq.'geographic') then
                call tri_center(center,tri(:,1),tri(:,2),tri(:,3))
                do iTri = 1,3
                    call lola2distaz(center(1),center(2),tri(1,iTri),tri(2,iTri),dist,az, &
                                'radians','radians',ierr)
                    if (ierr.ne.0) then
                        call usage('calc_disp_gfs: error calculating dist and az (usage:none)')
                    endif
                    dist = dist*radius_earth_m
                    tri_new(1,iTri) = dist*dsin(az)
                    tri_new(2,iTri) = dist*dcos(az)
                    tri_new(3,iTri) = tri(3,iTri)
                enddo
                call lola2distaz(center(1),center(2),sta(1),sta(2),dist,az,'radians','radians',ierr)
                if (ierr.ne.0) then
                    call usage('calc_disp_gfs: error calculating dist and az (usage:none)')
                endif
                dist = dist*radius_earth_m
                sta_new(1) = dist*dsin(az)
                sta_new(2) = dist*dcos(az)
            elseif (coord_type.eq.'cartesian') then
                tri_new = tri
                sta_new = sta
            endif

            ! Calculate displacements for both rake angles
            slip(1) = -slip_mag*cos(rak1*d2r)
            slip(2) = -slip_mag*sin(rak1*d2r)
            slip(3) = 0.0d0
            call tri_disloc_disp(disp1,sta_new,tri_new,poisson,slip)
            disp1(3) = -disp1(3) ! Returns displacement with positive z down, flip it
            slip(1) = -slip_mag*cos(rak2*d2r)
            slip(2) = -slip_mag*sin(rak2*d2r)
            slip(3) = 0.0d0
            call tri_disloc_disp(disp2,sta_new,tri_new,poisson,slip)
            disp2(3) = -disp2(3) ! Returns displacement with positive z down, flip it

        else
            call usage('calc_disp_gfs: no Greens function model named "'//trim(gf_model)//'"'//&
                       ' (usage:gf)')
        endif
        if (debug) then
            write(stderr,*) 'disp1:',disp1
            write(stderr,*) 'disp2:',disp2
        endif

        ! Load displacements into Green's function arrays
        if (displacement%file.ne.'none'.and.i.le.ndsp) then
            gf_disp%array(i       ,j     ) = disp1(1)
            gf_disp%array(i+1*ndsp,j     ) = disp1(2)
            gf_disp%array(i+2*ndsp,j     ) = disp1(3)
            gf_disp%array(i       ,j+nflt) = disp2(1)
            gf_disp%array(i+1*ndsp,j+nflt) = disp2(2)
            gf_disp%array(i+2*ndsp,j+nflt) = disp2(3)
        endif

        if (los%file.ne.'none'.and.i.gt.ndsp) then
            disp1(1) = disp1(1)*cos(los%array(i-ndsp,6)*d2r)*sin(los%array(i-ndsp,5)*d2r) + &
                       disp1(2)*cos(los%array(i-ndsp,6)*d2r)*cos(los%array(i-ndsp,5)*d2r) - &
                       disp1(3)*sin(los%array(i-ndsp,6)*d2r)
            disp2(1) = disp2(1)*cos(los%array(i-ndsp,6)*d2r)*sin(los%array(i-ndsp,5)*d2r) + &
                       disp2(2)*cos(los%array(i-ndsp,6)*d2r)*cos(los%array(i-ndsp,5)*d2r) - &
                       disp2(3)*sin(los%array(i-ndsp,6)*d2r)
            gf_los%array(i-ndsp,j     ) = disp1(1)
            gf_los%array(i-ndsp,j+nflt) = disp2(1)
        endif
    enddo

    if (verbosity.ge.1) then
        call progress_indicator(i,ndsp+nlos,'calc_disp_gfs',ierr)
        if (ierr.gt.0) then
            call usage('calc_disp_gfs: error in progress indicator subroutine (usage:none)')
        endif
    endif
enddo

if (verbosity.ge.1) then
    write(stdout,*) 'calc_disp_gfs: finished'
endif
if (verbosity.ge.3) then
    write(stdout,*) 'Three-component displacement GFs:'
    write(stdout,'(2A4,6A14)') 'sta','flt','dgf_ss_x','dgf_ss_y','dgf_ss_z',&
                                           'dgf_ds_x','dgf_ds_y','dgf_ds_z'
    do i = 1,ndsp
        do j = 1,nflt
            write(stdout,'(2I4,6F14.6)') i,j, &
                                         gf_disp%array(i       ,j), &
                                         gf_disp%array(i+1*ndsp,j), &
                                         gf_disp%array(i+2*ndsp,j), &
                                         gf_disp%array(i       ,j+nflt), &
                                         gf_disp%array(i+1*ndsp,j+nflt), &
                                         gf_disp%array(i+2*ndsp,j+nflt)
        enddo
    enddo
    write(stdout,*) 'Line-of-sight displacement GFs:'
    write(stdout,'(2A4,2A14)') 'sta','flt','dgf_ss_los','dgf_ds_los'
    do i = 1,nlos
        do j = 1,nflt
            write(stdout,'(2I4,2F14.6)') i,j,gf_los%array(i,j),gf_los%array(i,j+nflt)
        enddo
    enddo
endif
if (verbosity.ge.1) then
    write(stdout,*)
endif

return
end subroutine calc_disp_gfs




!--------------------------------------------------------------------------------------------------!
!----------------------------------- SHEAR TRACTION GFS -------------------------------------------!
!--------------------------------------------------------------------------------------------------!


subroutine calc_stress_gfs()
!----
! Calculate the Green's functions for shear tractions produced by each fault on every other fault
! (gf_stress), given the specified model (gf_model).
!----

use io, only: stdout, verbosity, progress_indicator
use trig, only: d2r, r2d
use algebra, only: rotate_matrix_angle_axis
use earth, only: radius_earth_m
use elast, only: strain2stress, stress2traction, traction_components
use geom, only: lola2distaz, strdip2normal
use okada92, only: o92_pt_strain, o92_rect_strain
use tri_disloc, only: tri_center, tri_geometry, tri_geo2cart, tri_disloc_strain

use fltinv, only: inversion_mode, &
                  fault, &
                  rake_constraint, &
                  prestress, &
                  gf_model, &
                  gf_stress, &
                  poisson, &
                  shearmod, &
                  lame, &
                  coord_type

implicit none

! Local variables
integer :: i, j, ierr, iTri
double precision, parameter :: eps = 1.0d-6
double precision :: evlo, evla, evdp, str, dip, wid, len, area, tri(3,4), tri_new(3,4), center(3)
double precision :: sta(3), sta_new(3), dist, az, dx, dy, pt1(3), pt2(3), pt3(3)
double precision :: slip_mag, rak1, rak2, slip(3), mom(4)
double precision :: strain1(3,3), strain2(3,3), stress1(3,3), stress2(3,3)
double precision :: nvec(3), tvec1(3), tvec2(3), tnor, tstr1, tupd1, tstr2, tupd2


if (verbosity.ge.1) then
    write(stdout,*) 'calc_stress_gfs: starting'
endif

if (prestress%file.eq.'none'.and.inversion_mode.ne.'anneal-psc') then
    if (verbosity.ge.1) then
        write(stdout,*) 'calc_stress_gfs: shear tractions are not used in this inversion'
        write(stdout,*)
    endif
    return
endif

if (gf_model.eq.'precomputed') then
    if (verbosity.ge.1) then
        write(stdout,*) 'calc_stress_gfs: stress Greens functions already provided'
        write(stdout,*)
    endif
    return
endif


! Set unit slip magnitude for the Green's functions calculation
slip_mag = 1.0d0


! Calculate stress Green's function for each fault-fault pair
do i = 1,fault%nrows

    ! Station (target fault) coordinates
    if (gf_model.eq.'okada_rect') then
        sta(1) = fault%array(i,1)
        sta(2) = fault%array(i,2)
        sta(3) = fault%array(i,3)
    elseif (gf_model.eq.'okada_pt') then
        call usage('calc_stress_gfs: cannot compute slip from shear tractions using okada_pt '//&
                   'GFs because there is a stress singularity at point source (usage:gf)')
    elseif (gf_model.eq.'triangle') then
        call tri_center(sta,fault%array(i,1:3),fault%array(i,4:6),fault%array(i,7:9))
    endif
    sta_new(3) = sta(3)

    ! Normal to plane at target fault
    if (gf_model.eq.'okada_rect') then
        call strdip2normal(fault%array(i,4),fault%array(i,5),nvec)
    elseif (gf_model.eq.'triangle') then
        if (coord_type.eq.'geographic') then
            call tri_geo2cart(pt1,pt2,pt3,fault%array(i,1:3),fault%array(i,4:6), &
                              fault%array(i,7:9),'m')
            call tri_geometry(nvec,tvec1,tvec2,pt1,pt2,pt3)
        elseif (coord_type.eq.'cartesian') then
            call tri_geometry(nvec,tvec1,tvec2,fault%array(i,1:3),fault%array(i,4:6), &
                              fault%array(i,7:9))
        endif
    endif

    do j = 1,fault%nrows

        ! Rake angles for computing Green's functions
        ! Default is strike-slip and dip-slip sources
        rak1 = 0.0d0
        rak2 = 90.0d0
        if (rake_constraint%file.ne.'none') then
            if (inversion_mode.eq.'lsqr') then
                ! lsqr: constraints change the rake for computing Green's functions
                rak1 = rake_constraint%array(j,1)
                if (rake_constraint%ncols.eq.2) then
                    rak2 = rake_constraint%array(j,2)
                endif
            elseif (inversion_mode.eq.'anneal-psc') then
                ! anneal-psc: constraint changes the rake for computing Green's functions
                rak1 = rake_constraint%array(j,1)
            endif
        endif

        ! The details of calculating the Green's functions differ depending on the model, but
        ! the general steps are the same:
        !     1. Calculate the location of the station relative to the center of the fault
        !     2. Calculate the shear tractions for unit slip in both rake directions

        if (gf_model.eq.'okada_rect') then
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str = fault%array(j,4)
            dip = fault%array(j,5)
            wid = fault%array(j,6)
            len = fault%array(j,7)

            ! Calculate location of station relative to fault center
            if (coord_type.eq.'geographic') then
                call lola2distaz(evlo,evla,sta(1),sta(2),dist,az,'radians','degrees',ierr)
                if (ierr.ne.0) then
                    call usage('calc_stress_gfs: error calculating dist and az (usage:none)')
                endif
                dist = dist*radius_earth_m
            elseif (coord_type.eq.'cartesian') then
                dx = sta(1) - evlo
                dy = sta(2) - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)*r2d
            endif

            ! Rotate coordinate axes so that x=strike, y=horizontal updip
            sta_new(1) = dist*( dcos((az-str)*d2r))
            sta_new(2) = dist*(-dsin((az-str)*d2r))

            ! Avoid singular solutions with measurement point lying on fault
            if (dabs(sta_new(1)).lt.eps) then
                sta_new(1) = sta_new(1) + eps
            endif
            if (dabs(sta_new(2)).lt.eps) then
                sta_new(2) = sta_new(2) + eps
            endif
            if (dabs(evdp-sta_new(3)).lt.eps) then
                sta_new(3) = evdp + eps
            endif

            ! Calculate strains for both rake angles
            slip(1) = slip_mag*cos(rak1*d2r)
            slip(2) = slip_mag*sin(rak1*d2r)
            slip(3) = 0.0d0
            call o92_rect_strain(strain1,sta_new,evdp,dip,slip,wid,len,lame,shearmod)
            slip(1) = slip_mag*cos(rak2*d2r)
            slip(2) = slip_mag*sin(rak2*d2r)
            slip(3) = 0.0d0
            call o92_rect_strain(strain2,sta_new,evdp,dip,slip,wid,len,lame,shearmod)

            ! Rotate strains back to ENZ
            call rotate_matrix_angle_axis(strain1,90.0d0-str,'z',strain1,ierr)
            call rotate_matrix_angle_axis(strain2,90.0d0-str,'z',strain2,ierr)

        elseif (gf_model.eq.'okada_pt') then
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str = fault%array(j,4)
            dip = fault%array(j,5)
            area = fault%array(j,6)

            ! Calculate location of station relative to fault center
            if (coord_type.eq.'geographic') then
                call lola2distaz(evlo,evla,sta(1),sta(2),dist,az,'radians','degrees',ierr)
                if (ierr.ne.0) then
                    call usage('calc_stress_gfs: error calculating dist and az (usage:none)')
                endif
                dist = dist*radius_earth_m
            elseif (coord_type.eq.'cartesian') then
                dx = sta(1) - evlo
                dy = sta(2) - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)*r2d
            endif

            ! Rotate coordinate axes so that x=strike, y=horizontal updip
            sta_new(1) = dist*( dcos((az-str)*d2r))
            sta_new(2) = dist*(-dsin((az-str)*d2r))

            ! Avoid singular solutions with measurement point lying on fault
            if (dabs(sta_new(1)).lt.eps) then
                sta_new(1) = sta_new(1) + eps
            endif
            if (dabs(sta_new(2)).lt.eps) then
                sta_new(2) = sta_new(2) + eps
            endif
            if (dabs(evdp-sta_new(3)).lt.eps) then
                sta_new(3) = evdp + eps
            endif

            ! Calculate strains for both rake angles
            mom(1) = slip_mag*cos(rak1*d2r)*area*shearmod
            mom(2) = slip_mag*sin(rak1*d2r)*area*shearmod
            mom(3:4) = 0.0d0
            call o92_pt_strain(strain1,sta_new,evdp,dip,mom,lame,shearmod)
            mom(1) = slip_mag*cos(rak2*d2r)*area*shearmod
            mom(2) = slip_mag*sin(rak2*d2r)*area*shearmod
            mom(3:4) = 0.0d0
            call o92_pt_strain(strain2,sta_new,evdp,dip,mom,lame,shearmod)

            ! Rotate strains back to ENZ
            call rotate_matrix_angle_axis(strain1,90.0d0-str,'z',strain1,ierr)
            call rotate_matrix_angle_axis(strain2,90.0d0-str,'z',strain2,ierr)

        elseif (gf_model.eq.'triangle') then
            tri(1,1) = fault%array(j,1)
            tri(2,1) = fault%array(j,2)
            tri(3,1) = fault%array(j,3)
            tri(1,2) = fault%array(j,4)
            tri(2,2) = fault%array(j,5)
            tri(3,2) = fault%array(j,6)
            tri(1,3) = fault%array(j,7)
            tri(2,3) = fault%array(j,8)
            tri(3,3) = fault%array(j,9)
            tri(:,4) = 0.0d0

            ! Calculate location of station relative to fault center
            if (coord_type.eq.'geographic') then
                call tri_center(center,tri(:,1),tri(:,2),tri(:,3))
                call lola2distaz(center(1),center(2),sta(1),sta(2),dist,az,'radians','radians',ierr)
                if (ierr.ne.0) then
                    call usage('calc_stress_gfs: error calculating dist and az (usage:none)')
                endif
                dist = dist*radius_earth_m
                sta_new(1) = dist*dsin(az)
                sta_new(2) = dist*dcos(az)
                do iTri = 1,3
                    call lola2distaz(center(1),center(2),tri(1,iTri),tri(2,iTri),dist,az,'radians','radians',ierr)
                    if (ierr.ne.0) then
                        call usage('calc_stress_gfs: error calculating dist and az (usage:none)')
                    endif
                    dist = dist*radius_earth_m
                    tri_new(1,iTri) = dist*dsin(az)
                    tri_new(2,iTri) = dist*dcos(az)
                    tri_new(3,iTri) = tri(3,iTri)
                enddo
            elseif (coord_type.eq.'cartesian') then
                tri_new = tri
                sta_new = sta
            endif

            ! Calculate strains for both rake angles
            slip(1) = -slip_mag*cos(rak1*d2r)
            slip(2) = -slip_mag*sin(rak1*d2r)
            slip(3) = 0.0d0
            call tri_disloc_strain(strain1,sta_new,tri_new,poisson,slip)
            strain1(1,3) = -strain1(1,3)
            strain1(3,1) = -strain1(1,3)
            strain1(2,3) = -strain1(2,3)
            strain1(3,2) = -strain1(2,3)

            slip(1) = -slip_mag*cos(rak2*d2r)
            slip(2) = -slip_mag*sin(rak2*d2r)
            slip(3) = 0.0d0
            call tri_disloc_strain(strain2,sta_new,tri_new,poisson,slip)
            strain2(1,3) = -strain2(1,3)
            strain2(3,1) = -strain2(1,3)
            strain2(2,3) = -strain2(2,3)
            strain2(3,2) = -strain2(2,3)

        else
            call usage('calc_stress_gfs: no Greens function model named "'//trim(gf_model)//'"'//&
                       '(usage:gf)')
        endif

        ! Calculate shear tractions on fault from strain tensor

        ! Resolve traction vector for strike-slip source
        call strain2stress(strain1,lame,shearmod,stress1)
        call stress2traction(stress1,nvec,tvec1)
        call traction_components(tvec1,nvec,tnor,tstr1,tupd1)

        ! Resolve traction vector for dip-slip source
        call strain2stress(strain2,lame,shearmod,stress2)
        call stress2traction(stress2,nvec,tvec2)
        call traction_components(tvec2,nvec,tnor,tstr2,tupd2)

        ! Load traction components into Green's function arrays
        ! We use the opposite sign because we are interested in the fault slip that reduces the
        ! traction resolved on the structure to zero, i.e.:
        !     trac_from_slip + applied_traction = 0
        gf_stress%array(i            ,j            ) = -tstr1
        gf_stress%array(i+fault%nrows,j            ) = -tupd1
        gf_stress%array(i            ,j+fault%nrows) = -tstr2
        gf_stress%array(i+fault%nrows,j+fault%nrows) = -tupd2
    enddo

    if (verbosity.ge.1) then
        call progress_indicator(i,fault%nrows,'calc_stress_gfs',ierr)
        if (ierr.gt.0) then
            call usage('calc_stress_gfs: error in progress indicator subroutine (usage:none)')
        endif
    endif
enddo

if (verbosity.ge.1) then
    write(stdout,*) 'calc_stress_gfs: finished'
endif
if (verbosity.ge.3) then
    write(stdout,*) 'Fault-fault shear traction GFs'
    write(stdout,'(2A4,4A14)') 'sta','flt','sgf_ss_str','sgf_ss_dip','sgf_ds_str','sgf_ds_dip'
    do i = 1,fault%nrows
        do j = 1,fault%nrows
            write(stdout,'(2I4,1P4E14.6)') i,j, &
                                           gf_stress%array(i            ,j            ), &
                                           gf_stress%array(i+fault%nrows,j            ), &
                                           gf_stress%array(i            ,j+fault%nrows), &
                                           gf_stress%array(i+fault%nrows,j+fault%nrows)
        enddo
    enddo
endif
if (verbosity.ge.1) then
    write(stdout,*)
endif

return
end subroutine calc_stress_gfs


!--------------------------------------------------------------------------------------------------!
!----------------------------------- RIGID BODY ROTATION GFS --------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine calc_euler_gfs()
!----
! Calculate the Green's functions for rigid body rotations at selected observation points.
!----

use io, only: verbosity, stdout
use trig, only: d2r
use earth, only: deg2km

use fltinv, only: euler_file, &
                  npoles, &
                  rigid_pt_array_disp, &
                  rigid_pt_array_los, &
                  gf_euler, &
                  coord_type, &
                  displacement, &
                  los

implicit none

! Local variables
integer :: i, j, iblock, ndsp
double precision :: angular_velocity, lon, lat, dep, r(3)
double precision, parameter :: spy = 60.0d0*60.0d0*24.0d0*365.25d0

if (verbosity.ge.1) then
    write(stdout,*) 'calc_euler_gfs: starting'
endif

if (euler_file.eq.'none') then
    if (verbosity.ge.1) then
        write(stdout,*) 'calc_euler_gfs: rigid body rotations are not determined in this inversion'
        write(stdout,*)
    endif
    return
endif

if (coord_type.eq.'cartesian') then
    write(stdout,*) 'calc_euler_gfs: rigid body rotations are not used in Cartesian coordinate mode'
    write(stdout,*)
    return
endif


! Set unit angular velocity for the Green's functions calculation
angular_velocity = 1.0d0


if (displacement%file.ne.'none') then
    ndsp = displacement%nrows
    do i = 1,ndsp

        iblock = rigid_pt_array_disp(i)

        if (iblock.ne.0) then
            lon = displacement%array(i,1)
            lat = displacement%array(i,2)
            dep = displacement%array(i,3)
            if (abs(dep).gt.1.0d0) then
                call usage('calc_euler_gfs: rigid plate motions only work for surface '//&
                           'observations (usage:input)')
            endif

            r(1) = dcos(lat*d2r)*dcos(lon*d2r)
            r(2) = dcos(lat*d2r)*dsin(lon*d2r)
            r(3) = dsin(lat*d2r)

            ! Angular velocity units of deg/Ma translates to GF units of deg/Ma
            gf_euler(i     ,iblock         ) = -r(3)*cos(lon*d2r)
            gf_euler(i+ndsp,iblock         ) =  r(3)*sin(lon*d2r)*sin(lat*d2r) + r(2)*cos(lat*d2r)
            gf_euler(i     ,iblock+  npoles) = -r(3)*sin(lon*d2r)
            gf_euler(i+ndsp,iblock+  npoles) = -r(3)*cos(lon*d2r)*sin(lat*d2r) - r(1)*cos(lat*d2r)
            gf_euler(i     ,iblock+2*npoles) =  r(2)*sin(lon*d2r)              + r(1)*cos(lon*d2r)
            gf_euler(i+ndsp,iblock+2*npoles) =  r(2)*cos(lon*d2r)*sin(lat*d2r) - &
                                                                     r(1)*sin(lon*d2r)*sin(lat*d2r)

            ! Convert Green's functions to units of m/s (i.e., one deg/Ma produces these GFs in m/s)
            gf_euler(i     ,iblock         ) = deg2km*gf_euler(i     ,iblock         )/1.0d3/spy
            gf_euler(i+ndsp,iblock         ) = deg2km*gf_euler(i+ndsp,iblock         )/1.0d3/spy
            gf_euler(i     ,iblock+  npoles) = deg2km*gf_euler(i     ,iblock+  npoles)/1.0d3/spy
            gf_euler(i+ndsp,iblock+  npoles) = deg2km*gf_euler(i+ndsp,iblock+  npoles)/1.0d3/spy
            gf_euler(i     ,iblock+2*npoles) = deg2km*gf_euler(i     ,iblock+2*npoles)/1.0d3/spy
            gf_euler(i+ndsp,iblock+2*npoles) = deg2km*gf_euler(i+ndsp,iblock+2*npoles)/1.0d3/spy
        endif
    enddo
endif


if (verbosity.ge.1) then
    write(stdout,*) 'calc_euler_gfs: finished'
endif
if (verbosity.ge.3) then
    write(stdout,*) 'Rigid body rotation GFs'
    write(stdout,'(2A4,6A14)') 'sta','pol','ang_x_ve','ang_x_vn','ang_y_ve','ang_y_vn','ang_z_ve', &
                               'ang_z_vn'
    do i = 1,ndsp
        iblock = rigid_pt_array_disp(i)
        do j = 1,npoles
            write(stdout,'(2I4,1P6E14.6)') i,j, &
                                           gf_euler(i     ,iblock         ), &
                                           gf_euler(i+ndsp,iblock         ), &
                                           gf_euler(i     ,iblock+  npoles), &
                                           gf_euler(i+ndsp,iblock+  npoles), &
                                           gf_euler(i     ,iblock+2*npoles), &
                                           gf_euler(i+ndsp,iblock+2*npoles)
        enddo
    enddo
endif
if (verbosity.ge.1) then
    write(stdout,*)
endif


return
end subroutine calc_euler_gfs
