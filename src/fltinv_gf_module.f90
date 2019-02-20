module gf_module

!--------------------------------------------------------------------------------------------------!

contains

!--------------------------------------------------------------------------------------------------!
! RECTANGULAR DISLOCATIONS IN AN ELASTIC HALF-SPACE
!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_disp_okada_rect()
    use io_module, only: stderr, verbosity
    use variable_module, only: inversion_mode, coord_type, &
                               displacement, fault, gf_disp, halfspace, rake_constraint
    implicit none
    ! Local variables
    integer :: i, j
    double precision :: slip, stlo, stla, stdp, evlo, evla, evdp, str, dip, rak, wid, len
    double precision :: dx, dy, dist, az, x, y, ux, uy, uz, uxp, uyp, theta, uhor, vp, vs, dens
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r = pi/180.0d0

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_disp_okada_rect says: starting'
    endif

    ! Unit slip for displacement GF computation
    slip = 1.0d0

    ! Half-space array contains vp, vs, density
    vp = halfspace%array(1,1)
    vs = halfspace%array(1,2)
    dens = halfspace%array(1,3)

    ! Displacement Green's function for each fault-station pair
    do i = 1,displacement%nrecords
        stlo = displacement%array(i,1)
        stla = displacement%array(i,2)
        stdp = displacement%array(i,3)
        if (verbosity.ge.1) then
            write(0,'(A1,A,I3,A)',advance='no') achar(13), &
                'calc_gf_disp_okada_rect progress: ',100*i/displacement%nrecords,'%'
            if (i.eq.displacement%nrecords) then
                write(0,*)
            endif
        endif

        do j = 1,fault%nrecords
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str  = fault%array(j,4)
            dip  = fault%array(j,5)
            wid  = fault%array(j,6)
            len  = fault%array(j,7)

            ! Origin is at fault coordinate, x-axis points horizontal up-dip
            if (coord_type.eq.'geographic') then
                call ddistaz(dist,az,evlo,evla,stlo,stla)
                dist = dist*6.371d6
            else
                dx = stlo - evlo
                dy = stla - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)
            endif
            x = dist*( dcos(az-d2r*str))
            y = dist*(-dsin(az-d2r*str))

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 0.0d0
            endif

            ! Strike-slip Green's function
            call o92rect(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,wid,len,slip,vp,vs,dens)
            ! Rotate back to original x-y coordinates
            theta = datan2(uyp,uxp)
            uhor = dsqrt(uxp*uxp+uyp*uyp)
            theta = d2r*str - theta
            ux = uhor*dsin(theta)
            uy = uhor*dcos(theta)
            gf_disp%array(i,                        j) = ux
            gf_disp%array(i+1*displacement%nrecords,j) = uy
            gf_disp%array(i+2*displacement%nrecords,j) = uz

            ! Check for rake constraints; if none, calculate dip-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
                    .and.rake_constraint%nfields.eq.2) then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,2)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,2)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 90.0d0
            endif

            ! Dip-slip (or second rake) Green's function
            call o92rect(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,wid,len,slip,vp,vs,dens)
            ! Rotate back to original x-y coordinates
            theta = datan2(uyp,uxp)
            uhor = dsqrt(uxp*uxp+uyp*uyp)
            theta = d2r*str - theta
            ux = uhor*dsin(theta)
            uy = uhor*dcos(theta)
            gf_disp%array(i,                        j+fault%nrecords) = ux
            gf_disp%array(i+1*displacement%nrecords,j+fault%nrecords) = uy
            gf_disp%array(i+2*displacement%nrecords,j+fault%nrecords) = uz
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_disp_okada_rect says: finished'
    endif
    if (verbosity.ge.3) then
        write(stderr,'(2A4,6A14)') 'sta','flt','dgf_ss_x','dgf_ss_y','dgf_ss_z',&
                                               'dgf_ds_x','dgf_ds_y','dgf_ds_z'
        do i = 1,displacement%nrecords
            do j = 1,fault%nrecords
                write(stderr,'(2I4,6F14.6)') i,j, &
                    gf_disp%array(i                        ,j), &
                    gf_disp%array(i+1*displacement%nrecords,j), &
                    gf_disp%array(i+2*displacement%nrecords,j), &
                    gf_disp%array(i                        ,j+fault%nrecords), &
                    gf_disp%array(i+1*displacement%nrecords,j+fault%nrecords), &
                    gf_disp%array(i+2*displacement%nrecords,j+fault%nrecords)
            enddo
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine calc_gf_disp_okada_rect

!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_stress_okada_rect()
    use io_module, only: stderr, verbosity
    use variable_module, only: inversion_mode, coord_type, fault, rake_constraint, halfspace, &
                               gf_stress, sts_dist
    use elast
    implicit none
    ! Local variables
    double precision :: slip, stlo, stla, stdp, sta_str, sta_dip
    double precision :: evlo, evla, evdp, str, dip, rak, wid, len
    double precision :: dist, az, dx, dy, x, y, vp, vs, dens
    double precision :: unit_normal(3), unit_strike(3), unit_updip(3)
    double precision :: strain(3,3), stress(3,3), traction(3), traction_components(3)
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r=pi/180.0d0
    integer :: i, j, prog

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_stress_okada_rect says: starting'
    endif
    prog = 0

    ! Unit slip for GF computation
    slip = 1.0d0

    ! Half-space array contains vp, vs, density
    vp = halfspace%array(1,1)
    vs = halfspace%array(1,2)
    dens = halfspace%array(1,3)

    ! Shear stress Green's functions for each fault-fault pair
    do i = 1,fault%nrecords
        if (verbosity.ge.1) then
            write(0,'(A1,A,I3,A)',advance='no') achar(13), &
                'calc_gf_stress_okada_rect progress: ',100*i/fault%nrecords,'%'
            if (i.eq.fault%nrecords) then
                write(0,*)
            endif
        endif
        stlo = fault%array(i,1)
        stla = fault%array(i,2)
        stdp = fault%array(i,3)
        sta_str = fault%array(i,4)
        sta_dip = fault%array(i,5)
        call calc_plane_unit_vectors(sta_str,sta_dip,unit_normal,unit_strike,unit_updip)

        do j = 1,fault%nrecords
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str  = fault%array(j,4)
            dip  = fault%array(j,5)
            wid  = fault%array(j,6)
            len  = fault%array(j,7)

            ! Rotate coordinates to x=along-strike, y=horizontal up-dip
            if (coord_type.eq.'geographic') then
                call ddistaz(dist,az,evlo,evla,stlo,stla)
                dist = dist*6.371d6
            else
                dx = stlo - evlo
                dy = stla - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)
            endif
             x = dist*( dcos(az-d2r*str))
             y = dist*(-dsin(az-d2r*str))

            ! Avoid singular solutions with measurement point lying on fault
            if (dabs(x).lt.1.0d0) then
                x = x + 1.0d0
            endif
            if (dabs(y).lt.1.0d0) then
                y = y + 1.0d0
            endif
            if (dabs(evdp-stdp).lt.1.0d0) then
                stdp = evdp + 1.0d0
            endif

            ! Set resolved tractions to zero if distance is larger than threshold
            dist = dsqrt(x*x+ &
                         y*y+ &
                        (stdp-evdp)*(stdp-evdp))
            if (dist.gt.sts_dist) then
                gf_stress%array(i,j) = 0.0d0
                gf_stress%array(i,j+fault%nrecords) = 0.0d0
                gf_stress%array(i+fault%nrecords,j) = 0.0d0
                gf_stress%array(i+fault%nrecords,j+fault%nrecords) = 0.0d0
                cycle
            endif

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 0.0d0
            endif

            ! Shear stress Green's function produced by strike-slip or fixed-rake source
            call o92rectstn(strain,x,y,stdp,evdp,dip,rak,wid,len,slip,vp,vs,dens)
            ! Rotate back to original coordinates
            call rotate_strain(strain,str)
            call calc_strain_to_stress(strain,vp,vs,dens,stress)
            call calc_traction(stress,unit_normal,traction)
            call calc_traction_components(traction,unit_normal,unit_strike,unit_updip, &
                                          traction_components)
            ! Load stress Green's functions array
            gf_stress%array(i               ,j               ) = -traction_components(2)
            gf_stress%array(i+fault%nrecords,j               ) = -traction_components(3)

            ! Shear stress Green's function produced by dip-slip source
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
                    .and.rake_constraint%nfields.eq.2) then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,2)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,2)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 90.0d0
            endif

            call o92rectstn(strain,x,y,stdp,evdp,dip,rak,wid,len,slip,vp,vs,dens)
            ! Rotate back to original coordinates
            call rotate_strain(strain,str)
            call calc_strain_to_stress(strain,vp,vs,dens,stress)
            call calc_traction(stress,unit_normal,traction)
            call calc_traction_components(traction,unit_normal,unit_strike,unit_updip, &
                                          traction_components)
            ! Load stress Green's functions array
            gf_stress%array(i               ,j+fault%nrecords) = -traction_components(2)
            gf_stress%array(i+fault%nrecords,j+fault%nrecords) = -traction_components(3)
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') "calc_gf_stress_okada_rect says: finished"
    endif
    if (verbosity.ge.3) then
        write(stderr,'(2A4,4A14)') 'sta','flt','sgf_ss_ss','sgf_ss_ds','sgf_ds_ss','sgf_ds_ds'
        do i = 1,fault%nrecords
            do j = 1,fault%nrecords
                write(stderr,'(2I4,1P4E14.6)') i,j, &
                    gf_stress%array(i               ,j), &
                    gf_stress%array(i+fault%nrecords,j), &
                    gf_stress%array(i               ,j+fault%nrecords), &
                    gf_stress%array(i+fault%nrecords,j+fault%nrecords)
            enddo
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine calc_gf_stress_okada_rect

!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_los_okada_rect()
    use io_module, only: stderr, verbosity
    use variable_module, only: inversion_mode, coord_type, &
                               los, fault, gf_los, halfspace, rake_constraint
    implicit none
    ! Local variables
    integer :: i, j
    double precision :: slip, stlo, stla, stdp, evlo, evla, evdp, str, dip, rak, wid, len
    double precision :: dx, dy, dist, az, x, y, ux, uy, uz, uxp, uyp, theta, uhor, vp, vs, dens
    double precision :: lookaz, lookinc
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r = pi/180.0d0

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_los_okada_rect says: starting'
    endif

    ! Unit slip for displacement GF computation
    slip = 1.0d0

    ! Half-space array contains vp, vs, density
    vp = halfspace%array(1,1)
    vs = halfspace%array(1,2)
    dens = halfspace%array(1,3)

    ! LOS displacement Green's function for each fault-station pair
    do i = 1,los%nrecords
        stlo = los%array(i,1)
        stla = los%array(i,2)
        stdp = los%array(i,3)
        lookaz   = los%array(i,5)*d2r
        lookinc  = los%array(i,6)*d2r
        do j = 1,fault%nrecords
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str  = fault%array(j,4)
            dip  = fault%array(j,5)
            wid  = fault%array(j,6)
            len  = fault%array(j,7)

            ! Origin is at fault coordinate, x-axis points horizontal up-dip
            if (coord_type.eq.'geographic') then
                call ddistaz(dist,az,evlo,evla,stlo,stla)
                dist = dist*6.371d6
            else
                dx = stlo - evlo
                dy = stla - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)
            endif
             x = dist*( dcos(az-d2r*str))
             y = dist*(-dsin(az-d2r*str))

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 0.0d0
            endif

            ! Strike-slip Green's function
            call o92rect(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,wid,len,slip,vp,vs,dens)
            ! Rotate back to original x-y coordinates
            theta = datan2(uyp,uxp)
            uhor = dsqrt(uxp*uxp+uyp*uyp)
            theta = d2r*str - theta
            ux = uhor*dsin(theta)
            uy = uhor*dcos(theta)
            gf_los%array(i,j) = &
                         ux*sin(lookaz)*cos(lookinc) + uy*cos(lookaz)*cos(lookinc) - uz*sin(lookinc)

            ! Check for rake constraints; if none, calculate dip-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
                    .and.rake_constraint%nfields.eq.2) then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,2)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,2)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 90.0d0
            endif

            call o92rect(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,wid,len,slip,vp,vs,dens)
            ! Rotate back to original x-y coordinates
            theta = datan2(uyp,uxp)
            uhor = dsqrt(uxp*uxp+uyp*uyp)
            theta = d2r*str - theta
            ux = uhor*dsin(theta)
            uy = uhor*dcos(theta)
            gf_los%array(i,j+fault%nrecords) = &
                         ux*sin(lookaz)*cos(lookinc) + uy*cos(lookaz)*cos(lookinc) - uz*sin(lookinc)
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_los_okada_rect says: finished'
    endif
    if (verbosity.ge.3) then
        write(stderr,'(2A4,6A14)') 'sta','flt','losgf_ss','losgf_ds'
        do i = 1,los%nrecords
            do j = 1,fault%nrecords
                write(stderr,'(2I4,6F14.6)') i,j, &
                    gf_los%array(i,j), gf_los%array(i,j+fault%nrecords)
            enddo
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine calc_gf_los_okada_rect

!--------------------------------------------------------------------------------------------------!
! POINT DISLOCATIONS IN AN ELASTIC HALF-SPACE
!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_disp_okada_pt()
    use io_module, only: stderr, verbosity
    use variable_module, only: inversion_mode, coord_type, &
                               displacement, fault, gf_disp, halfspace, rake_constraint
    implicit none
    ! Local variables
    integer :: i, j
    double precision :: slip, stlo, stla, stdp, evlo, evla, evdp, str, dip, rak, area
    double precision :: dx, dy, dist, az, x, y, ux, uy, uz, uxp, uyp, theta, uhor, vp, vs, dens
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r = pi/180.0d0

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_disp_okada_pt says: starting'
    endif

    ! Unit slip for displacement GF computation
    slip = 1.0d0

    ! Half-space array contains vp, vs, density
    vp = halfspace%array(1,1)
    vs = halfspace%array(1,2)
    dens = halfspace%array(1,3)

    ! Displacement Green's function for each fault-station pair
    do i = 1,displacement%nrecords
        stlo = displacement%array(i,1)
        stla = displacement%array(i,2)
        stdp = displacement%array(i,3)
        do j = 1,fault%nrecords
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str  = fault%array(j,4)
            dip  = fault%array(j,5)
            area  = fault%array(j,6)

            ! Origin is at fault coordinate, x-axis points horizontal up-dip
            if (coord_type.eq.'geographic') then
                call ddistaz(dist,az,evlo,evla,stlo,stla)
                dist = dist*6.371d6
            else
                dx = stlo - evlo
                dy = stla - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)
            endif
             x = dist*( dcos(az-d2r*str))
             y = dist*(-dsin(az-d2r*str))

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call usage('!! calc_gf_disp_okada_pt: incorrect number of rake '// &
                                     'constraints')
                endif
            else
                rak = 0.0d0
            endif

            ! Strike-slip Green's function
            call o92pt(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,area,slip,vp,vs,dens)
            ! Rotate back to original x-y coordinates
            theta = datan2(uyp,uxp)
            uhor = dsqrt(uxp*uxp+uyp*uyp)
            theta = d2r*str - theta
            ux = uhor*dsin(theta)
            uy = uhor*dcos(theta)
            gf_disp%array(i,                        j) = ux
            gf_disp%array(i+1*displacement%nrecords,j) = uy
            gf_disp%array(i+2*displacement%nrecords,j) = uz

            ! Check for rake constraints; if none, calculate dip-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
                    .and.rake_constraint%nfields.eq.2) then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,2)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,2)
                else
                    call usage('!! calc_gf_disp_okada_pt: incorrect number of rake '// &
                                     'constraints')
                endif
            else
                rak = 90.0d0
            endif

            ! Dip-slip (or second rake) Green's function
            call o92pt(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,area,slip,vp,vs,dens)
            ! Rotate back to original x-y coordinates
            theta = datan2(uyp,uxp)
            uhor = dsqrt(uxp*uxp+uyp*uyp)
            theta = d2r*str - theta
            ux = uhor*dsin(theta)
            uy = uhor*dcos(theta)
            gf_disp%array(i,                        j+fault%nrecords) = ux
            gf_disp%array(i+1*displacement%nrecords,j+fault%nrecords) = uy
            gf_disp%array(i+2*displacement%nrecords,j+fault%nrecords) = uz
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_disp_okada_pt says: finished'
    endif
    if (verbosity.ge.3) then
        write(stderr,'(2A4,6A14)') 'sta','flt','dgf_ss_x','dgf_ss_y','dgf_ss_z',&
                                               'dgf_ds_x','dgf_ds_y','dgf_ds_z'
        do i = 1,displacement%nrecords
            do j = 1,fault%nrecords
                write(stderr,'(2I4,6F14.6)') i,j, &
                    gf_disp%array(i                        ,j), &
                    gf_disp%array(i+1*displacement%nrecords,j), &
                    gf_disp%array(i+2*displacement%nrecords,j), &
                    gf_disp%array(i                        ,j+fault%nrecords), &
                    gf_disp%array(i+1*displacement%nrecords,j+fault%nrecords), &
                    gf_disp%array(i+2*displacement%nrecords,j+fault%nrecords)
            enddo
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine calc_gf_disp_okada_pt

!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_stress_okada_pt()
    use io_module, only: stderr, verbosity
    use variable_module, only: inversion_mode, coord_type, fault, rake_constraint, halfspace, &
                               gf_stress
    use elast
    implicit none
    ! Local variables
    double precision :: slip, stlo, stla, stdp, sta_str, sta_dip
    double precision :: evlo, evla, evdp, str, dip, rak, area
    double precision :: dist, az, dx, dy, x, y, vp, vs, dens
    double precision :: unit_normal(3), unit_strike(3), unit_updip(3)
    double precision :: strain(3,3), stress(3,3), traction(3), traction_components(3)
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r=pi/180.0d0
    integer :: i, j, prog

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_stress_okada_pt says: starting'
    endif
    prog = 0

    ! Unit slip for GF computation
    slip = 1.0d0

    ! Half-space array contains vp, vs, density
    vp = halfspace%array(1,1)
    vs = halfspace%array(1,2)
    dens = halfspace%array(1,3)

    ! Shear stress Green's functions for each fault-fault pair
    do i = 1,fault%nrecords
        if (i.gt.prog) then
            if (verbosity.ge.2) then
                write(0,'(A,I3,A)') 'calc_gf_stress_okada_pt says:',100*prog/fault%nrecords,&
                               '% complete'
            endif
            prog = prog + fault%nrecords/100
        endif
        stlo = fault%array(i,1)
        stla = fault%array(i,2)
        stdp = fault%array(i,3)
        sta_str = fault%array(i,4)
        sta_dip = fault%array(i,5)
        call calc_plane_unit_vectors(sta_str,sta_dip,unit_normal,unit_strike,unit_updip)

        do j = 1,fault%nrecords
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str  = fault%array(j,4)
            dip  = fault%array(j,5)
            area  = fault%array(j,6)

            ! Rotate coordinates to x=along-strike, y=horizontal up-dip
            if (coord_type.eq.'geographic') then
                call ddistaz(dist,az,evlo,evla,stlo,stla)
                dist = dist*6.371d6
            else
                dx = stlo - evlo
                dy = stla - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)
            endif
             x = dist*( dcos(az-d2r*str))
             y = dist*(-dsin(az-d2r*str))

            ! Avoid singular solutions with measurement point lying on fault
            if (dabs(x).lt.1.0d0) then
                x = x + 1.0d0
            endif
            if (dabs(y).lt.1.0d0) then
                y = y + 1.0d0
            endif
            if (dabs(evdp-stdp).lt.1.0d0) then
                stdp = evdp + 1.0d0
            endif

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call usage('!! calc_gf_stress_okada_pt: incorrect number of rake '// &
                                     'constraints')
                endif
            else
                rak = 0.0d0
            endif

            ! Shear stress Green's function produced by strike-slip or fixed-rake source
            call o92ptstn(strain,x,y,stdp,evdp,dip,rak,area,slip,vp,vs,dens)
            ! Rotate back to original coordinates
            call rotate_strain(strain,str)
            call calc_strain_to_stress(strain,vp,vs,dens,stress)
            call calc_traction(stress,unit_normal,traction)
            call calc_traction_components(traction,unit_normal,unit_strike,unit_updip, &
                                          traction_components)
            ! Load stress Green's functions array
            gf_stress%array(i               ,j               ) = -traction_components(2)
            gf_stress%array(i+fault%nrecords,j               ) = -traction_components(3)

            ! Shear stress Green's function produced by dip-slip source
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
                    .and.rake_constraint%nfields.eq.2) then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,2)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,2)
                else
                    call usage('!! calc_gf_stress_okada_pt: incorrect number of rake '// &
                                     'constraints')
                endif
            else
                rak = 90.0d0
            endif

            call o92ptstn(strain,x,y,stdp,evdp,dip,rak,area,slip,vp,vs,dens)
            ! Rotate back to original coordinates
            call rotate_strain(strain,str)
            call calc_strain_to_stress(strain,vp,vs,dens,stress)
            call calc_traction(stress,unit_normal,traction)
            call calc_traction_components(traction,unit_normal,unit_strike,unit_updip, &
                                          traction_components)
            ! Load stress Green's functions array
            gf_stress%array(i               ,j+fault%nrecords) = -traction_components(2)
            gf_stress%array(i+fault%nrecords,j+fault%nrecords) = -traction_components(3)
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') "calc_gf_stress_okada_pt says: finished"
    endif
    if (verbosity.ge.3) then
        write(stderr,'(2A4,4A14)') 'sta','flt','sgf_ss_ss','sgf_ss_ds','sgf_ds_ss','sgf_ds_ds'
        do i = 1,fault%nrecords
            do j = 1,fault%nrecords
                write(stderr,'(2I4,1P4E14.6)') i,j, &
                    gf_stress%array(i               ,j), &
                    gf_stress%array(i+fault%nrecords,j), &
                    gf_stress%array(i               ,j+fault%nrecords), &
                    gf_stress%array(i+fault%nrecords,j+fault%nrecords)
            enddo
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine calc_gf_stress_okada_pt

!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_los_okada_pt()
    use io_module, only: stderr, verbosity
    use variable_module, only: inversion_mode, coord_type, &
                               los, fault, gf_los, halfspace, rake_constraint
    implicit none
    ! Local variables
    integer :: i, j
    double precision :: slip, stlo, stla, stdp, evlo, evla, evdp, str, dip, rak, area
    double precision :: dx, dy, dist, az, x, y, ux, uy, uz, uxp, uyp, theta, uhor, vp, vs, dens
    double precision :: lookaz, lookinc
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r = pi/180.0d0

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_los_okada_pt says: starting'
    endif

    ! Unit slip for displacement GF computation
    slip = 1.0d0

    ! Half-space array contains vp, vs, density
    vp = halfspace%array(1,1)
    vs = halfspace%array(1,2)
    dens = halfspace%array(1,3)

    ! LOS displacement Green's function for each fault-station pair
    do i = 1,los%nrecords
        stlo = los%array(i,1)
        stla = los%array(i,2)
        stdp = los%array(i,3)
        lookaz   = los%array(i,5)*d2r
        lookinc  = los%array(i,6)*d2r
        do j = 1,fault%nrecords
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str  = fault%array(j,4)
            dip  = fault%array(j,5)
            area  = fault%array(j,6)

            ! Origin is at fault coordinate, x-axis points horizontal up-dip
            if (coord_type.eq.'geographic') then
                call ddistaz(dist,az,evlo,evla,stlo,stla)
                dist = dist*6.371d6
            else
                dx = stlo - evlo
                dy = stla - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)
            endif
             x = dist*( dcos(az-d2r*str))
             y = dist*(-dsin(az-d2r*str))

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call usage('!! calc_gf_los_okada_pt: incorrect number of rake '// &
                                     'constraints')
                endif
            else
                rak = 0.0d0
            endif

            ! Strike-slip Green's function
            call o92pt(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,area,slip,vp,vs,dens)
            ! Rotate back to original x-y coordinates
            theta = datan2(uyp,uxp)
            uhor = dsqrt(uxp*uxp+uyp*uyp)
            theta = d2r*str - theta
            ux = uhor*dsin(theta)
            uy = uhor*dcos(theta)
            gf_los%array(i,j) = &
                         ux*sin(lookaz)*cos(lookinc) + uy*cos(lookaz)*cos(lookinc) - uz*sin(lookinc)

            ! Check for rake constraints; if none, calculate dip-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
                    .and.rake_constraint%nfields.eq.2) then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,2)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,2)
                else
                    call usage('!! calc_gf_los_okada_pt: incorrect number of rake '// &
                                     'constraints')
                endif
            else
                rak = 90.0d0
            endif

            call o92pt(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,area,slip,vp,vs,dens)
            ! Rotate back to original x-y coordinates
            theta = datan2(uyp,uxp)
            uhor = dsqrt(uxp*uxp+uyp*uyp)
            theta = d2r*str - theta
            ux = uhor*dsin(theta)
            uy = uhor*dcos(theta)
            gf_los%array(i,j+fault%nrecords) = &
                         ux*sin(lookaz)*cos(lookinc) + uy*cos(lookaz)*cos(lookinc) - uz*sin(lookinc)
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_los_okada_pt says: finished'
    endif
    if (verbosity.ge.3) then
        write(stderr,'(2A4,6A14)') 'sta','flt','losgf_ss','losgf_ds'
        do i = 1,los%nrecords
            do j = 1,fault%nrecords
                write(stderr,'(2I4,6F14.6)') i,j, &
                    gf_los%array(i,j), gf_los%array(i,j+fault%nrecords)
            enddo
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine calc_gf_los_okada_pt

!--------------------------------------------------------------------------------------------------!
! TRIANGLE DISLOCATIONS IN AN ELASTIC HALF-SPACE
!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_disp_tri()
    use io_module, only: stderr, verbosity
    use variable_module, only: inversion_mode, coord_type, &
                               displacement, fault, gf_disp, halfspace, rake_constraint
    use tri_disloc_module, only: tri_disloc_disp, tri_center
    implicit none
    ! Local variables
    integer :: i, j, iTri
    double precision :: slip(3), slip_magnitude, rak, sta_coord(3), tri_coord(3,4)
    double precision :: sta_coord_new(3), tri_coord_new(3,4), disp(3), center(3), dist, az
    double precision :: lambda, shear_modulus, poisson
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r = pi/180.0d0

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_disp_tri says: starting'
    endif

    ! Initialize rake angle
    rak = 0.0d0

    ! Unit slip for displacement GF computation
    slip_magnitude = 1.0d0

    ! Half-space array contains vp, vs, density
    shear_modulus = halfspace%array(1,2)*halfspace%array(1,2)*halfspace%array(1,3)
    lambda = halfspace%array(1,1)*halfspace%array(1,1)*halfspace%array(1,3) - 2.0d0*shear_modulus
    poisson = lambda/2.0d0/(lambda+shear_modulus)

    ! Displacement Green's function for each fault-station pair
    do i = 1,displacement%nrecords
        sta_coord(1) = displacement%array(i,1)
        sta_coord(2) = displacement%array(i,2)
        sta_coord(3) = displacement%array(i,3)
        if (verbosity.ge.1) then
            write(0,'(A1,A,I3,A)',advance='no') achar(13), &
                'calc_gf_disp_tri progress: ',100*i/displacement%nrecords,'%'
            if (i.eq.displacement%nrecords) then
                write(0,*)
            endif
        endif

        do j = 1,fault%nrecords
            tri_coord(1,1) = fault%array(j,1)
            tri_coord(2,1) = fault%array(j,2)
            tri_coord(3,1) = fault%array(j,3)
            tri_coord(1,2) = fault%array(j,4)
            tri_coord(2,2) = fault%array(j,5)
            tri_coord(3,2) = fault%array(j,6)
            tri_coord(1,3) = fault%array(j,7)
            tri_coord(2,3) = fault%array(j,8)
            tri_coord(3,3) = fault%array(j,9)
            tri_coord(:,4) = 0.0d0

            if (coord_type.eq.'geographic') then
                call tri_center(center,tri_coord(:,1),tri_coord(:,2),tri_coord(:,3))
                call ddistaz(dist,az,center(1),center(2),sta_coord(1),sta_coord(2))
                dist = dist*6.371d6
                sta_coord_new(1) = dist*dsin(az)
                sta_coord_new(2) = dist*dcos(az)
                sta_coord_new(3) = sta_coord(3)
                do iTri = 1,3
                    call ddistaz(dist,az,center(1),center(2),tri_coord(1,iTri),tri_coord(2,iTri))
                    dist = dist*6.371d6
                    tri_coord_new(1,iTri) = dist*dsin(az)
                    tri_coord_new(2,iTri) = dist*dcos(az)
                    tri_coord_new(3,iTri) = tri_coord(3,iTri)
                enddo
            else
                sta_coord_new = sta_coord
                tri_coord_new = tri_coord
            endif

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 0.0d0
            endif
            slip(1) = -slip_magnitude*dcos(rak*d2r)   ! strike-slip, right lateral=positive
            slip(2) = -slip_magnitude*dsin(rak*d2r)   ! dip-slip, normal=positive
            slip(3) = 0.0d0                           ! tensile slip

            ! Strike-slip Green's function
            call tri_disloc_disp(disp, sta_coord_new, tri_coord_new, poisson, slip)
            ! Flip z-component
            disp(3) = -disp(3)
            gf_disp%array(i,                        j) = disp(1)
            gf_disp%array(i+1*displacement%nrecords,j) = disp(2)
            gf_disp%array(i+2*displacement%nrecords,j) = disp(3)

            ! Check for rake constraints; if none, calculate dip-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
                    .and.rake_constraint%nfields.eq.2) then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,2)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,2)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 90.0d0
            endif
            slip(1) = -slip_magnitude*dcos(rak*d2r)   ! strike-slip, right lateral=positive
            slip(2) = -slip_magnitude*dsin(rak*d2r)   ! dip-slip, normal=positive
            slip(3) = 0.0d0                           ! tensile slip

            ! Dip-slip (or second rake) Green's function
            call tri_disloc_disp(disp, sta_coord_new, tri_coord_new, poisson, slip)
            ! Flip z-component
            disp(3) = -disp(3)
            gf_disp%array(i,                        j+fault%nrecords) = disp(1)
            gf_disp%array(i+1*displacement%nrecords,j+fault%nrecords) = disp(2)
            gf_disp%array(i+2*displacement%nrecords,j+fault%nrecords) = disp(3)
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_disp_tri says: finished'
    endif
    if (verbosity.ge.3) then
        write(stderr,'(2A4,6A14)') 'sta','flt','dgf_ss_x','dgf_ss_y','dgf_ss_z',&
                                               'dgf_ds_x','dgf_ds_y','dgf_ds_z'
        do i = 1,displacement%nrecords
            do j = 1,fault%nrecords
                write(stderr,'(2I4,6F14.6)') i,j, &
                    gf_disp%array(i                        ,j), &
                    gf_disp%array(i+1*displacement%nrecords,j), &
                    gf_disp%array(i+2*displacement%nrecords,j), &
                    gf_disp%array(i                        ,j+fault%nrecords), &
                    gf_disp%array(i+1*displacement%nrecords,j+fault%nrecords), &
                    gf_disp%array(i+2*displacement%nrecords,j+fault%nrecords)
            enddo
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine calc_gf_disp_tri

!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_stress_tri()
    use io_module, only: stderr, verbosity
    use variable_module, only: inversion_mode, coord_type, fault, rake_constraint, halfspace, &
                               gf_stress, sts_dist
    use elast
    use tri_disloc_module, only: tri_disloc_strain, tri_center, tri_geometry, tri_geo2cart
    implicit none
    ! Local variables
    double precision :: slip(3), slip_magnitude, rak, sta_coord(3), tri_coord(3,4), &
                        tri_coord_new(3,4), sta_coord_new(3), center(3)
    double precision :: shear_modulus, lambda, poisson
    double precision :: unit_normal(3), unit_strike(3), unit_updip(3)
    double precision :: strain(3,3), stress(3,3), traction(3), traction_components(3)
    double precision :: dist, az, pt1(3), pt2(3), pt3(3)
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r=pi/180.0d0
    integer :: i, j, k, prog, iTri

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_stress_tri says: starting'
    endif
    prog = 0

    ! Initialize rake angle
    rak = 0.0d0

    ! Unit slip for GF computation
    slip_magnitude = 1.0d0

    ! Half-space array contains vp, vs, density
    shear_modulus = halfspace%array(1,2)*halfspace%array(1,2)*halfspace%array(1,3)
    lambda = halfspace%array(1,1)*halfspace%array(1,1)*halfspace%array(1,3) - 2.0d0*shear_modulus
    poisson = lambda/2.0d0/(lambda+shear_modulus)

    ! Shear stress Green's functions for each fault-fault pair
    do i = 1,fault%nrecords
        if (verbosity.ge.1) then
            write(0,'(A1,A,I3,A)',advance='no') achar(13), &
                'calc_gf_stress_tri progress: ',100*i/fault%nrecords,'%'
            if (i.eq.fault%nrecords) then
                write(0,*)
            endif
        endif
        call tri_center(sta_coord,fault%array(i,1:3),fault%array(i,4:6),fault%array(i,7:9))
        if (coord_type.eq.'cartesian') then
            call tri_geometry(unit_normal,unit_strike,unit_updip,&
                              fault%array(i,1:3),fault%array(i,4:6),fault%array(i,7:9))
        elseif (coord_type.eq.'geographic') then
            ! Triangle points: lon lat dep(m) to x y z
            call tri_geo2cart(pt1,pt2,pt3,fault%array(i,1:3),fault%array(i,4:6), &
                              fault%array(i,7:9),'m')
            call tri_geometry(unit_normal,unit_strike,unit_updip,pt1,pt2,pt3)
        endif
        ! write(0,*) 'unit_normal:',unit_normal
        ! write(0,*) 'unit_strike:',unit_strike
        ! write(0,*) 'unit_updip:',unit_updip

        do j = 1,fault%nrecords
            ! write(0,*) i,j
            tri_coord(1,1) = fault%array(j,1)
            tri_coord(2,1) = fault%array(j,2)
            tri_coord(3,1) = fault%array(j,3)
            tri_coord(1,2) = fault%array(j,4)
            tri_coord(2,2) = fault%array(j,5)
            tri_coord(3,2) = fault%array(j,6)
            tri_coord(1,3) = fault%array(j,7)
            tri_coord(2,3) = fault%array(j,8)
            tri_coord(3,3) = fault%array(j,9)
            tri_coord(:,4) = 0.0d0

            if (coord_type.eq.'geographic') then
                call tri_center(center,tri_coord(:,1),tri_coord(:,2),tri_coord(:,3))
                call ddistaz(dist,az,center(1),center(2),sta_coord(1),sta_coord(2))
                dist = dist*6.371d6
                sta_coord_new(1) = dist*dsin(az)
                sta_coord_new(2) = dist*dcos(az)
                sta_coord_new(3) = sta_coord(3)
                do iTri = 1,3
                    call ddistaz(dist,az,center(1),center(2),tri_coord(1,iTri),tri_coord(2,iTri))
                    dist = dist*6.371d6
                    tri_coord_new(1,iTri) = dist*dsin(az)
                    tri_coord_new(2,iTri) = dist*dcos(az)
                    tri_coord_new(3,iTri) = tri_coord(3,iTri)
                enddo
            else
                sta_coord_new = sta_coord
                tri_coord_new = tri_coord
            endif

            ! Set resolved tractions to zero if distance is larger than threshold
            dist = dsqrt(sta_coord_new(1)*sta_coord_new(1)+ &
                         sta_coord_new(2)*sta_coord_new(2)+ &
                        (sta_coord_new(3)-center(3))*(sta_coord(3)-center(3)))
            if (dist.gt.sts_dist) then
                gf_stress%array(i,j) = 0.0d0
                gf_stress%array(i,j+fault%nrecords) = 0.0d0
                gf_stress%array(i+fault%nrecords,j) = 0.0d0
                gf_stress%array(i+fault%nrecords,j+fault%nrecords) = 0.0d0
                cycle
            endif

            ! Avoid singular solutions with measurement point lying on fault
            do k = 1,3
                if (dabs(center(k)-sta_coord_new(k)).lt.1.0d-3) then
                    sta_coord_new(k) = center(k) + 1.0d-3
                endif
            enddo

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 0.0d0
            endif
            slip(1) = -slip_magnitude*dcos(rak*d2r)   ! strike-slip, right lateral=positive
            slip(2) = -slip_magnitude*dsin(rak*d2r)   ! dip-slip, normal=positive
            slip(3) = 0.0d0                           ! tensile slip

            ! Shear stress Green's function produced by strike-slip or fixed-rake source
            ! write(0,*) 'sta_coord:',sta_coord_new
            ! write(0,*) 'tri_coord1:',tri_coord_new(:,1)
            ! write(0,*) 'tri_coord2:',tri_coord_new(:,2)
            ! write(0,*) 'tri_coord3:',tri_coord_new(:,3)
            ! write(0,*) 'poisson:',poisson
            ! write(0,*) 'slip:',slip
            call tri_disloc_strain(strain, sta_coord_new, tri_coord_new, poisson, slip)
            ! Flip some strain components
            strain(1,3) = -strain(1,3)
            strain(3,1) = -strain(1,3)
            strain(2,3) = -strain(2,3)
            strain(3,2) = -strain(2,3)
            ! write(0,*) strain(1,1),strain(2,2),strain(3,3)
            ! write(0,*) strain(1,2),strain(1,3),strain(2,3)
            ! Calculate tractions
            call calc_strain_to_stress(strain,halfspace%array(1,1),halfspace%array(1,2), &
                                       halfspace%array(1,3),stress)
            call calc_traction(stress,unit_normal,traction)
            call calc_traction_components(traction,unit_normal,unit_strike,unit_updip, &
                                          traction_components)
            ! Load stress Green's functions array
            gf_stress%array(i               ,j               ) = -traction_components(2)
            gf_stress%array(i+fault%nrecords,j               ) = -traction_components(3)
            ! write(0,*) 'trac_ss:',-traction_components(2)
            ! write(0,*) ' trac_ds:',-traction_components(3)

            ! Shear stress Green's function produced by dip-slip source
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
                    .and.rake_constraint%nfields.eq.2) then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,2)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,2)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 90.0d0
            endif
            slip(1) = -slip_magnitude*dcos(rak*d2r)   ! strike-slip, right lateral=positive
            slip(2) = -slip_magnitude*dsin(rak*d2r)   ! dip-slip, normal=positive
            slip(3) = 0.0d0                           ! tensile slip

            ! Shear stress Green's function produced by strike-slip or fixed-rake source
            ! write(0,*) 'sta_coord:',sta_coord_new
            ! write(0,*) 'tri_coord1:',tri_coord_new(:,1)
            ! write(0,*) 'tri_coord2:',tri_coord_new(:,2)
            ! write(0,*) 'tri_coord3:',tri_coord_new(:,3)
            ! write(0,*) 'poisson:',poisson
            ! write(0,*) 'slip:',slip
            call tri_disloc_strain(strain, sta_coord_new, tri_coord_new, poisson, slip)
            strain(1,3) = -strain(1,3)
            strain(3,1) = -strain(1,3)
            strain(2,3) = -strain(2,3)
            strain(3,2) = -strain(2,3)
            ! write(0,*) strain(1,1),strain(2,2),strain(3,3)
            ! write(0,*) strain(1,2),strain(1,3),strain(2,3)
            ! Calculate tractions
            call calc_strain_to_stress(strain,halfspace%array(1,1),halfspace%array(1,2), &
                                       halfspace%array(1,3),stress)
            call calc_traction(stress,unit_normal,traction)
            call calc_traction_components(traction,unit_normal,unit_strike,unit_updip, &
                                          traction_components)
            ! Load stress Green's functions array
            gf_stress%array(i               ,j+fault%nrecords) = -traction_components(2)
            gf_stress%array(i+fault%nrecords,j+fault%nrecords) = -traction_components(3)
            ! write(0,*) 'trac_ss:',-traction_components(2)
            ! write(0,*) 'trac_ds:',-traction_components(3)
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') "calc_gf_stress_tri says: finished"
    endif
    if (verbosity.ge.3) then
        write(stderr,'(2A4,4A14)') 'sta','flt','sgf_ss_ss','sgf_ss_ds','sgf_ds_ss','sgf_ds_ds'
        do i = 1,fault%nrecords
            do j = 1,fault%nrecords
                write(stderr,'(2I4,1P4E14.6)') i,j, &
                    gf_stress%array(i               ,j), &
                    gf_stress%array(i+fault%nrecords,j), &
                    gf_stress%array(i               ,j+fault%nrecords), &
                    gf_stress%array(i+fault%nrecords,j+fault%nrecords)
            enddo
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine calc_gf_stress_tri

!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_los_tri()
    ! use io_module, only: stderr, verbosity
    ! use variable_module, only: inversion_mode, &
    !                            los, fault, gf_los, halfspace, rake_constraint
    implicit none
    ! ! Local variables
    ! integer :: i, j
    ! double precision :: slip, stlo, stla, stdp, evlo, evla, evdp, str, dip, rak, wid, len
    ! double precision :: dx, dy, dist, az, x, y, ux, uy, uz, uxp, uyp, theta, uhor, vp, vs, dens
    ! double precision :: lookaz, lookinc
    ! double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r = pi/180.0d0
	!
    ! if (verbosity.ge.2) then
    !     write(stderr,'(A)') 'calc_gf_los_okada_rect says: starting'
    ! endif
	!
    ! ! Unit slip for displacement GF computation
    ! slip = 1.0d0
	!
    ! ! Half-space array contains vp, vs, density
    ! vp = halfspace%array(1,1)
    ! vs = halfspace%array(1,2)
    ! dens = halfspace%array(1,3)
	!
    ! ! LOS displacement Green's function for each fault-station pair
    ! do i = 1,los%nrecords
    !     stlo = los%array(i,1)
    !     stla = los%array(i,2)
    !     stdp = los%array(i,3)
    !     lookaz   = los%array(i,5)*d2r
    !     lookinc  = los%array(i,6)*d2r
    !     do j = 1,fault%nrecords
    !         evlo = fault%array(j,1)
    !         evla = fault%array(j,2)
    !         evdp = fault%array(j,3)
    !         str  = fault%array(j,4)
    !         dip  = fault%array(j,5)
    !         wid  = fault%array(j,6)
    !         len  = fault%array(j,7)
	!
    !         ! Origin is at fault coordinate, x-axis points horizontal up-dip
    !         dx = stlo - evlo
    !         dy = stla - evla
    !         dist = dsqrt(dx*dx+dy*dy)
    !         az = datan2(dx,dy)
    !         x = dist*( dcos(az-d2r*str))
    !         y = dist*(-dsin(az-d2r*str))
	!
    !         ! Check for rake constraints; if none, calculate strike-slip GFs
    !         if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
    !             if (rake_constraint%nrecords.eq.1) then
    !                 rak = rake_constraint%array(1,1)
    !             elseif (rake_constraint%nrecords.eq.fault%nrecords) then
    !                 rak = rake_constraint%array(j,1)
    !             else
    !                 call usage('!! Error: incorrect number of rake constraints')
    !             endif
    !         else
    !             rak = 0.0d0
    !         endif
	!
    !         ! Strike-slip Green's function
    !         call o92rect(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,wid,len,slip,vp,vs,dens)
    !         ! Rotate back to original x-y coordinates
    !         theta = datan2(uyp,uxp)
    !         uhor = dsqrt(uxp*uxp+uyp*uyp)
    !         theta = d2r*str - theta
    !         ux = uhor*dsin(theta)
    !         uy = uhor*dcos(theta)
    !         gf_los%array(i,j) = &
    !                      ux*sin(lookaz)*cos(lookinc) + uy*cos(lookaz)*cos(lookinc) - uz*sin(lookinc)
	!
    !         ! Check for rake constraints; if none, calculate dip-slip GFs
    !         if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
    !                 .and.rake_constraint%nfields.eq.2) then
    !             if (rake_constraint%nrecords.eq.1) then
    !                 rak = rake_constraint%array(1,2)
    !             elseif (rake_constraint%nrecords.eq.fault%nrecords) then
    !                 rak = rake_constraint%array(j,2)
    !             else
    !                 call usage('!! Error: incorrect number of rake constraints')
    !             endif
    !         else
    !             rak = 90.0d0
    !         endif
	!
    !         call o92rect(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,wid,len,slip,vp,vs,dens)
    !         ! Rotate back to original x-y coordinates
    !         theta = datan2(uyp,uxp)
    !         uhor = dsqrt(uxp*uxp+uyp*uyp)
    !         theta = d2r*str - theta
    !         ux = uhor*dsin(theta)
    !         uy = uhor*dcos(theta)
    !         gf_los%array(i,j+fault%nrecords) = &
    !                      ux*sin(lookaz)*cos(lookinc) + uy*cos(lookaz)*cos(lookinc) - uz*sin(lookinc)
    !     enddo
    ! enddo
	!
    ! if (verbosity.ge.2) then
    !     write(stderr,'(A)') 'calc_gf_los_okada_rect says: finished'
    ! endif
    ! if (verbosity.ge.3) then
    !     write(stderr,'(2A4,6A14)') 'sta','flt','losgf_ss','losgf_ds'
    !     do i = 1,los%nrecords
    !         do j = 1,fault%nrecords
    !             write(stderr,'(2I4,6F14.6)') i,j, &
    !                 gf_los%array(i,j), gf_los%array(i,j+fault%nrecords)
    !         enddo
    !     enddo
    ! endif
    ! if (verbosity.ge.2) then
    !     write(stderr,*)
    ! endif
	!
    return
    end subroutine calc_gf_los_tri

end module gf_module
