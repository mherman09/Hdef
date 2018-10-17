module okada_module

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_disp_okada_rect()
    use io_module, only: stderr, verbosity
    use variable_module, only: inversion_mode, &
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
        do j = 1,fault%nrecords
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str  = fault%array(j,4)
            dip  = fault%array(j,5)
            wid  = fault%array(j,6)
            len  = fault%array(j,7)

            ! Origin is at fault coordinate, x-axis points horizontal up-dip
            dx = stlo - evlo
            dy = stla - evla
            dist = dsqrt(dx*dx+dy*dy)
            az = datan2(dx,dy)
            x = dist*( dcos(az-d2r*str))
            y = dist*(-dsin(az-d2r*str))

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call print_usage('!! Error: incorrect number of rake constraints')
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
                    call print_usage('!! Error: incorrect number of rake constraints')
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
    use variable_module, only: inversion_mode, fault, rake_constraint, halfspace, gf_stress
    use elast_module
    implicit none
    ! Local variables
    double precision :: slip, stlo, stla, stdp, sta_str, sta_dip
    double precision :: evlo, evla, evdp, str, dip, rak, wid, len
    double precision :: delx, dely, dist, az, x, y, vp, vs, dens
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
        if (i.gt.prog) then
            if (verbosity.ge.2) then
                write(0,'(A,I3,A)') 'calc_gf_stress_okada_rect says:',100*prog/fault%nrecords,&
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
            wid  = fault%array(j,6)
            len  = fault%array(j,7)

            ! Avoid singular solutions with measurement point lying on fault
            if (dabs(evlo-stlo).lt.1.0d0) then
                stlo = evlo + 1.0d0
            endif
            if (dabs(evla-stla).lt.1.0d0) then
                stla = evla + 1.0d0
            endif
            if (dabs(evdp-stdp).lt.1.0d0) then
                stdp = evdp + 1.0d0
            endif

            ! Rotate coordinates to x=along-strike, y=horizontal up-dip
             delx = stlo - evlo
             dely = stla - evla
             dist = dsqrt(delx*delx+dely*dely)
             az   = datan2(delx,dely) ! clockwise from north
             x = dist*( dcos(az-d2r*str))
             y = dist*(-dsin(az-d2r*str))

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call print_usage('!! Error: incorrect number of rake constraints')
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
                    call print_usage('!! Error: incorrect number of rake constraints')
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
    use variable_module, only: inversion_mode, &
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
            dx = stlo - evlo
            dy = stla - evla
            dist = dsqrt(dx*dx+dy*dy)
            az = datan2(dx,dy)
            x = dist*( dcos(az-d2r*str))
            y = dist*(-dsin(az-d2r*str))

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call print_usage('!! Error: incorrect number of rake constraints')
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

            ! Dip-slip Green's function
            rak  = 90.0d0
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

end module okada_module
