      PROGRAM main
C----
C Compute relative plate velocities given (a) a pair of plates from a
C published study or (b) an Euler pole
C----
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      CHARACTER*80 ifile,ofile,polec
      REAL*8 pole(3),stlo,stla,r(3),v(3),lo,la,ve,vn,vz
      INTEGER model,verbos

C----
C Get command line arguments
C----
      call gcmdln(ifile,ofile,polec,model,verbos)
      if (verbos.eq.1) then
          write(0,*) 'platemotion running in verbose mode'
          write(0,*)
          write(0,*) 'Command line variables'
          write(0,*) 'ifile  : ',trim(ifile)
          write(0,*) 'ofile  : ',trim(ofile)
          write(0,*) 'polec  : ',trim(polec)
          write(0,*) 'model  : ',model
      endif

C----
C Get pole (or list available plates if requested)
C----
      call getpole(pole,polec,model,verbos)

C----
C Calculate predicted velocities
C----
      if (ifile.ne.'none') open(unit=11,file=ifile,status='old')
  102 if (ifile.eq.'none') then
          read(*,*,end=101) stlo,stla
      else
          read(11,*,end=101) stlo,stla
      endif
      la = stla*d2r
      lo = stlo*d2r
      r(1) = dcos(stla*d2r)*dcos(stlo*d2r)
      r(2) = dcos(stla*d2r)*dsin(stlo*d2r)
      r(3) = dsin(stla*d2r)
      call cross(v,pole,r)
      vn = -(v(1)*dsin(la)*dcos(lo)+v(2)*dsin(la)*dsin(lo)
     1                                                   -v(3)*dcos(la))
      ve =  -v(1)         *dsin(lo)+v(2)         *dcos(lo)
      vz =   v(1)*dcos(la)*dcos(lo)+v(2)*dcos(la)*dsin(lo)
     1                                                    +v(3)*dsin(la)
      vn = vn*111.19d0
      ve = ve*111.19d0
      write(*,*) stlo,stla,ve,vn
      goto 102
  101 continue
      END

C----------------------------------------------------------------------C

      SUBROUTINE getpole(pole,polec,model,verbos)
      IMPLICIT none
      CHARACTER*80 polec,model_name
      REAL*8 pole(3),x(3)
      INTEGER i,isave,model,verbos

      if (verbos.eq.1) then
          write(0,*)
          write(0,*) 'getpole says: Getting Euler pole'
      endif

      ! Print available plates if requested
      if (model.eq.-1) then
          call morvelplates()
      elseif (model.eq.-2) then
          call nuvel1aplates()
      elseif (model.eq.-3) then
          call morvel56plates()
      elseif (model.lt.0) then
          write(0,*) '!! Error: no plate model number ',abs(model)
          call usage('!! Exiting')
      endif

      ! Parse the pole from the command line
      if (verbos.eq.1) then
          write(0,*) 'getpole says: Parsing plates/pole from input: ',
     1               trim(polec)
      endif
      if (polec.eq.'none') call usage('!! Error: no plate/pole defined')

      ! Input format is either "P1/P2" or "lon/lat/vel"
      ! Check for first "/"
      i = index(polec,'/')
      if (i.gt.0) then
          polec(i:i) = ' '
          isave = i
      else
          call usage('!! Error: input must contain at least one "/"')
      endif

      ! Check for second "/" to determine operation mode
      i = index(polec,'/')
      if (i.gt.0) then
          polec(i:i) = ' '
          i = -1
      else
          i = 1
      endif

      ! Compute Cartesian Euler pole
      if (i.gt.0) then
          if (i.eq.1) then
              model_name = 'MORVEL'
          elseif (i.eq.2) then
              model_name = 'NUVEL-1A'
          elseif (i.eq.3) then
              model_name = 'MORVEL56'
          else
              call usage('!! Error: No such model!')
          endif
          if (verbos.eq.1) then
              write(0,*) 'getpole says: Using plate circuit model: ',
     1                   trim(model_name)
          endif
          call model2pole(pole,polec,model,verbos)
      elseif (i.lt.0) then
          if (verbos.eq.1) then
              write(0,*) 'getpole says: Reading pole from command line'
          endif
          read(polec,*) pole(1),pole(2),pole(3)
      else
          call usage('!! Error in parsing pole inputs')
      endif
      call lalomag2xyz(x(1),x(2),x(3),pole(2),pole(1),pole(3))
      pole(1) = x(1)
      pole(2) = x(2)
      pole(3) = x(3)
      if (verbos.eq.1) then
          write(0,1003) pole(1),pole(2),pole(3)
      endif
 1003 format(' getpole says: Euler pole is: ',3(1PE14.6))

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE lalomag2xyz(x,y,z,lat,lon,mag)
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),d2r=pi/180.0d0)
      REAL*8 lat,lon,mag,x,y,z
      x = mag*cos(lat*d2r)*cos(lon*d2r)
      y = mag*cos(lat*d2r)*sin(lon*d2r)
      z = mag*sin(lat*d2r)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE cross(ovect,v1,v2)
      IMPLICIT none
      REAL*8 ovect(3),v1(3),v2(3)
      ovect(1) = v1(2)*v2(3) - v1(3)*v2(2)
      ovect(2) = v1(3)*v2(1) - v1(1)*v2(3)
      ovect(3) = v1(1)*v2(2) - v1(2)*v2(1)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE model2pole(pole,polec,model,verbos)
      IMPLICIT none
      CHARACTER*80 polec,plate(2)
      REAL*8 pole(3),x(4,3)
      INTEGER i,model,verbos

      ! Read plate pair
      if (verbos.eq.1) then
          write(0,*)
          write(0,*) 'model2pole says: Using plate pair: ',trim(polec)
      endif
      read(polec,*) plate(1),plate(2)

      do 101 i = 1,2
          ! MORVEL
          if (model.eq.1) then
              if (plate(i).eq.'AM') then
                  x(i,1) = -82.7
                  x(i,2) = 65.9
                  x(i,3) = 0.929
              elseif (plate(i).eq.'AN') then
                  x(i,1) = -78.5
                  x(i,2) = 65.9
                  x(i,3) = 0.887
              elseif (plate(i).eq.'AR') then
                  x(i,1) = -33.2
                  x(i,2) = 60.0
                  x(i,3) = 1.159
              elseif (plate(i).eq.'AU') then
                  x(i,1) = 6.3
                  x(i,2) = 60.1
                  x(i,3) = 1.079
              elseif (plate(i).eq.'CA') then
                  x(i,1) = -77.5
                  x(i,2) = 55.8
                  x(i,3) = 0.905
              elseif (plate(i).eq.'CO') then
                  x(i,1) = -112.8
                  x(i,2) = 42.2
                  x(i,3) = 1.676
              elseif (plate(i).eq.'CP') then
                  x(i,1) = -10.1
                  x(i,2) = 62.3
                  x(i,3) = 1.139
              elseif (plate(i).eq.'EU') then
                  x(i,1) = -78.9
                  x(i,2) = 61.3
                  x(i,3) = 0.856
              elseif (plate(i).eq.'IN') then
                  x(i,1) = -31.2
                  x(i,2) = 61.4
                  x(i,3) = 1.141
              elseif (plate(i).eq.'JF') then
                  x(i,1) = 37.8
                  x(i,2) = -0.6
                  x(i,3) = 0.625
              elseif (plate(i).eq.'LW') then
                  x(i,1) = -66.9
                  x(i,2) = 60.0
                  x(i,3) = 0.932
              elseif (plate(i).eq.'MQ') then
                  x(i,1) = -8.0
                  x(i,2) = 59.2
                  x(i,3) = 1.686
              elseif (plate(i).eq.'NA') then
                  x(i,1) = -71.7
                  x(i,2) = 48.9
                  x(i,3) = 0.750
              elseif (plate(i).eq.'NB') then
                  x(i,1) = -66.6
                  x(i,2) = 58.7
                  x(i,3) = 0.935
              elseif (plate(i).eq.'NZ') then
                  x(i,1) = -87.8
                  x(i,2) = 55.9
                  x(i,3) = 1.311
              elseif (plate(i).eq.'PS') then
                  x(i,1) = -41.9
                  x(i,2) = -4.6
                  x(i,3) = 0.890
              elseif (plate(i).eq.'RI') then
                  x(i,1) = -104.8
                  x(i,2) = 25.7
                  x(i,3) = 4.966
              elseif (plate(i).eq.'SA') then
                  x(i,1) = -77.0
                  x(i,2) = 56.0
                  x(i,3) = 0.653
              elseif (plate(i).eq.'SC') then
                  x(i,1) = -78.0
                  x(i,2) = 57.8
                  x(i,3) = 0.755
              elseif (plate(i).eq.'SM') then
                  x(i,1) = -73.5
                  x(i,2) = 59.3
                  x(i,3) = 0.980
              elseif (plate(i).eq.'SR') then
                  x(i,1) = -75.8
                  x(i,2) = 55.7
                  x(i,3) = 0.636
              elseif (plate(i).eq.'SU') then
                  x(i,1) = -78.0
                  x(i,2) = 59.8
                  x(i,3) = 0.973
              elseif (plate(i).eq.'SW') then
                  x(i,1) = -42.4
                  x(i,2) = -3.8
                  x(i,3) = 1.444
              elseif (plate(i).eq.'YZ') then
                  x(i,1) = -82.4
                  x(i,2) = 65.5
                  x(i,3) = 0.968
              elseif (plate(i).eq.'ITRF') then
                  x(i,1) = -68.2
                  x(i,2) = 63.4
                  x(i,3) = 0.677
              elseif (plate(i).eq.'PA') then
                  x(i,1) = 0.0
                  x(i,2) = 0.0
                  x(i,3) = 0.0
              else
                  write(0,*) '!! Error: no plate named ',trim(plate(i)),
     1                       ' in model ',model
                  call usage('')
              endif
          ! NUVEL-1A
          elseif (model.eq.2) then
              if (plate(i).eq.'AF') then
                  x(i,1) = -73.174
                  x(i,2) = 59.160
                  x(i,3) = 0.9270
              elseif (plate(i).eq.'AN') then
                  x(i,1) = -83.984
                  x(i,2) = 64.315
                  x(i,3) = 0.8695
              elseif (plate(i).eq.'AR') then
                  x(i,1) = -33.193
                  x(i,2) = 59.658
                  x(i,3) = 1.1107
              elseif (plate(i).eq.'AU') then
                  x(i,1) = 1.742
                  x(i,2) = 60.080
                  x(i,3) = 1.0744
              elseif (plate(i).eq.'CA') then
                  x(i,1) = -80.802
                  x(i,2) = 54.195
                  x(i,3) = 0.8160
              elseif (plate(i).eq.'CO') then
                  x(i,1) = 251.371
                  x(i,2) = 36.823
                  x(i,3) = 1.9975
              elseif (plate(i).eq.'EU') then
                  x(i,1) = -85.819
                  x(i,2) = 61.066
                  x(i,3) = 0.8591
              elseif (plate(i).eq.'IN') then
                  x(i,1) = -30.403
                  x(i,2) = 60.494
                  x(i,3) = 1.1034
              elseif (plate(i).eq.'JF') then
                  x(i,1) = 29.3
                  x(i,2) = 28.3
                  x(i,3) = 0.520
              elseif (plate(i).eq.'NA') then
                  x(i,1) = -78.167
                  x(i,2) = 48.709
                  x(i,3) = 0.7486
              elseif (plate(i).eq.'NZ') then
                  x(i,1) = -90.096
                  x(i,2) = 55.578
                  x(i,3) = 1.3599
              elseif (plate(i).eq.'PS') then
                  x(i,1) = 0.96
                  x(i,2) = -45.8
                  x(i,3) = 0.116
              elseif (plate(i).eq.'RI') then
                  x(i,1) = 257.6
                  x(i,2) = 31.0
                  x(i,3) = 2.45
              elseif (plate(i).eq.'SA') then
                  x(i,1) = -85.752
                  x(i,2) = 54.999
                  x(i,3) = 0.6365
              elseif (plate(i).eq.'SC') then
                  x(i,1) = -81.4
                  x(i,2) = 49.1
                  x(i,3) = 0.66
              elseif (plate(i).eq.'NNR') then
                  x(i,1) = -72.6
                  x(i,2) = 63.0
                  x(i,3) = 0.6411
              elseif (plate(i).eq.'PA') then
                  x(i,1) = 0.0
                  x(i,2) = 0.0
                  x(i,3) = 0.0
              else
                  write(0,*) '!! Error: no plate named ',trim(plate(i)),
     1                       ' in model ',model
                  call usage('')
              endif
          ! MORVEL56
          elseif (model.eq.3) then
              if (plate(i).eq.'AM') then
                  x(i,1) = -82.68
                  x(i,2) = 65.92
                  x(i,3) = 0.929
              elseif (plate(i).eq.'AN') then
                  x(i,1) = -78.53
                  x(i,2) = 65.92
                  x(i,3) = 0.887
              elseif (plate(i).eq.'AR') then
                  x(i,1) = -33.23
                  x(i,2) = 60.02
                  x(i,3) = 1.159
              elseif (plate(i).eq.'AU') then
                  x(i,1) = 6.33
                  x(i,2) = 60.08
                  x(i,3) = 1.079
              elseif (plate(i).eq.'CA') then
                  x(i,1) = -77.48
                  x(i,2) = 55.76
                  x(i,3) = 0.905
              elseif (plate(i).eq.'CO') then
                  x(i,1) = -112.78
                  x(i,2) = 42.18
                  x(i,3) = 1.676
              elseif (plate(i).eq.'CP') then
                  x(i,1) = -10.12
                  x(i,2) = 62.34
                  x(i,3) = 1.139
              elseif (plate(i).eq.'EU') then
                  x(i,1) = -78.87
                  x(i,2) = 61.27
                  x(i,3) = 0.856
              elseif (plate(i).eq.'IN') then
                  x(i,1) = -31.21
                  x(i,2) = 61.39
                  x(i,3) = 1.141
              elseif (plate(i).eq.'JF') then
                  x(i,1) = 37.84
                  x(i,2) = -0.62
                  x(i,3) = 0.625
              elseif (plate(i).eq.'LW') then
                  x(i,1) = -66.90
                  x(i,2) = 60.03
                  x(i,3) = 0.932
              elseif (plate(i).eq.'MQ') then
                  x(i,1) = -7.98
                  x(i,2) = 59.21
                  x(i,3) = 1.686
              elseif (plate(i).eq.'NA') then
                  x(i,1) = -71.71
                  x(i,2) = 48.89
                  x(i,3) = 0.750
              elseif (plate(i).eq.'NB') then
                  x(i,1) = -66.57
                  x(i,2) = 58.68
                  x(i,3) = 0.935
              elseif (plate(i).eq.'NZ') then
                  x(i,1) = -87.76
                  x(i,2) = 55.86
                  x(i,3) = 1.311
              elseif (plate(i).eq.'PS') then
                  x(i,1) = -41.87
                  x(i,2) = -4.63
                  x(i,3) = 0.890
              elseif (plate(i).eq.'RI') then
                  x(i,1) = -104.80
                  x(i,2) = 25.69
                  x(i,3) = 4.966
              elseif (plate(i).eq.'SA') then
                  x(i,1) = -77.03
                  x(i,2) = 55.97
                  x(i,3) = 0.653
              elseif (plate(i).eq.'SC') then
                  x(i,1) = -78.02
                  x(i,2) = 57.84
                  x(i,3) = 0.755
              elseif (plate(i).eq.'SM') then
                  x(i,1) = -73.55
                  x(i,2) = 59.27
                  x(i,3) = 0.980
              elseif (plate(i).eq.'SR') then
                  x(i,1) = -75.77
                  x(i,2) = 55.69
                  x(i,3) = 0.636
              elseif (plate(i).eq.'SU') then
                  x(i,1) = -77.96
                  x(i,2) = 59.81
                  x(i,3) = 0.973
              elseif (plate(i).eq.'SW') then
                  x(i,1) = -42.36
                  x(i,2) = -3.84
                  x(i,3) = 1.444
              elseif (plate(i).eq.'YZ') then
                  x(i,1) = -82.38
                  x(i,2) = 65.45
                  x(i,3) = 0.968
              elseif (plate(i).eq.'AS') then
                  x(i,1) = -70.76
                  x(i,2) = 74.36
                  x(i,3) = 0.648
              elseif (plate(i).eq.'AP') then
                  x(i,1) = -77.01
                  x(i,2) = 34.56
                  x(i,3) = 0.929
              elseif (plate(i).eq.'AT') then
                  x(i,1) = 9.12
                  x(i,2) = 54.82
                  x(i,3) = 1.667
              elseif (plate(i).eq.'BR') then
                  x(i,1) = -111.00
                  x(i,2) = 45.90
                  x(i,3) = 0.200
              elseif (plate(i).eq.'BS') then
                  x(i,1) = 122.56
                  x(i,2) = 13.34
                  x(i,3) = 2.248
              elseif (plate(i).eq.'BH') then
                  x(i,1) = 88.39
                  x(i,2) = 11.59
                  x(i,3) = 0.346
              elseif (plate(i).eq.'BU') then
                  x(i,1) = -76.63
                  x(i,2) = 7.86
                  x(i,3) = 2.523
              elseif (plate(i).eq.'CL') then
                  x(i,1) = -27.64
                  x(i,2) = 0.99
                  x(i,3) = 0.199
              elseif (plate(i).eq.'CR') then
                  x(i,1) = 174.43
                  x(i,2) = -12.56
                  x(i,3) = 3.609
              elseif (plate(i).eq.'EA') then
                  x(i,1) = 66.32
                  x(i,2) = 28.04
                  x(i,3) = 11.420
              elseif (plate(i).eq.'FT') then
                  x(i,1) = -178.82
                  x(i,2) = -10.12
                  x(i,3) = 4.847
              elseif (plate(i).eq.'GP') then
                  x(i,1) = 79.43
                  x(i,2) = 8.94
                  x(i,3) = 5.307
              elseif (plate(i).eq.'JZ') then
                  x(i,1) = 70.11
                  x(i,2) = 35.77
                  x(i,3) = 22.532
              elseif (plate(i).eq.'KE') then
                  x(i,1) = -1.83
                  x(i,2) = 47.61
                  x(i,3) = 2.832
              elseif (plate(i).eq.'MN') then
                  x(i,1) = 150.46
                  x(i,2) = -3.04
                  x(i,3) = 51.300
              elseif (plate(i).eq.'MO') then
                  x(i,1) = 79.96
                  x(i,2) = 57.43
                  x(i,3) = 0.918
              elseif (plate(i).eq.'MA') then
                  x(i,1) = 144.24
                  x(i,2) = 39.19
                  x(i,3) = 1.319
              elseif (plate(i).eq.'MS') then
                  x(i,1) = -56.78
                  x(i,2) = 10.54
                  x(i,3) = 3.915
              elseif (plate(i).eq.'NH') then
                  x(i,1) = -12.00
                  x(i,2) = 13.00
                  x(i,3) = 2.700
              elseif (plate(i).eq.'NI') then
                  x(i,1) = 169.62
                  x(i,2) = 6.95
                  x(i,3) = 3.248
              elseif (plate(i).eq.'ND') then
                  x(i,1) = -80.24
                  x(i,2) = 59.68
                  x(i,3) = 0.716
              elseif (plate(i).eq.'NBi') then
                  x(i,1) = 139.00
                  x(i,2) = -4.00
                  x(i,3) = 0.330
              elseif (plate(i).eq.'OK') then
                  x(i,1) = -76.20
                  x(i,2) = 55.81
                  x(i,3) = 0.844
              elseif (plate(i).eq.'ON') then
                  x(i,1) = 141.58
                  x(i,2) = 49.30
                  x(i,3) = 2.743
              elseif (plate(i).eq.'PM') then
                  x(i,1) = -88.73
                  x(i,2) = 55.66
                  x(i,3) = 0.906
              elseif (plate(i).eq.'SL') then
                  x(i,1) = -92.39
                  x(i,2) = 65.24
                  x(i,3) = 0.870
              elseif (plate(i).eq.'SS') then
                  x(i,1) = 133.82
                  x(i,2) = 19.26
                  x(i,3) = 1.509
              elseif (plate(i).eq.'SB') then
                  x(i,1) = -32.99
                  x(i,2) = 10.61
                  x(i,3) = 8.440
              elseif (plate(i).eq.'TI') then
                  x(i,1) = 113.28
                  x(i,2) = 15.62
                  x(i,3) = 1.629
              elseif (plate(i).eq.'TO') then
                  x(i,1) = 2.57
                  x(i,2) = 28.82
                  x(i,3) = 9.303
              elseif (plate(i).eq.'WL') then
                  x(i,1) = 131.23
                  x(i,2) = 21.81
                  x(i,3) = 1.578
              elseif (plate(i).eq.'ITRF') then
                  x(i,1) = -65.30
                  x(i,2) = 63.58
                  x(i,3) = 0.651
              elseif (plate(i).eq.'PA') then
                  x(i,1) = 0
                  x(i,2) = 0
                  x(i,3) = 0
              else
                  write(0,*) '!! Error: no plate named ',trim(plate(i)),
     1                       ' in model ',model
                  call usage('')
              endif
          endif
  101 continue
      if (verbos.eq.1) then
          write(0,1001) trim(plate(1)),x(1,1),x(1,2),x(1,3)
          write(0,1001) trim(plate(2)),x(2,1),x(2,2),x(2,3)
      endif
 1001 format(' model2pole says: pole for PA/',A,' is: ',F10.3,'E',F10.3,
     1       'N',F10.5,'deg/Ma')

      ! Convert pole lon/lat/vel to Cartesian
      call lalomag2xyz(x(3,1),x(3,2),x(3,3),x(1,2),x(1,1),x(1,3))
      call lalomag2xyz(x(4,1),x(4,2),x(4,3),x(2,2),x(2,1),x(2,3))

      ! Add poles to get pole between plates of interest
      x(1,1) = x(4,1) - x(3,1)
      x(1,2) = x(4,2) - x(3,2)
      x(1,3) = x(4,3) - x(3,3)

      ! Convert back to lon/lat/vel
      pole(3) = dsqrt(x(1,1)*x(1,1)+x(1,2)*x(1,2)+x(1,3)*x(1,3))
      pole(2) = dacos(x(1,3)/pole(3))
      pole(2) = 90.0d0 - pole(2)*1.8d2/3.14159265d0
      pole(1) = datan2(x(1,2),x(1,1))*1.8d2/3.14159265d0
      
      if (verbos.eq.1) then
          write(0,1002) trim(plate(1)),trim(plate(2)),pole(1),pole(2),
     1                  pole(3)
      endif
 1002 format(' model2pole says: pole for ',A,'/',A,' is: ',F10.3,'E',
     1       F10.3,'N',F10.5,'deg/Ma')

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE morvelplates()
      IMPLICIT none
      write(0,*) 'MORVEL plates (DeMets et al., 2010)'
      write(0,*) '    AM: Amur'
      write(0,*) '    AN: Antarctica'
      write(0,*) '    AR: Arabia'
      write(0,*) '    AU: Australia'
      write(0,*) '    CA: Caribbean'
      write(0,*) '    CO: Cocos'
      write(0,*) '    CP: Capricorn'
      write(0,*) '    EU: Eurasia'
      write(0,*) '    IN: India'
      write(0,*) '    ITRF: International reference frame'
      write(0,*) '    JF: Juan de Fuca'
      write(0,*) '    LW: Lwandle'
      write(0,*) '    MQ: Macquarie'
      write(0,*) '    NA: North America'
      write(0,*) '    NB: Nubia'
      write(0,*) '    NZ: Nazca'
      write(0,*) '    PA: Pacific'
      write(0,*) '    PS: Philippine Sea'
      write(0,*) '    RI: Rivera'
      write(0,*) '    SA: South America'
      write(0,*) '    SC: Scotia'
      write(0,*) '    SM: Somalia'
      write(0,*) '    SR: Sur'
      write(0,*) '    SU: Sunda'
      write(0,*) '    SW: Sandwich'
      write(0,*) '    YZ: Yangtze'
      STOP
      END

C----------------------------------------------------------------------C

      SUBROUTINE morvel56plates()
      IMPLICIT none
      write(0,*) 'MORVEL56 plates (Argus et al., 2011)'
      write(0,*) '    AM: Amur'
      write(0,*) '    AN: Antarctica'
      write(0,*) '    AR: Arabia'
      write(0,*) '    AS: Aegean Sea'
      write(0,*) '    AP: Altiplano'
      write(0,*) '    AT: Anatolia'
      write(0,*) '    AU: Australia'
      write(0,*) '    BH: Birds Head'
      write(0,*) '    BR: Balmoral Reef'
      write(0,*) '    BS: Banda Sea'
      write(0,*) '    BU: Burma'
      write(0,*) '    CA: Caribbean'
      write(0,*) '    CL: Caroline'
      write(0,*) '    CO: Cocos'
      write(0,*) '    CP: Capricorn'
      write(0,*) '    CR: Conway Reef'
      write(0,*) '    EA: Easter'
      write(0,*) '    EU: Eurasia'
      write(0,*) '    FT: Futuna'
      write(0,*) '    GP: Galapagos'
      write(0,*) '    IN: India'
      write(0,*) '    ITRF: International reference frame'
      write(0,*) '    JF: Juan de Fuca'
      write(0,*) '    JZ: Juan Fernandez'
      write(0,*) '    KE: Kermadec'
      write(0,*) '    LW: Lwandle'
      write(0,*) '    MA: Mariana'
      write(0,*) '    MN: Manus'
      write(0,*) '    MO: Maoke'
      write(0,*) '    MQ: Macquarie'
      write(0,*) '    MS: Molucca Sea'
      write(0,*) '    NA: North America'
      write(0,*) '    NB: Nubia'
      write(0,*) '    NBi: North Bismarck'
      write(0,*) '    ND: North Andes'
      write(0,*) '    NH: New Hebrides'
      write(0,*) '    NI: Niuafo’ou'
      write(0,*) '    NZ: Nazca'
      write(0,*) '    OK: Okhotsk'
      write(0,*) '    ON: Okinawa'
      write(0,*) '    PA: Pacific'
      write(0,*) '    PM: Panama'
      write(0,*) '    PS: Philippine Sea'
      write(0,*) '    RI: Rivera'
      write(0,*) '    SA: South America'
      write(0,*) '    SB: South Bismarck'
      write(0,*) '    SC: Scotia'
      write(0,*) '    SL: Shetland'
      write(0,*) '    SM: Somalia'
      write(0,*) '    SR: Sur'
      write(0,*) '    SS: Solomon Sea'
      write(0,*) '    SU: Sunda'
      write(0,*) '    SW: Sandwich'
      write(0,*) '    TI: Timor'
      write(0,*) '    TO: Tonga'
      write(0,*) '    WL: Woodlark'
      write(0,*) '    YZ: Yangtze'
      STOP
      END

C----------------------------------------------------------------------C

      SUBROUTINE nuvel1aplates()
      IMPLICIT none
      write(0,*) 'NUVEL-1A plates (DeMets et al., 1994)'
      write(0,*) '    AF: Africa'
      write(0,*) '    AN: Antarctica'
      write(0,*) '    AR: Arabia'
      write(0,*) '    AU: Australia'
      write(0,*) '    CA: Caribbean'
      write(0,*) '    CO: Cocos'
      write(0,*) '    EU: Eurasia'
      write(0,*) '    IN: India'
      write(0,*) '    NNR: No net rotation'
      write(0,*) '    JF: Juan de Fuca'
      write(0,*) '    NA: North America'
      write(0,*) '    NZ: Nazca'
      write(0,*) '    PA: Pacific'
      write(0,*) '    PS: Philippine Sea'
      write(0,*) '    RI: Rivera'
      write(0,*) '    SA: South America'
      write(0,*) '    SC: Scotia'
      STOP
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(ifile,ofile,pole,model,verbos)
      IMPLICIT none
      CHARACTER*80 tag,ifile,ofile,pole
      INTEGER i,narg,model,verbos
      ifile = 'none'
      ofile = 'none'
      pole  = 'none'
      model = 1
      verbos = 0
      narg = iargc()
      if (narg.eq.0) call usage('')
      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
          call getarg(i,tag)
          if (tag(1:5).eq.'-pole'.or.tag(1:6).eq.'-plate') then
              i = i + 1
              call getarg(i,pole)
          elseif (tag(1:6).eq.'-model') then
              if(index(tag,'?').ne.0) then
                  model = -1
              endif
              i = i + 1
              call getarg(i,tag)
              if (tag.eq.'MORVEL'.or.tag.eq.'1') then
                  model = 1*model
              elseif (tag.eq.'NUVEL1A'.or.tag.eq.'NUVEL-1A'.or.
     1                tag.eq.'2') then
                  model = 2*model
              elseif (tag.eq.'MORVEL56'.or.tag.eq.'3') then
                  model = 3*model
              endif
          elseif (tag(1:4).eq.'-ver') then
              verbos = 1
          elseif (tag(1:2).eq.'-h') then
              call usage('')
          elseif (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          endif
      goto 101
  102 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage(str)
      IMPLICIT none
      CHARACTER str*(*)
      if(trim(str).ne.'') write(*,*) trim(str)
      write(*,*)
     1 'Usage: platemotion -plate P1/P2 | -pole lon/lat/vel ',
     2                       '[-model[?] OPTION] [-f IFILE] ',
     3                       '[-o OFILE] [-verbos] [-h]'
      write(*,*)
     1 '  -plate P1/P2      Define plates (P2 moving w.r.t. fixed P1)'
      write(*,*)
     1 '  -pole LON/LAT/VEL Define Euler pole (velocity in deg/Ma)'
      write(*,*)
     1 '  -model[?] OPTION  Define which plate velocity model to use ',
     2                      '(default: MORVEL)'
      write(*,*)
     1 '                        1. MORVEL (DeMets et al., 2010)'
      write(*,*)
     1 '                        2. NUVEL-1A  (DeMets et al., 1994)'
      write(*,*)
     1 '                        3. MORVEL56 (Argus et al., 2011)'
      write(*,*)
     1 '                        To see available plates, use -model?'
      write(*,*)
     1 '  -f IFILE          Input file (default: run with stdin)'
      write(*,*)
     1 '                        lon lat'
      write(*,*)
     1 '  -o OFILE          Output file (default: print to stdout)'
      write(*,*)
     1 '                        v_east v_north (mm/yr)'
      write(*,*)
     1 '  -verbose          Turn on verbose mode (for debugging)'
      write(*,*)
     1 '  -h                Online help'
      write(*,*)
      STOP
      END

