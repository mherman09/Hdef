      PROGRAM main
!----
! Use empirical relations to infer fault size/slip from magnitude/moment
!----
      IMPLICIT none
      CHARACTER*200 model
      REAL*8 mom,mag
      INTEGER verbos

!----
! Parse command line inputs
!----
      call gcmdln(mom,mag,model,verbos)
      if (verbos.eq.1) then
          write(0,*) 'mom: ',mom
          write(0,*) 'mag: ',mag
          write(0,*) 'model: ',trim(model)
      endif
      if (model.eq.'none') then
          call usage('!! Error: no scaling relation defined')
      endif
      if (mom.lt.-1.0d9.and.mag.lt.-1.0d9) then
          call usage('!! Error: specify event magnitude or moment')
      endif
      if (mom.gt.-1.0d9) then
          mag = 2.0d0/3.0d0*dlog10(mom*1e7) - 10.7
          if (verbos.eq.1) then
              write(0,*)
              write(0,*) 'Converting moment to magnitude'
              write(0,*) 'Moment = ',mom,'---> Magnitude ',mag
          endif
      endif
      mom = 10**(1.5d0*(mag+10.7))/1e7

!----
! Use requested relation
!----
      if (model.eq.'WELLSCOPP') then
          call wellscopp(mom,mag)
      elseif (model.eq.'MAIBEROZA') then
          call maiberoza(mom,mag)
      elseif (model.eq.'ALLENHAYES') then
          call allenhayes(mom,mag)
      elseif (model.eq.'YENMA') then
          call yenma(mom,mag)
      else
          call usage('!! Error: no relationship named '//trim(model))
      endif

      END

!----------------------------------------------------------------------!

      SUBROUTINE wellscopp(mom,mag)
! Wells, D.L., Coppersmith, K.J. (1994). New empirical relationships
! among magnitude, rupture length, rupture width, rupture area, and
! surface displacement. Bulletin of the Seismological Society of America
! 84, 974-1002.
      IMPLICIT none
      REAL*8 mom,mag,a(4),b(4),surlen,totlen,width,area,maxslp,avgslp
      CHARACTER*3 typ(4)
      INTEGER i
      write(*,1000) mom,mag
 1000 format('For moment =',X,1PE10.3,X,'Nm,'4X,'magnitude =',X,0PF6.2)

      typ(1) = 'ss'
      typ(2) = 'rev'
      typ(3) = 'nor'
      typ(4) = 'all'
      ! Surface rupture length
      a(1) = -3.55d0  ! +/-0.37
      a(2) = -2.86d0  ! 0.55
      a(3) = -2.01d0  ! 0.65
      a(4) = -3.22d0  ! 0.27
      b(1) = 0.74d0   ! 0.05
      b(2) = 0.63d0   ! 0.08
      b(3) = 0.50d0   ! 0.10
      b(4) = 0.69d0   ! 0.04
      write(*,1001)
      do 101 i = 1,4
          surlen = 10**(a(i)+b(i)*mag)
          write(*,1002) typ(i),surlen
  101 continue
      ! Total subsurface rupture length
      a(1) = -2.57d0  ! 0.12
      a(2) = -2.42d0  ! 0.21
      a(3) = -1.88d0  ! 0.37
      a(4) = -2.44d0  ! 0.11
      b(1) = 0.62d0   ! 0.02
      b(2) = 0.58d0   ! 0.03
      b(3) = 0.50d0   ! 0.06
      b(4) = 0.59d0   ! 0.02
      write(*,1003)
      do 102 i = 1,4
          totlen = 10**(a(i)+b(i)*mag)
          write(*,1002) typ(i),totlen
  102 continue
      ! Downdip rupture width
      a(1) = -0.76d0  ! 0.12
      a(2) = -1.61d0  ! 0.20
      a(3) = -1.14d0  ! 0.28
      a(4) = -1.01d0  ! 0.10
      b(1) = 0.27d0   ! 0.02
      b(2) = 0.41d0   ! 0.03
      b(3) = 0.35d0   ! 0.05
      b(4) = 0.32d0   ! 0.02
      write(*,1004)
      do 103 i = 1,4
          width = 10**(a(i)+b(i)*mag)
          write(*,1002) typ(i),width
  103 continue
      ! Total rupture area
      a(1) = -3.42d0  ! 0.18
      a(2) = -3.99d0  ! 0.36
      a(3) = -2.87d0  ! 0.50
      a(4) = -3.49d0  ! 0.16
      b(1) = 0.90d0   ! 0.03
      b(2) = 0.98d0   ! 0.06
      b(3) = 0.82d0   ! 0.08
      b(4) = 0.91d0   ! 0.03
      write(*,1005)
      do 104 i = 1,4
          area = 10**(a(i)+b(i)*mag)
          write(*,1006) typ(i),area
  104 continue
      ! Maximum slip
      a(1) = -7.03d0  ! 0.55
      a(2) = -1.84d0  ! 1.14
      a(3) = -5.90d0  ! 1.18
      a(4) = -5.46d0  ! 0.51
      b(1) = 1.03d0   ! 0.08
      b(2) = 0.29d0   ! 0.17
      b(3) = 0.89d0   ! 0.18
      b(4) = 0.82d0   ! 0.08
      write(*,1007)
      do 105 i = 1,4
          maxslp = 10**(a(i)+b(i)*mag)
          write(*,1008) typ(i),maxslp
  105 continue
      ! Average slip
      a(1) = -6.32d0  ! 0.61
      a(2) = -0.74d0  ! 1.40
      a(3) = -4.45d0  ! 1.59
      a(4) = -4.80d0  ! 0.57
      b(1) = 0.90d0   ! 0.09
      b(2) = 0.08d0   ! 0.21
      b(3) = 0.63d0   ! 0.24
      b(4) = 0.69d0   ! 0.08
      write(*,1009)
      do 106 i = 1,4
          avgslp = 10**(a(i)+b(i)*mag)
          write(*,1008) typ(i),avgslp
  106 continue
 1001 format('Surface rupture length:')
 1002 format(4X,A6,':',F12.3,X,'km')
 1003 format('Total subsurface rupture length:')
 1004 format('Downdip rupture width:')
 1005 format('Total rupture area:')
 1006 format(4X,A6,':',F12.3,X,'km^2')
 1007 format('Maximum slip:')
 1008 format(4X,A6,':',F12.3,X,'m')
 1009 format('Average slip:')
      RETURN
      END

!----------------------------------------------------------------------!

      SUBROUTINE maiberoza(mom,mag)
! Mai, P.M., Beroza, G.C. (2000). Source scaling properties from finite-
! fault-rupture models. Bulletin of the Seismological Society of America
! 90, 604-615.
      IMPLICIT none
      REAL*8 mom,mag,a(3),b(3),efflen,effwid,effarea,effslp
      CHARACTER*3 typ(3)
      INTEGER i
      write(*,1100) mom,mag
 1100 format('For moment =',X,1PE10.3,X,'Nm,'4X,'magnitude =',X,0PF6.2)

      typ(1) = 'ss'
      typ(2) = 'ds'
      typ(3) = 'all'
      ! Effective rupture length
      a(1) = -6.31d0  ! +/-1.15
      a(2) = -6.39d0  ! 1.04
      a(3) = -6.13d0  ! 0.74
      b(1) = 0.40d0   ! 0.06
      b(2) = 0.40d0   ! 0.05
      b(3) = 0.39d0   ! 0.04
      write(*,1101)
      do 101 i = 1,3
          efflen = 10**(a(i)+b(i)*dlog10(mom))
          write(*,1102) typ(i),efflen
  101 continue
      ! Effective rupture width
      a(1) = -2.18d0  ! 1.09
      a(2) = -5.51d0  ! 0.81
      a(3) = -5.05d0  ! 0.74
      b(1) = 0.17d0   ! 0.06
      b(2) = 0.35d0   ! 0.04
      b(3) = 0.32d0   ! 0.04
      write(*,1103)
      do 102 i = 1,3
          effwid = 10**(a(i)+b(i)*dlog10(mom))
          write(*,1102) typ(i),effwid
  102 continue
      ! Effective rupture area
      a(1) = -8.49d0  ! 1.85
      a(2) = -11.90d0  ! 1.67
      a(3) = -11.18d0  ! 1.18
      b(1) = 0.57d0   ! 0.10
      b(2) = 0.75d0   ! 0.09
      b(3) = 0.72d0   ! 0.06
      write(*,1104)
      do 103 i = 1,3
          effarea = 10**(a(i)+b(i)*dlog10(mom))
          write(*,1105) typ(i),effarea
  103 continue
      ! Effective slip
      a(1) = -6.03d0  ! 1.85
      a(2) = -2.62d0  ! 1.67
      a(3) = -3.34d0  ! 1.18
      b(1) = 0.43d0   ! 0.10
      b(2) = 0.25d0   ! 0.09
      b(3) = 0.29d0   ! 0.06
      write(*,1106)
      do 104 i = 1,3
          effslp = 10**(a(i)+b(i)*dlog10(mom))
          effslp = effslp/1.0d2 ! relation is in cm, convert to m
          write(*,1107) typ(i),effslp
  104 continue
 1101 format('Effective rupture length:')
 1102 format(4X,A6,':',F12.3,X,'km')
 1103 format('Effective rupture width:')
 1104 format('Effective rupture area:')
 1105 format(4X,A6,':',F12.3,X,'km^2')
 1106 format('Effective slip:')
 1107 format(4X,A6,':',F12.3,X,'m')
      RETURN
      END

!----------------------------------------------------------------------!

      SUBROUTINE yenma(mom,mag)
! Yen, Y.T., Ma, K.F. (2011). Source-scaling relationship for M 4.6-8.9
! earthquakes, specifically for earthquakes in the collision zone of
! Taiwan. Bulletin of the Seismological Society of America 101, 464–481.
      IMPLICIT none
      REAL*8 mom,mag,a(3),b(3),len,wid,area,slip
      CHARACTER*3 typ(3)
      INTEGER i
      write(*,1000) mom,mag
 1000 format('For moment =',X,1PE10.3,X,'Nm,'4X,'magnitude =',X,0PF6.2)

      typ(1) = 'all'
      typ(2) = 'ds'
      typ(3) = 'ss'
      ! rupture length
      a(1) = -7.46 ! +/- 0.77
      a(2) = -6.66 ! 1.05
      a(3) = -8.11 ! 1.20
      b(1) = 0.47 ! 0.04
      b(2) = 0.42 ! 0.06
      b(3) = 0.50 ! 0.07
      write(*,1301)
      do 131 i = 1,3
          len = 10**(a(i)+b(i)*dlog10(mom))
          write(*,1302) typ(i),len
  131 continue
      ! rupture width
      a(1) = -6.30 ! 0.90
      a(2) = -5.76 ! 1.30
      a(3) = -6.67 ! 1.38
      b(1) = 0.40 ! 0.05
      b(2) = 0.37 ! 0.07
      b(3) = 0.42 ! 0.08
      write(*,1303)
      do 132 i = 1,3
          wid = 10**(a(i)+b(i)*dlog10(mom))
          write(*,1302) typ(i),wid
  132 continue
      ! rupture area
      a(1) = -13.79 ! 1.63
      a(2) = -12.45 ! 2.32
      a(3) = -14.77 ! 2.46
      b(1) = 0.87 ! 0.09
      b(2) = 0.80 ! 0.13
      b(3) = 0.92 ! 0.14
      write(*,1304)
      do 133 i = 1,3
          area = 10**(a(i)+b(i)*dlog10(mom))
          write(*,1305) typ(i),area
  133 continue
      ! average slip
      a(1) = -0.65 ! 1.64
      a(2) = -2.01 ! 2.32
      a(3) = 0.36 ! 2.47
      b(1) = 0.13 ! 0.09
      b(2) = 0.20 ! 0.13
      b(3) = 0.08 ! 0.14
      write(*,1306)
      do 134 i = 1,3
          slip = 10**(a(i)+b(i)*dlog10(mom))
          slip = slip/1.0d2 ! relation is in cm, convert to m
          write(*,1307) typ(i),slip
  134 continue
 1301 format('Rupture length:')
 1302 format(4X,A6,':',F12.3,X,'km')
 1303 format('Rupture width:')
 1304 format('Rupture area:')
 1305 format(4X,A6,':',F12.3,X,'km^2')
 1306 format('Average slip:')
 1307 format(4X,A6,':',F12.3,X,'m')
      RETURN
      END

!----------------------------------------------------------------------!

      SUBROUTINE allenhayes(mom,mag)
! Allen, T.I., Hayes, G.P. (2017). Alternative rupture‐scaling
! relationships for subduction interface and other offshore
! environments. Bulletin of the Seismological Society of America 107,
! 1240–1253.
      IMPLICIT none
      REAL*8 mom,mag,a(1),b(1),len,wid1,wid2,area1,area2,maxslp,avgslp
      CHARACTER*3 typ(1)
      INTEGER i
      write(*,1000) mom,mag
 1000 format('For moment =',X,1PE10.3,X,'Nm,'4X,'magnitude =',X,0PF6.2)

      typ(1) = 'thr'
      ! rupture length
      a(1) = -2.90d0
      b(1) = 0.63d0
      write(*,1201)
      do 101 i = 1,1
          len = 10**(a(i)+b(i)*mag)
          write(*,1202) typ(i),len
  101 continue
      ! rupture width (linear regression)
      a(1) = -0.86d0
      b(1) = 0.35d0
      write(*,1203)
      do 102 i = 1,1
          wid1 = 10**(a(i)+b(i)*mag)
          write(*,1202) typ(i),wid1
  102 continue
      ! rupture width (bilinear regression w saturation)
      if (mag.le.8.67) then
          a(1) = -1.91d0
          b(1) = 0.48d0
      else
          a(1) = 2.29d0
          b(1) = 0.0d0
      endif
      write(*,1204)
      do 103 i = 1,1
          wid2 = 10**(a(i)+b(i)*mag)
          write(*,1202) typ(i),wid2
  103 continue
      ! rupture area (linear regression)
      a(1) = -3.63d0
      b(1) = 0.96d0
      write(*,1205)
      do 104 i = 1,1
          area1 = 10**(a(i)+b(i)*mag)
          write(*,1206) typ(i),area1
  104 continue
      ! rupture area (bilinear regression w saturation)
      if (mag.le.8.67) then
          a(1) = -5.62d0
          b(1) = 1.22d0
      else
          a(1) = 2.23d0
          b(1) = 0.31d0
      endif
      write(*,1207)
      do 105 i = 1,1
          area2 = 10**(a(i)+b(i)*mag)
          write(*,1206) typ(i),area2
  105 continue
      ! maximum slip
      a(1) = -4.94d0
      b(1) = 0.71d0
      write(*,1208)
      do 106 i = 1,1
          maxslp = 10**(a(i)+b(i)*mag)
          write(*,1209) typ(i),maxslp
  106 continue
      ! average slip
      a(1) = -5.05d0
      b(1) = 0.66d0
      write(*,1210)
      do 107 i = 1,1
          avgslp = 10**(a(i)+b(i)*mag)
          write(*,1209) typ(i),avgslp
  107 continue
 1201 format('Rupture length:')
 1202 format(4X,A6,':',F12.3,X,'km')
 1203 format('Rupture width (linear):')
 1204 format('Rupture width (bilinear: saturates at Mw 8.67):')
 1205 format('Rupture area (linear):')
 1206 format(4X,A6,':',F12.3,X,'km^2')
 1207 format('Rupture area (bilinear: saturates at Mw 8.67):')
 1208 format('Maximum slip:')
 1209 format(4X,A6,':',F12.3,X,'m')
 1210 format('Average slip:')
      RETURN
      END

!----------------------------------------------------------------------!

      SUBROUTINE gcmdln(mom,mag,model,verbos)
      IMPLICIT none
      CHARACTER*200 model,tag
      INTEGER i,narg,verbos
      REAL*8 mom,mag
      model = 'none'
      mom = -1.0d10
      mag = -1.0d10
      verbos = 0
      narg = iargc()
      if (narg.eq.0) call usage('')
      i = 0
  901 i = i + 1
      if (i.gt.narg) goto 902
          call getarg(i,tag)
          if (tag(1:4).eq.'-mom') then
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F20.0)') mom
          elseif (tag(1:4).eq.'-mag') then
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F20.0)') mag
          elseif (tag(1:6).eq.'-model') then
              i = i + 1
              call getarg(i,model)
          elseif(tag(1:2).eq.'-v') then
              verbos = 1
          else
              call usage('!! Error: no option '//trim(tag))
          endif
          goto 901
  902 continue
      RETURN
      END

!----------------------------------------------------------------------!

      SUBROUTINE usage(str)
      IMPLICIT none
      CHARACTER str*(*)
      if (trim(str).ne.'') then
          write(0,*) trim(str)
          write(0,*)
      endif
      write(0,*)
     1 'Usage: eqempirical -mom|-mag ARG -model OPTION'
      write(0,*)
     1 '-mom|-mag ARG      Event size: ARG is either moment (Nm) ',
     2                     'or magnitude'
      write(0,*)
     1 '-model OPTION      Define which empirical relation to use'
      write(0,*)
     1 '                       WELLSCOPP (Wells and Coppersmith, ',
     2                         '1994): global continental ',
     3                         'interplate and intraplate (Mw 4.5+)'
      write(0,*)
     1 '                       MAIBEROZA (Mai and Beroza, 2000): ',
     2                         'global finite fault models (Mw 5.5-',
     3                         '8.0)'
!      write(0,*)
!     1 '                       HANKSBAKUN (Hanks and Bakun, 2008)'
      write(0,*)
     1 '                       *BLASERETAL (Blaser et al., 2010): ',
     2                         'global subduction zone (Mw 4.8-',
     3                         '9.5)'
      write(0,*)
     1 '                       *STRASSERETAL (Strasser et al., 2010): ',
     2                         'global subduction zone (Mw 6.0-',
     3                         '9.5)'
      write(0,*)
     1 '                       YENMA (Yen and Ma, 2011): Taiwan ',
     2                         'collision zone (Mw 4.6-8.9)'
!      write(0,*)
!     1 '                       STIRLINGETAL (Stirling et al., 2013)'
      write(0,*)
     1 '                       ALLENHAYES (Allen and Hayes, 2017): ',
     2                         'global finite fault models, ',
     3                         'subduction interface (Mw 7.1-9.2)'
      write(0,*)
     1 '                       * not yet implemented'
      STOP
      END
