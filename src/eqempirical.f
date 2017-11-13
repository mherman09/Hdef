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
      else
          call usage('!! Error: no relationship named '//trim(model))
      endif

      END

!----------------------------------------------------------------------!

      SUBROUTINE wellscopp(mom,mag)
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
          write(*,1002) typ(i),maxslp
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
      write(*,1008)
      do 106 i = 1,4
          avgslp = 10**(a(i)+b(i)*mag)
          write(*,1002) typ(i),avgslp
  106 continue
 1001 format('Surface rupture length:')
 1002 format(4X,A6,':',F12.3,X,'m')
 1003 format('Total subsurface rupture length:')
 1004 format('Downdip rupture width:')
 1005 format('Total rupture area:')
 1006 format(4X,A6,':',F12.3,X,'m^2')
 1007 format('Maximum slip:')
 1008 format('Average slip:')
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
     3                         '8.0 events)'
!      write(0,*)
!     1 '                       HANKSBAKUN (Hanks and Bakun, 2008)'
      write(0,*)
     1 '                       BLASERETAL (Blaser et al., 2010): ',
     2                         'global subduction zone events (Mw 4.8-',
     3                         '9.5)'
      write(0,*)
     1 '                       STRASSERETAL (Strasser et al., 2010): ',
     2                         'global subduction zone events (Mw 6.0-',
     3                         '9.5)'
      write(0,*)
     1 '                       YENMA (Yen and Ma, 2011): collision ',
     2                         'zone (Taiwan) events (Mw 4.6-8.9)'
!      write(0,*)
!     1 '                       STIRLINGETAL (Stirling et al., 2013)'
      write(0,*)
     1 '                       HAYES (Hayes, 2017): global finite ',
     2                         'fault models of Mw 7.5-9.2 events'
      STOP
      END
