      PROGRAM main
      IMPLICIT none
      REAL*8 x,y,rlon4,rlat4,rx4,ry4
      INTEGER iway,utmzon
      call gcmdln(iway,utmzon)
      if (utmzon.eq.0) then
          call usage('!! Error: must specify UTM zone')
      endif
  101 read (*,*,end=102) x,y
          if (iway.eq.1) then
              rx4 = x
              ry4 = y
              call utmgeo(rlon4,rlat4,rx4,ry4,utmzon,iway)
              print *,rlon4,rlat4
          elseif (iway.eq.2) then
              rlon4 = x
              rlat4 = y
              call utmgeo(rlon4,rlat4,rx4,ry4,utmzon,iway)
              print *,rx4,ry4
          endif
          goto 101
  102 continue
      END

      SUBROUTINE gcmdln(iway,utmzon)
      IMPLICIT NONE
      CHARACTER*40 tag
      INTEGER i,narg,iway,utmzon
      iway = 0
      utmzon = 0
      narg = iargc()
      if (narg.eq.0) then
          call usage('!! Error: no command line arguments specified')
      endif
      i = 0
  201 i = i + 1
      if (i.gt.narg) goto 202
          call getarg(i,tag)
          if (tag(1:8).eq.'-utm2geo') then
              iway = 1
          elseif (tag(1:8).eq.'-geo2utm') then
              iway = 2
          elseif (tag(1:7).eq.'-utmzon') then
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,I10)') utmzon
          else
              call usage('!! Error: No option '//tag)
          endif
          goto 201
  202 continue
      RETURN
      END

      SUBROUTINE usage(str)
      IMPLICIT NONE
      INTEGER lstr
      CHARACTER str*(*)
      if (str.ne.' ') then
          lstr = len(str)
          write(*,*) str(1:lstr)
          write(*,*)
      endif
      write(*,*)
     1 'Usage: utm2geo -utm2geo|-geo2utm -utmzon ZONE'
      write(*,*)
      STOP
      END
