      PROGRAM MAIN
C----
C Use empirical relations to convert from a point source
C given earthquake magnitude to a finite rectangular source
C with length, width, and slip.
C---- 
      IMPLICIT none
      INTEGER stin,stout,sterr
      PARAMETER (stin=5,stout=6,sterr=0)
      INTEGER i,nflt,typ,pf,batch
      REAL*8 lon,lat,dep,str,dip,rak,wid,len,slip,mag,mom,area,mod
      REAL*8 shrmod,ans(20,3)
      CHARACTER*30 ifile,ofile,typc
      INTEGER u,p
      LOGICAL ex

      call gcmdln(ifile,ofile,shrmod,u,p)
C     check for shear modulus
      if (shrmod.le.0.0d0) then
          write(*,*) '!! Error: No shear modulus specified'
          call usage('!! Use -mod SHR_MOD to enter shear modulus')
      endif
      shrmod = shrmod*1e9
      if (u.eq.1) goto 101
C     check for input file
      if (ifile.eq.'none') then
          write(*,*) '!! Error: Input file unspecified'
          call usage('!! Use -f IFILE to specify input file')
      else
          inquire(file=ifile,EXIST=ex)
          if (.not.ex) call usage('!! Error: no input file: '//ifile)
      endif
C     check for output file
      if (ofile.eq.'none'.and.p.eq.0) then
          write(*,*) '!! Error: No output file specified'
          call usage('!! Use -o OFILE to specify output file')
      endif

C----
C Read magnitude, fault type from file, convert to finite fault params
C----
      open(unit=11,file=ifile,status='old')
      if (p.eq.0) then
          open(unit=12,file=ofile,status='unknown')
      else
          write(sterr,9998) 'length(km)','width(km)','slip(m)'
      endif
  102 read(11,*,end=103) mag,typc
          if (typc.eq.'ss') then
              typ = 1
          elseif (typc.eq.'no') then
              typ = 3
          elseif (typc.eq.'th'.or.typc.eq.'rv') then
              typ = 2
          else
              write(*,*) '!! Error: fault type incorrect; check input'
              call usage('!! file format (mag typ)')
          endif
          call wellscoppersmith(wid,len,mag,typ)
          mom = (10.0d0**(1.5d0*(mag+10.7d0)))/1.0d7 ! N*m
          area = wid*len*1.0d6 ! m^3
          slip = mom/(area*shrmod)
          if (p.eq.0) then
              write(12,9999) len,wid,slip
          else
              write(stout,9999) len,wid,slip
          endif
          goto 102
  103 continue
 9998 format(2A12,A14)
 9999 format(2F12.3,E14.4)
      stop

C----
C User enters point source parameters
C----
  101 continue
      print *,'How many faults?'
      read *,nflt
      do 104 i = 1,nflt
          write(stout,8991) i
          read(stin,*), mag
  106     write(stout,8992) i
          read(stin,*) typc
          if (typc.eq.'ss') then
              typ = 1
          elseif (typc.eq.'no') then
              typ = 3
          elseif (typc.eq.'th'.or.typc.eq.'rv') then
              typ = 2
          else
              write(*,*) '!! Error: fault type incorrect; ss, no, th'
              goto 106
          endif
          call wellscoppersmith(wid,len,mag,typ)
          mom = (10.0d0**(1.5d0*(mag+10.7d0)))/1.0d7 ! N*m
          area = wid*len*1.0d6 ! m^3
          slip = mom/(area*shrmod)
          ans(i,1) = len
          ans(i,2) = wid
          ans(i,3) = slip
  104 continue
      write(sterr,9998) 'length(km)','width(km)','slip(m)'
      do 105 i = 1,nflt
          write(stout,9999) ans(i,1),ans(i,2),ans(i,3)
  105 continue
 8991 format('Enter fault #',I2,' magnitude')
 8992 format('Enter fault #',I2,' type ("ss", "no", "th")')

      END

C======================================================================C

      SUBROUTINE wellscoppersmith(wid,len,mag,typ)
      IMPLICIT none
      INTEGER typ
      REAL*8 wid,len,mag

      if (typ.eq.1) then
          len = 10.0d0**(-2.57d0+0.62d0*mag)
          wid = 10.0d0**(-0.76d0+0.27d0*mag)
      elseif (typ.eq.2) then
          len = 10.0d0**(-2.42d0+0.58d0*mag)
          wid = 10.0d0**(-1.61d0+0.41d0*mag)
      elseif (typ.eq.3) then
          len = 10.0d0**(-1.88d0+0.50d0*mag)
          wid = 10.0d0**(-1.14d0+0.35d0*mag)
      else
          len = 10.0d0**(-2.44d0+0.59d0*mag)
          wid = 10.0d0**(-1.01d0+0.32d0*mag)
      endif

      RETURN
      END

C----------------------------------------------------------------------C
      
      SUBROUTINE gcmdln(ifile,ofile,shrmod,u,p)
      IMPLICIT NONE
      CHARACTER*30 ifile,ofile,tag
      INTEGER i,narg,u,p
      REAL*8 shrmod
      ifile = 'none'
      ofile = 'none'
      u = 0
      p = 0
      shrmod = -1.0d0
      narg = iargc()
      if (narg.eq.0) call usage('!! Error: no command line arguments')
      i = 0
  901 i = i + 1
      if (i.gt.narg) goto 902
          call getarg(i,tag)
          if (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:4).eq.'-mod') then
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F10.0)') shrmod
          elseif (tag(1:2).eq.'-u') then
              u = 1
          elseif (tag(1:2).eq.'-p') then
              p = 1
          elseif (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
              call usage(' ')
          endif
          goto 901
  902 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage(str)
      IMPLICIT none
      INTEGER lstr
      CHARACTER str*(*)
      if (str.ne.' ') then
          lstr = len(str)
          write(*,*) str(1:lstr)
          write(*,*)
      endif
      write(*,*)
     1 'Usage: src2flt -f IFILE -o OFILE -mod SHR_MOD [-u] [-p] [-h]'
      write(*,*)
      write(*,*)
     1 '  -f IFILE        Input file: mag type("ss", "no", or "th")'
      write(*,*)
     1 '  -o OFILE        Output file: str_len(km) dip_wid(km) slip(m)'
      write (*,*)
     1 '  -mod SHR_MOD    Shear modulus of medium (GPa)'
      write (*,*)
     1 '  -u              Enter input point sources on command line'
      write(*,*)
     1 '  -p              Print to standard out instead of file'
      write (*,*)
     1 '  -h              Online help (this screen)'
      write (*,*)
      STOP
      END
