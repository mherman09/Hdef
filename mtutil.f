      PROGRAM main
C----
C Utility for manipulating seismic source data
C----
      IMPLICIT none
      INTEGER mode
      CHARACTER*180 input,output
C Get command line arguments
      call gcmdln(mode,input,output)
C Take stdin, command line, or file and convert to file
      call parsein(input)
C Compute output
      call work(input,output,mode)
C Clean up if necessary
      if (input(1:5).eq.'86.tmp') then
          open(unit=86,file='86.tmp',status='old')
          close(86,status='DELETE')
      endif
      END

C----------------------------------------------------------------------C

      SUBROUTINE work(input,output,mode)
C----
C Read the input and operate on the values using the defined mode.
C----
      IMPLICIT none
      INTEGER mode
      CHARACTER*180 input,output
      CHARACTER*180 line,soln
      line = ''
      soln = ''
C Open input file and output file if defined
      open(unit=11,file=input,status='old')
      if (output(1:5).ne.'print') then
          open(unit=12,file=output,status='unknown')
      endif
C Read inputs
  101 read(11,'(A)',end=102) line
          if (mode.eq.11) then
              call sdr2sdr(line,soln)
          elseif (mode.eq.12) then
              call sdr2mij(line,soln)
          elseif (mode.eq.13) then
              call sdr2pnt(line,soln)
          elseif (mode.eq.14) then
              call usage('!! Error: cannot compute magnitude from SDR')
          elseif (mode.eq.15) then
              call usage('!! Error: cannot compute moment from SDR')
          elseif (mode.eq.16) then
              call sdr2ter(line,soln)
          elseif (mode.eq.21) then
              call mij2sdr(line,soln)
          elseif (mode.eq.22) then
              call usage('!! Silly, you already have MT elements!')
          elseif (mode.eq.23) then
              call mij2pnt(line,soln)
          elseif (mode.eq.24) then
              call mij2mag(line,soln)
          elseif (mode.eq.25) then
              call mij2mom(line,soln)
          elseif (mode.eq.26) then
              call mij2ter(line,soln)
          elseif (mode.eq.31) then
              call pnt2sdr(line,soln)
          elseif (mode.eq.32) then
              call usage('!! PNT to MIJ not operational yet')
          elseif (mode.eq.33) then
              call usage('!! Silly, you already have PNT axes!')
          elseif (mode.eq.34) then
              call pnt2mag(line,soln)
          elseif (mode.eq.35) then
              call pnt2mom(line,soln)
          elseif (mode.eq.36) then
              call pnt2ter(line,soln)
          elseif (mode.eq.41) then
              call usage('!! Error: cannot compute SDR from MAG')
          elseif (mode.eq.42) then
              call usage('!! Error: cannot compute MIJ from MAG')
          elseif (mode.eq.43) then
              call usage('!! Error: cannot compute PNT from MAG')
          elseif (mode.eq.44) then
              call usage('!! Silly, you have the magnitude!')
          elseif (mode.eq.45) then
              call mag2mom(line,soln)
          elseif (mode.eq.46) then
              call usage('!! Error: cannot compute Ternary from MAG')
          elseif (mode.eq.51) then
              call usage('!! Error: cannot compute SDR from MOM')
          elseif (mode.eq.52) then
              call usage('!! Error: cannot compute MIJ from MOM')
          elseif (mode.eq.53) then
              call usage('!! Error: cannot compute PNT from MOM')
          elseif (mode.eq.54) then
              call mom2mag(line,soln)
          elseif (mode.eq.55) then
              call usage('!! Silly, you have the moment!')
          elseif (mode.eq.56) then
              call usage('!! Error: cannot compute Ternary from MAG')
          endif
C Write outputs
          if (output(1:5).eq.'print') then
              write(*,*) trim(soln)
          else
              write(12,*) trim(soln)
          endif
          goto 101
  102 continue
      close(11)
      if (output(1:5).ne.'print') close(12)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE parsein(input)
C----
C Determine whether input is standard input, on command line, or is a
C file name.
C----
      IMPLICIT none
      CHARACTER*180 input
      INTEGER i
      LOGICAL ex
C Standard input
      if (input(1:5).eq.'stdin') then
          open(unit=11,file='86.tmp',status='unknown')
  102     read(*,'(A)',end=103) input
              write(11,*) input
              goto 102
  103     continue
          close(11,status='KEEP')
          input = '86.tmp'
C File
      else
          i = index(input,',')
          if (i.eq.0) then
              inquire(FILE=input,EXIST=ex)
              if (.not.ex) call usage('!! Error: no input file '//
     1                                        trim(input))
C Command line
          else
  101         input(i:i) = ' '
              i = index(input,',')
              if (i.gt.0) goto 101
              open(unit=11,file='86.tmp',status='unknown')
              write(11,*) input
              close(11,status='KEEP')
              input = '86.tmp'
          endif
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(mode,input,output)
      IMPLICIT none
      INTEGER i,narg
      CHARACTER*180 tag,input,output
      INTEGER mode
      LOGICAL ex
      mode = 0
      input = ''
      output = ''
      narg = iargc()
      if (narg.lt.3) call usage('!! Error: three arguments required')
C Get input type
      i = 1
      call getarg(i,tag)
      if (tag(1:4).eq.'-sdr') then
          mode = 10
      elseif (tag(1:4).eq.'-mij') then
          mode = 20
      elseif (tag(1:4).eq.'-pnt') then
          mode = 30
      elseif (tag(1:4).eq.'-mag') then
          mode = 40
      elseif (tag(1:4).eq.'-mom') then
          mode = 50
      elseif (tag(1:4).eq.'-ter') then
          mode = 60
      else
          call usage('!! Error: no option '//tag)
      endif
      i = i + 1
      call getarg(i,input)
      if (mode/10.eq.4.or.mode/10.eq.5) then
          inquire(file=input,exist=ex)
          if(.not.ex) input = trim(input)//','
      endif
C Get output
      i = i + 1
      call getarg(i,tag)
      if (tag(1:4).eq.'-sdr') then
          mode = mode + 1
      elseif (tag(1:4).eq.'-mij') then
          mode = mode + 2
      elseif (tag(1:4).eq.'-pnt') then
          mode = mode + 3
      elseif (tag(1:4).eq.'-mag') then
          mode = mode + 4
      elseif (tag(1:4).eq.'-mom') then
          mode = mode + 5
      elseif (tag(1:4).eq.'-ter') then
          mode = mode + 6
      else
          call usage('!! Error: no option '//tag)
      endif
      i = i + 1
      if (i.gt.narg) then
          output = 'print'
      else
          call getarg(i,output)
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage(str)
      IMPLICIT NONE
      INTEGER lstr
      CHARACTER str*(*)
      if (str.ne.' ') then
          write(*,*) trim(str)
          write(*,*)
      endif
      write(*,*)
     1 'Usage: mtutil -opt INPUT -opt [OFILE]'
      write(*,*)
     1 '    I/O options:'
      write(*,*)
     1 '    -sdr      Strike, dip, and rake angles'
      write(*,*)
     1 '    -mij      Moment tensor components (mrr mtt mpp mrt mrp ',
     2                                                         'mtp)'
      write(*,*)
     1 '    -pnt      P, N, T axes (x,y,z,mag)'
      write(*,*)
     1 '    -mag      Magnitude'
      write(*,*)
     1 '    -mom      Seismic moment (Nm)'
      write(*,*)
      write(*,*)
     1 '    Output only options:'
      write(*,*)
     1 '    -ternary  fth fss fno'
      write(*,*)
     1 '    INPUT format:'
      write(*,*)
     1 '        file name: read from file INPUT'
      write(*,*)
     1 '        comma delimited list: read from command line'
      write(*,*)
     1 '        "stdin": read from standard input'
      write(*,*)
      STOP
      END

