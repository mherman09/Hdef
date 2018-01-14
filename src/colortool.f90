program main
    use colormodule
    implicit none
    real :: l,c,h,a,b,a1,a2,b1,b2
    real :: red,grn,blu
    character(len=200) :: rgbstring
    real :: cmax
    integer :: i,n
    integer :: clip(3)

    ! User-defined variables
    real :: hue1,hue2,light1,light2,chroma1,chroma2
    character(len=200) :: option
    integer :: ncolors,vrb

    ! Parse command line
    call gcmdln(option,hue1,hue2,light1,light2,chroma1,chroma2,ncolors,vrb)
    if (vrb.gt.0) then
        write(0,*) 'Color mapping option: ',trim(option)
        write(0,*) 'Hue range:            ',hue1,hue2
        write(0,*) 'Lightness range:      ',light1,light2
        write(0,*) 'Chroma range:         ',chroma1,chroma2
        write(0,*) 'Number of colors:     ',ncolors
    endif

    n = ncolors
    if (option.eq.'spiral') then
        do i = 1,n
            l = light1  + (light2  - light1 )*(i-1)/(n-1)
            c = chroma1 + (chroma2 - chroma1)*(i-1)/(n-1)
            h = hue1    + (hue2    - hue1   )*(i-1)/(n-1)
            call maxchroma(l,h,cmax)
            if (c.gt.cmax) c = cmax
            call lch2lab(c,h,a,b)
            call lab2rgb(l,a,b,red,grn,blu)
            call checkrgb(red,grn,blu,clip)
            call formatrgb(red,grn,blu,rgbstring)
            print *,trim(rgbstring)
        enddo
    elseif (option.eq.'linear') then
        call lch2lab(chroma1,hue1,a1,b1)
        call lch2lab(chroma2,hue2,a2,b2)
        do i = 1,n
            l = light1 + (light2 - light1)*(i-1)/(n-1)
            a = a1     + (a2     - a1    )*(i-1)/(n-1)
            b = b1     + (b2     - b1    )*(i-1)/(n-1)
            call lab2lch(a,b,c,h)
            call maxchroma(l,h,cmax)
            if (c.gt.cmax) c = cmax
            call lch2lab(c,h,a,b)
            call lab2rgb(l,a,b,red,grn,blu)
            call checkrgb(red,grn,blu,clip)
            call formatrgb(red,grn,blu,rgbstring)
            print *,trim(rgbstring)
        enddo
    else
        call usage('!! Error: no option named "'//trim(option)//'"')
    endif
end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

subroutine formatrgb(red,grn,blu,rgbstring)
    implicit none
    real :: red,grn,blu
    character(len=200) :: tmpstring
    character(len=*) :: rgbstring
    integer :: j,jj
    write(tmpstring,1001) red,grn,blu
1001 format(F10.3,'/',F10.3,'/'F10.3)
    jj = 1
    do j = 1,len(rgbstring)
        rgbstring(j:j) = ' '
    enddo
    do j = 1,len(tmpstring)
        if (tmpstring(j:j).ne.' ') then
            rgbstring(jj:jj) = tmpstring(j:j)
            jj = jj + 1
        endif
    enddo
return
end

subroutine checkrgb(r,g,b,clip)
! RGB colors only make up ~50% of the visible lightness-chroma-hue spectrum
! Check that RGB values are in range 0-255 and return whether they needed to be clipped
    implicit none
    real :: r,g,b
    integer :: clip(3)
    clip = 0
    if (r.lt.0.0) then
        r = 0.0
        clip(1) = -1
    elseif (r.gt.255.0) then
        r = 255.0
        clip(1) = 1
    endif
    if (g.lt.0.0) then
        g = 0.0
        clip(2) = -1
    elseif (g.gt.255.0) then
        g = 255.0
        clip(2) = 1
    endif
    if (b.lt.0.0) then
        b = 0.0
        clip(3) = -1
    elseif (b.gt.255.0) then
        b = 255.0
        clip(3) = 1
    endif
return
end

subroutine maxchroma(lin,hin,cmax)
! Look for largest possible chroma at lightness-hue coordinate
    implicit none
    real :: lin,hin,cmax,c,c1,c2,r,g,b
    integer :: clip(3)
    c1 = 0.0
    c2 = 128.0
    c = c1
    do while (c.le.c2)
        call lch2rgb(lin,c,hin,r,g,b)
        call checkrgb(r,g,b,clip)
        if (clip(1).eq.0.and.clip(2).eq.0.and.clip(3).eq.0) then
            cmax = c
        else
            !write(0,*) 'Using max chroma= ',cmax
            return
        endif
        c = c + 0.2
    enddo
return
end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Color space conversion subroutines
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

subroutine lch2lab(ciec,cieh,ciea,cieb)
! Lightness-chroma-hue to lightness-a-b
    implicit none
    real :: ciec,cieh,ciea,cieb
    real, parameter :: pi = 4.0*atan(1.0)
    ciea = ciec*cos(cieh*pi/180.0)
    cieb = ciec*sin(cieh*pi/180.0)
return
end

subroutine lab2lch(ciea,cieb,ciec,cieh)
! Lightness-a-b to lightness-chroma-hue
    implicit none
    real :: ciea,cieb,ciec,cieh
    real, parameter :: pi = 4.0*atan(1.0)
    cieh = atan2(cieb,ciea)
    if (cieh.gt.0) then
        cieh = cieh/pi*180.0
    else
        cieh = 360.0 - abs(cieh)/pi*180.0
    endif
    ciec = sqrt(ciea*ciea+cieb*cieb)
return
end

subroutine lab2xyz (ciel,ciea,cieb,x,y,z)
! Lightness-a-b to x-y-z tristimulus
    use colormodule
    implicit none
    real :: ciel,ciea,cieb,x,y,z
    y = (ciel+16.0)/116.0
    x = ciea/500.0 + y
    z = y - cieb/200.0
    if (x**3.gt.0.008856) then
        x = x**3
    else
        x = (x-16.0/116.0)/7.787
    endif
    if (y**3.gt.0.008856) then
        y = y**3
    else
        y = (y-16.0/116.0)/7.787
    endif
    if (z**3.gt.0.008856) then
        z = z**3
    else
        z = (z-16.0/116.0)/7.787
    endif
    x = refx*x
    y = refy*y
    z = refz*z
return
end

subroutine xyz2lab(xin,yin,zin,ciel,ciea,cieb)
! X-y-z tristimulus to lightness-a-b
    use colormodule
    implicit none
    real :: xin,yin,zin,x,y,z,ciel,ciea,cieb
    x = xin/refx
    y = yin/refy
    z = zin/refz
    if (x.gt.0.008856) then
        x = x**(1/3)
    else
        x = 7.787*x + 16.0/116.0
    endif
    if (y.gt.0.008856) then
        y = y**(1/3)
    else
        y = 7.787*y + 16.0/116.0
    endif
    if (z.gt.0.008856) then
        z = z**(1/3)
    else
        z = 7.787*z + 16.0/116.0
    endif
    ciel = 116.0*y - 16.0
    ciea = 500.0*(x-y)
    cieb = 200.0*(y-z)
return
end

subroutine rgb2xyz(rin,gin,bin,x,y,z)
! Red-green-blue to x-y-z tristimulus
    implicit none
    real :: rin,gin,bin,r,g,b,x,y,z
    r = rin/255.0
    g = gin/255.0
    b = bin/255.0
    if (r.gt.0.4045) then
        r = ((r+0.055)/1.055)**2.4
    else
        r = r/12.92
    endif
    if (g.gt.0.4045) then
        g = ((g+0.055)/1.055)**2.4
    else
        g = g/12.92
    endif
    if (b.gt.0.4045) then
        b = ((b+0.055)/1.055)**2.4
    else
        b = b/12.92
    endif
    r = r*100.0
    g = g*100.0
    b = b*100.0
    x = 0.4124*r + 0.3576*g + 0.1805*b
    y = 0.2126*r + 0.7152*g + 0.0722*b
    z = 0.0193*r + 0.1192*g + 0.9505*b
return
end

subroutine xyz2rgb(xin,yin,zin,r,g,b)
! X-y-z tristimulus to red-green-blue
    implicit none
    real :: xin,yin,zin,x,y,z,r,g,b
    x = xin/100.0
    y = yin/100.0
    z = zin/100.0
    r =  3.2406*x - 1.5372*y - 0.4986*z
    g = -0.9689*x + 1.8758*y + 0.0415*z
    b =  0.0557*x - 0.2040*y + 1.0570*z
    if (r.gt.0.0031308) then
        r = 1.055*(r**(1.0/2.4)) - 0.055
    else
        r = 12.92*r
    endif
    if (g.gt.0.0031308) then
        g = 1.055*(g**(1.0/2.4)) - 0.055
    else
        g = 12.92*g
    endif
    if (b.gt.0.0031308) then
        b = 1.055*(b**(1.0/2.4)) - 0.055
    else
        b = 12.92*b
    endif
    r = r*255.0
    g = g*255.0
    b = b*255.0
return
end

subroutine lab2rgb(lin,ain,bin,r,g,b)
! Lightness-a-b to red-green-blue
    implicit none
    real :: lin,ain,bin,r,g,b
    real :: x,y,z
    call lab2xyz(lin,ain,bin,x,y,z)
    call xyz2rgb(x,y,z,r,g,b)
return
end

subroutine lch2rgb(lin,cin,hin,r,g,b)
! Lightness-chroma-hue to red-green-blue
    implicit none
    real :: lin,cin,hin,r,g,b
    real :: a1,b1
    call lch2lab(cin,hin,a1,b1)
    call lab2rgb(lin,a1,b1,r,g,b)
return
end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! User interface subroutines
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

subroutine gcmdln(option,hue1,hue2,light1,light2,chroma1,chroma2,ncolors,vrb)
    implicit none
    integer :: i,j,narg
    character(len=200) :: tag
    character(len=*) :: option
    character(len=1) :: ch1,ch2
    real :: hue1,hue2,light1,light2,chroma1,chroma2
    integer :: ncolors,vrb
    option = ''
    hue1 = 0.0
    hue2 = 360.0
    light1 = 0.0
    light2 = 100.0
    chroma1 = 40.0
    chroma2 = 40.0
    ncolors = 10
    vrb = 0
    narg = iargc()
    if (narg.eq.0) call usage('')
    i = 1
    do while (i.le.narg)
        call getarg(i,tag)
        if (tag(1:5).eq.'-cmap') then
            i = i + 1
            call getarg(i,option)
        elseif (tag(1:4).eq.'-hue') then
            hue1 = -1.0
            hue2 = -1.0
            i = i + 1
            call getarg(i,tag)
            j = index(tag,',')
            tag(j:j) = ' '
            read(tag,*) ch1,ch2
            if (ch1.eq.'m') then
                hue1 = 0
            elseif (ch1.eq.'r') then
                hue1 = 30
            elseif (ch1.eq.'o') then
                hue1 = 50
            elseif (ch1.eq.'y') then
                hue1 = 90
            elseif (ch1.eq.'g') then
                hue1 = 150
            elseif (ch1.eq.'b') then
                hue1 = 270
            elseif (ch1.eq.'p') then
                hue1 = 330
            endif
            if (ch2.eq.'m') then
                hue2 = 0
            elseif (ch2.eq.'r') then
                hue2 = 30
            elseif (ch2.eq.'o') then
                hue2 = 50
            elseif (ch2.eq.'y') then
                hue2 = 90
            elseif (ch2.eq.'g') then
                hue2 = 150
            elseif (ch2.eq.'b') then
                hue2 = 270
            elseif (ch2.eq.'p') then
                hue2 = 330
            endif
            if (hue1.lt.-0.5.and.hue2.lt.-0.5) then
                read(tag,*) hue1,hue2
            elseif (hue1.lt.-0.5) then
                read(tag,*) hue1,ch2
            elseif (hue2.lt.-0.5) then
                read(tag,*) ch1,hue2
            endif
        elseif (tag(1:6).eq.'-light') then
            i = i + 1
            call getarg(i,tag)
            j = index(tag,',')
            tag(j:j) = ' '
            read(tag,*) light1,light2
        elseif (tag(1:7).eq.'-chroma') then
            i = i + 1
            call getarg(i,tag)
            j = index(tag,',')
            tag(j:j) = ' '
            read(tag,*) chroma1,chroma2
        elseif (tag(1:8).eq.'-ncolors') then
            i = i + 1
            call getarg(i,tag)
            read(tag,*) ncolors
        elseif (tag(1:2).eq.'-v') then
            vrb = 1
        elseif (tag(1:2).eq.'-d') then
            call usage('long')
        else
            call usage('!! Error: no option '//trim(tag))
        endif
        i = i + 1
    enddo
return
end

subroutine usage(str)
    implicit none
    character(len=*) :: str
    if (str.ne.''.and.str.ne.'long') then
        write(0,*) str
        write(0,*)
    endif
    write(0,*) 'Usage: colortool -cmap OPTION [-hue H1,H2] [-lightness L1,L2] [-chroma C1,C2]'
    write(0,*) '                 [-ncolors N] [-v verbose] [-d]'
    write(0,*)
    write(0,*) '-cmap OPTION            Generate a perceptually linear colormap'
    if (str.eq.'long') then
        write(0,*) '    spiral: rotate through hues in lightness-chroma-hue space (r->o->y->g->b->p)'
        write(0,*) '    linear: draw line through lightness-a-b space (a=chroma*cos(hue), b=chroma*sin(hue))'
        write(0,*)
    endif
    write(0,*) '-hue H1,H2              Starting and ending hues (default: 0,360)'
    if (str.eq.'long') then
        write(0,*) '    Can also define hues with the following names or first letters'
        write(0,*) '        0=360=magenta'
        write(0,*) '        30=red'
        write(0,*) '        50=orange'
        write(0,*) '        90=yellow'
        write(0,*) '        150=green'
        write(0,*) '        270=blue'
        write(0,*) '        330=purple'
        write(0,*)
    endif
    write(0,*) '-lightness L1,L2        Starting and ending hues (default: 0,100)'
    if (str.eq.'long') then
        write(0,*)
    endif
    write(0,*) '-chroma C1,C2           Starting and ending chroma (default: 40,40)'
    if (str.eq.'long') then
        write(0,*) '    Chroma is similar to saturation; its max value depends on hue and lightness'
        write(0,*) '    0=grey, max=full color'
        write(0,*) '    The max chroma will be limited by the program automatically'
        write(0,*)
    endif
    write(0,*) '-ncolors N              Number of colors (default: 10)'
    if (str.eq.'long') then
        write(0,*)
    endif
    write(0,*) '-verbose                Turn on verbose mode'
    if (str.eq.'long') then
        write(0,*)
    endif
    write(0,*) '-d                      Detailed help'
    stop
return
end



