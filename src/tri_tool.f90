!--------------------------------------------------------------------------------------------------!
! TRI_TOOL
!
! Utility for determining geometry and calculating characteristics of a triangular network, such as
! those produced by the Delaunay triangulation program triangle or extracted from a 2D cross-section
! through a 3D finite element model.
!
! References
! Shewchuck, J.R. (1996). Triangle: Engineering a 2D quality mesh generator and Delaunay
!     triangulator, in Applied Computational Geometry: Towards Geometric Engineering.
!--------------------------------------------------------------------------------------------------!

module tri_tool_mod

! Inputs
character(len=512) :: np_file                               ! Shewchuck triangle node file
character(len=512) :: ele_file                              ! Shewchuck triangle element file
character(len=512) :: tri_file                              ! Triangle vertices
character(len=16) :: input_mode                             ! triangle, vertices_1line, vertices_3line

! Outputs
character(len=512) :: outline_file                          ! Outline of triangle network
character(len=512) :: neighbor_file                         ! Triangles with shared sides

! Variables
integer :: nnodes                                           ! Number of unique nodes
integer :: nsegments                                        ! Number of unique line segments
integer :: ntriangles                                       ! Number of unique triangles
integer, allocatable :: nodes_in_segments(:,:)              ! Nodes making up each line segment
integer, allocatable :: nodes_in_triangles(:,:)             ! Nodes making up each triangle
integer, allocatable :: segments_in_triangles(:,:)          ! Line segments making up each triangle
integer, allocatable :: triangles_in_segments(:,:)          ! Triangles connected to each line segment
double precision, allocatable :: node_coords(:,:)           ! Node coordinates
double precision, allocatable :: tri_vert_coords(:,:,:)     ! Triangle vertex coordinates

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

subroutine deallocate_all_arrays()
implicit none
if (allocated(nodes_in_segments)) then
    deallocate(nodes_in_segments)
endif
if (allocated(nodes_in_triangles)) then
    deallocate(nodes_in_triangles)
endif
if (allocated(segments_in_triangles)) then
    deallocate(segments_in_triangles)
endif
if (allocated(triangles_in_segments)) then
    deallocate(triangles_in_segments)
endif
if (allocated(node_coords)) then
    deallocate(node_coords)
endif
if (allocated(tri_vert_coords)) then
    deallocate(tri_vert_coords)
endif
return
end subroutine deallocate_all_arrays


end module


!==================================================================================================!


program main

use io, only: verbosity, stdout

use tri_tool_mod, only: input_mode, &
                        tri_vert_coords, &
                        outline_file, &
                        neighbor_file, &
                        deallocate_all_arrays

implicit none



! Parse command line and check input values
call gcmdln()

if (verbosity.ge.1) then
    write(stdout,*) 'tri_tool: starting'
endif

call check_inputs()

if (verbosity.ge.1) then
    write(stdout,*) 'tri_tool: inputs look reasonable'
    write(stdout,*)
endif


! Clean up memory
call deallocate_all_arrays()


! Read triangle network information
if (verbosity.ge.1) then
    write(stdout,*) 'tri_tool: reading triangle network'
endif

if (input_mode.eq.'vertices_1line') then
    call read_vertices_1line()
    call build_node_list()
    call build_segment_list()
    deallocate(tri_vert_coords)
elseif (input_mode.eq.'vertices_3line') then
    call read_vertices_3line()
    call build_node_list()
    call build_segment_list()
    deallocate(tri_vert_coords)
elseif (input_mode.eq.'triangle') then
    call read_triangle()
else
    call usage('tri_tool: no input_mode named '//trim(input_mode))
endif

if (verbosity.ge.1) then
    write(stdout,*) 'tri_tool: finished reading triangle network'
    write(stdout,*)
endif


! Write specified output information
if (verbosity.ge.1) then
    write(stdout,*) 'tri_tool: writing output'
endif

if (outline_file.ne.'') then
    call write_outline_file()
endif
if (neighbor_file.ne.'') then
    call write_neighbor_file()
endif

if (verbosity.ge.1) then
    write(stdout,*) 'tri_tool: finished writing output'
    write(stdout,*)
endif


! Clean up memory
call deallocate_all_arrays()


end


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!-------------------------------------------- INPUTS ----------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine check_inputs()

use io, only: fileExists

use tri_tool_mod, only: np_file, &
                        ele_file, &
                        tri_file, &
                        input_mode

implicit none


! Make sure some input has been defined
if (input_mode.eq.'') then
    call usage('tri_tool: no input defined')
endif

! Make sure input mode and files have been defined and files exist
if (input_mode.eq.'triangle') then
    if (np_file.eq.'') then
        call usage('tri_tool: using triangle mode, no nodal point file defined')
    else
        if (.not.fileExists(np_file)) then
            call usage('tri_tool: no nodal point file found named '//trim(np_file))
        endif
    endif
    if (ele_file.eq.'') then
        call usage('tri_tool: using triangle mode, no element file defined')
    else
        if (.not.fileExists(ele_file)) then
            call usage('tri_tool: no element file found named '//trim(ele_file))
        endif
    endif

elseif (input_mode.eq.'vertices_1line') then
    if (tri_file.eq.'') then
        call usage('tri_tool: using vertices_1line mode, no triangle vertex file defined')
    else
        if (.not.fileExists(tri_file)) then
            call usage('tri_tool: no triangle vertex file found named '//trim(tri_file))
        endif
    endif

elseif (input_mode.eq.'vertices_3line') then
    if (tri_file.eq.'') then
        call usage('tri_tool: using vertices_3line mode, no triangle vertex file defined')
    else
        if (.not.fileExists(tri_file)) then
            call usage('tri_tool: no triangle vertex file found named '//trim(tri_file))
        endif
    endif

else
    call usage('tri_tool: no input mode named '//trim(input_mode))
endif


return
end subroutine check_inputs

!--------------------------------------------------------------------------------------------------!

subroutine read_vertices_1line()
!----
! Read a 2D triangular network from a file with the format:
!     x1 y1 x2 y2 x3 y3
!----

use io, only: fileExists, line_count, stderr, verbosity, stdout

use tri_tool_mod, only: tri_file, &
                        ntriangles, &
                        tri_vert_coords

implicit none

! Local variables
integer :: i, ierr
character(len=512) :: input_line


if (verbosity.ge.2) then
    write(stdout,*) 'read_vertices_1line: starting'
endif


! Open triangle network file and allocate memory to tri_vert_coords array
if (.not.fileExists(tri_file)) then
    call usage('tri_tool: no file found named "'//trim(tri_file)//'"')
endif
ntriangles = line_count(tri_file)
if (allocated(tri_vert_coords)) then
    deallocate(tri_vert_coords)
endif
allocate(tri_vert_coords(ntriangles,3,2))


! Read triangle vertex coordinates
open(unit=41,file=tri_file,status='old')
do i = 1,ntriangles
    read(41,'(A)',end=4101,err=4101,iostat=ierr) input_line
    read(input_line,*,end=4102,err=4102,iostat=ierr) tri_vert_coords(i,1,1:2), &
                                                     tri_vert_coords(i,2,1:2), &
                                                     tri_vert_coords(i,3,1:2)
enddo
4101 if (ierr.ne.0) then
    write(stderr,*) 'tri_tool: error reading tri_file at line',i,' of',ntriangles
    write(stderr,*)
    call usage('')
endif
4102 if (ierr.ne.0) then
    write(stderr,*) 'tri_tool: error parsing line',i,' of',ntriangles, &
                    ' as "x1 y1 x2 y2 x3 y3"'
    call usage('Offending line: '//trim(input_line))
endif


! Finished: close file
close(41)


if (verbosity.ge.3) then
    write(stdout,*) 'tri_vert_coords'
    write(stdout,4191) 'itri','x1','y1','x2','y2','x3','y3'
    4191 format(X,2X,A4,6(12X,A2))
    do i = 1,ntriangles
        write(stdout,4192) i,tri_vert_coords(i,1,1:2), &
                             tri_vert_coords(i,2,1:2), &
                             tri_vert_coords(i,3,1:2)
        4192 format(X,I6,6F14.6)
    enddo
endif
if (verbosity.ge.2) then
    write(stdout,*) 'read_vertices_1line: finished'
endif

return
end subroutine read_vertices_1line


!--------------------------------------------------------------------------------------------------!


subroutine read_vertices_3line()
!----
! Read a 2D triangular network from a file with the format:
!     x1 y1
!     x2 y2
!     x3 y3
!----

use io, only: fileExists, line_count, stderr, verbosity, stdout

use tri_tool_mod, only: tri_file, &
                        ntriangles, &
                        tri_vert_coords

implicit none

! Local variables
integer :: i, ierr
character(len=512) :: input_line


if (verbosity.ge.2) then
    write(stdout,*) 'read_vertices_3line: starting'
endif


! Open triangle network file and allocate memory to tri_vert_coords array
if (.not.fileExists(tri_file)) then
    call usage('tri_tool: no file found named "'//trim(tri_file)//'"')
endif
ntriangles = line_count(tri_file)/3
if (allocated(tri_vert_coords)) then
    deallocate(tri_vert_coords)
endif
allocate(tri_vert_coords(ntriangles,3,2))


! Read triangle vertex coordinates
open(unit=42,file=tri_file,status='old')
do i = 1,ntriangles
    read(42,'(A)',end=4201,err=4201,iostat=ierr) input_line
    read(input_line,*,end=4202,err=4202,iostat=ierr) tri_vert_coords(i,1,1:2)
    read(42,'(A)',end=4201,err=4201,iostat=ierr) input_line
    read(input_line,*,end=4202,err=4202,iostat=ierr) tri_vert_coords(i,2,1:2)
    read(42,'(A)',end=4201,err=4201,iostat=ierr) input_line
    read(input_line,*,end=4202,err=4202,iostat=ierr) tri_vert_coords(i,3,1:2)
enddo
4201 if (ierr.ne.0) then
    write(stderr,*) 'tri_tool: error reading tri_file at line',i,' of',ntriangles
    write(stderr,*)
    call usage('')
endif
4202 if (ierr.ne.0) then
    write(stderr,*) 'tri_tool: error parsing line',i,' of',ntriangles, &
                    ' as "x y"'
    call usage('Offending line: '//trim(input_line))
endif


! Finished: close file
close(42)


if (verbosity.ge.3) then
    write(stdout,*) 'tri_vert_coords'
    write(stdout,4291) 'itri','x1','y1','x2','y2','x3','y3'
    4291 format(X,2X,A4,6(12X,A2))
    do i = 1,ntriangles
        write(stdout,4292) i,tri_vert_coords(i,1,1:2), &
                             tri_vert_coords(i,2,1:2), &
                             tri_vert_coords(i,3,1:2)
        4292 format(X,I6,6F14.6)
    enddo
endif
if (verbosity.ge.2) then
    write(stdout,*) 'read_vertices_3line: finished'
endif

return
end subroutine read_vertices_3line


!--------------------------------------------------------------------------------------------------!


subroutine read_triangle()

use io, only: fileExists, stderr

use tri_tool_mod, only: np_file, &
                        ele_file, &
                        nnodes, &
                        ntriangles, &
                        nodes_in_triangles, &
                        node_coords

implicit none

! Local variables
integer :: inode, itri, ios
character(len=2) :: ch
character(len=512) :: input_line
logical :: isComment


! Double check that files exist
if (.not.fileExists(np_file)) then
    call usage('read_triangle: no triangle node file found named "'//trim(np_file)//'"')
endif
if (.not.fileExists(ele_file)) then
    call usage('read_triangle: no triangle element file found named "'//trim(ele_file)//'"')
endif


ios = 0


! Read node file
open(51,file=np_file,status='old')

! Read number of nodes, ignoring comment lines (start with "#")
isComment = .true.
do while (isComment)
    read(51,'(A)',end=5101,err=5101,iostat=ios) input_line
    if (index(input_line,'#').eq.1) then
        isComment = .true.
    else
        isComment = .false.
        read(input_line,*,err=5102,end=5102,iostat=ios) nnodes
    endif
enddo
5101 if (ios.ne.0) then
    call usage('read_triangle: reached end of file trying to read number of nodes')
endif
5102 if (ios.ne.0) then
    write(stderr,*) 'read_triangle: error parsing number of nodes from line'
    call usage(trim(input_line))
endif

! Allocate memory to node coordinates
if (allocated(node_coords)) then
    deallocate(node_coords)
endif
allocate(node_coords(nnodes,2))

! Read node coordinates, ignoring comment lines (start with "#")
inode = 0
do while (inode.lt.nnodes)
    read(51,'(A)',err=5103,end=5103,iostat=ios) input_line
    if (index(input_line,'#').ne.1) then
        inode = inode + 1
        read(input_line,*,err=5104,end=5104,iostat=ios) ch, node_coords(inode,1:2)
    endif
enddo
5103 if (ios.ne.0) then
    call usage('read_triangle: reached end of file while reading node coordinates')
endif
5104 if (ios.ne.0) then
    write(stderr,*) 'read_triangle: error parsing node coordinates from line'
    call usage(trim(input_line))
endif

! Finished with file
close(51)


! Read element file with nodes-in-triangles information
open(52,file=ele_file,status='old')

! Read number of triangles, ignoring comment lines (start with "#")
isComment = .true.
do while (isComment)
    read(51,'(A)',end=5105,err=5105,iostat=ios) input_line
    if (index(input_line,'#').eq.1) then
        isComment = .true.
    else
        isComment = .false.
        read(input_line,*,err=5106,end=5106,iostat=ios) ntriangles
    endif
enddo
5105 if (ios.ne.0) then
    call usage('read_triangle: reached end of file trying to read number of triangles')
endif
5106 if (ios.ne.0) then
    write(stderr,*) 'read_triangle: error parsing number of triangles from line'
    call usage(trim(input_line))
endif

! Allocate memory to node coordinates
if (allocated(nodes_in_triangles)) then
    deallocate(nodes_in_triangles)
endif
allocate(nodes_in_triangles(ntriangles,3))

! Read nodes in triangles
itri = 0
do while (itri.lt.ntriangles)
    read(51,'(A)',err=5107,end=5107,iostat=ios) input_line
    if (index(input_line,'#').ne.1) then
        itri = itri + 1
        read(input_line,*,err=5108,end=5108,iostat=ios) ch, nodes_in_triangles(itri,1:3)
    endif
enddo
5107 if (ios.ne.0) then
    call usage('read_triangle: reached end of file while reading nodes in triangles')
endif
5108 if (ios.ne.0) then
    write(stderr,*) 'read_triangle: error parsing nodes in triangles from line'
    call usage(trim(input_line))
endif

! Finished with file
close(52)


! Build line segment list and segments in triangles
call build_segment_list()


return
end subroutine read_triangle


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----------------------------------- TRIANGLE NETWORK TOOLS ---------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine build_node_list()
!----
! Build a list of unique nodes from the vertices defined in the triangle network.
! Define the variables:
!     nnodes                  Number of unique nodes
!     nodes_in_triangles      Nodes making up each triangle
!     node_coords             Node coordinates
!----

use io, only: stdout, verbosity

use tri_tool_mod, only: nnodes, &
                        ntriangles, &
                        nodes_in_triangles, &
                        node_coords, &
                        tri_vert_coords

implicit none

! Local variables
integer :: inode, itri, j
double precision :: x, y, dx, dy
double precision, allocatable :: node_coords_tmp(:,:)
logical :: nodeInList


! Allocate memory to node list arrays
if (.not.allocated(node_coords_tmp)) then
    allocate(node_coords_tmp(ntriangles*3,2))
endif
if (.not.allocated(nodes_in_triangles)) then
    allocate(nodes_in_triangles(ntriangles,3))
endif


if (verbosity.ge.2) then
    write(stdout,*) 'build_node_list: starting'
endif

! Initialize array values
nodes_in_triangles = 0
node_coords_tmp = 0.0d0


! For each node in each triangle, check whether it has been added to the node list already. If it
! is not already in the list, update the node count and add it to the list. Save the node to the
! nodes_in_triangles array.

! Start the node counter
nnodes = 0

! Loop through all triangles
do itri = 1,ntriangles

    ! Loop through all three nodes in each triangle
    do j = 1,3

        ! Initialize node as being not already in list
        nodeInList = .false.

        ! Get coordinates of current point
        x = tri_vert_coords(itri,j,1)
        y = tri_vert_coords(itri,j,2)

        ! Check whether this point coincides spatially with a previous point
        do inode = 1,nnodes
            dx = x-node_coords_tmp(inode,1)
            dy = y-node_coords_tmp(inode,2)
            if (abs(dx).lt.1.0d-6.and.abs(dy).lt.1.0d-6) then
                nodeInList = .true.
                exit
            endif
        enddo

        ! If this is a new node, update the count and add it to node list
        if (.not.nodeInList) then
            nnodes = nnodes + 1
            node_coords_tmp(nnodes,1) = x
            node_coords_tmp(nnodes,2) = y
        endif

        ! Indicate the global node number is in this triangle
        nodes_in_triangles(itri,j) = inode
    enddo
enddo


! Allocate the node_coords array to the correct size and load it
if (allocated(node_coords)) then
    deallocate(node_coords)
endif
allocate(node_coords(nnodes,2))
node_coords = node_coords_tmp(1:nnodes,:)
deallocate(node_coords_tmp)


if (verbosity.ge.3) then
    write(stdout,*) 'node_coords'
    write(stdout,4591) 'inode','x','y'
    4591 format(X,1X,A5,2(13X,A1))
    do inode = 1,nnodes
        write(stdout,4592) inode,node_coords(inode,1:2)
        4592 format(X,I6,2F14.6)
    enddo
    write(stdout,*) 'nodes_in_triangles'
    write(stdout,4593) 'itri','inode1','inode2','inode3'
    4593 format(X,2X,A4,3(4X,A6))
    do itri = 1,ntriangles
        write(stdout,4594) itri,nodes_in_triangles(itri,1:3)
        4594 format(X,I6,3I10)
    enddo
endif
if (verbosity.ge.2) then
    write(stdout,*) 'build_node_list: finished'
endif

return
end subroutine build_node_list


!--------------------------------------------------------------------------------------------------!


subroutine build_segment_list()
!----
! Build a list of unique line segments from the vertices defined in the triangle network.
! Define the variables:
!     nsegments               Number of unique line segments
!     nodes_in_segments       Nodes making up each line segment
!     segments_in_triangles   Line segments making up each triangle
!     triangles_in_segments   Triangles attached to each line segment
!----

use io, only: stdout, verbosity

use tri_tool_mod, only: nsegments, &
                        ntriangles, &
                        nodes_in_segments, &
                        nodes_in_triangles, &
                        segments_in_triangles, &
                        triangles_in_segments

implicit none

! Local variables
integer :: iseg, itri, j, j1, j2
logical :: segInList, segReverse


if (verbosity.ge.2) then
    write(stdout,*) 'build_segment_list: starting'
endif

! Allocate memory to line segment list arrays
if (allocated(nodes_in_segments)) then
    deallocate(nodes_in_segments)
endif
if (allocated(segments_in_triangles)) then
    deallocate(segments_in_triangles)
endif
allocate(nodes_in_segments(ntriangles*3,2))
allocate(segments_in_triangles(ntriangles,3))


! For each node pair in each triangle, check whether it has been added to the line segment list
! already. If it is not already in the list, update the line segment count and add it to the list.
! Save the node pair to the nodes_in_segments array and save the line segment to the
! segments_in_triangles array. Determine which triangles are attached to each line segment

! Start the line segment counter
nsegments = 0

! Loop through all triangles
do itri = 1,ntriangles

    ! Loop through all three nodes in each triangle
    do j = 1,3

        ! Initialize line segment as being not already in list and being correctly oriented
        segInList = .false.
        segReverse = .false.

        ! Get node pair
        j1 = nodes_in_triangles(itri,j)
        j2 = nodes_in_triangles(itri,mod(j,3)+1)

        ! Check whether this node pair is already defined as a line segment
        ! The orientation can be 1->2 or 2->1, so also check whether it is reversed
        do iseg = 1,nsegments
            if (nodes_in_segments(iseg,1).eq.j1.and.nodes_in_segments(iseg,2).eq.j2) then
                segInList = .true.
                segReverse = .false.
                exit
            endif
            if (nodes_in_segments(iseg,1).eq.j2.and.nodes_in_segments(iseg,2).eq.j1) then
                segInList = .true.
                segReverse = .true.
                exit
            endif
        enddo

        ! If this is a new line segment, update the count and add it to segment list
        if (.not.segInList) then
            nsegments = nsegments + 1
            nodes_in_segments(nsegments,1) = j1
            nodes_in_segments(nsegments,2) = j2
        endif

        ! Indicate the global line segment number is in this triangle
        ! Indicate that its orientation is reversed with a negative value
        segments_in_triangles(itri,j) = iseg
        if (segReverse) then
            segments_in_triangles(itri,j) = -segments_in_triangles(itri,j)
        endif
    enddo
enddo


! Find the triangles attached to each line segment
if (allocated(triangles_in_segments)) then
    deallocate(triangles_in_segments)
endif
allocate(triangles_in_segments(nsegments,3))
triangles_in_segments = 0

do itri = 1,ntriangles
    do j = 1,3

        ! Get the segment number
        iseg = abs(segments_in_triangles(itri,j))

        ! Indicate that this triangle is connected to this segment
        triangles_in_segments(iseg,1) = triangles_in_segments(iseg,1) + 1
        j1 = triangles_in_segments(iseg,1)
        if (j1.gt.2) then
            call usage('build_segment_list: you found >2 triangles connected to this segment')
        endif
        triangles_in_segments(iseg,1+j1) = itri
    enddo
enddo


if (verbosity.ge.3) then
    write(stdout,*) 'nodes_in_segments'
    write(stdout,4691) 'iseg','inode1','inode2'
    4691 format(X,2X,A4,2(4X,A6))
    do iseg = 1,nsegments
        write(stdout,4692) iseg,nodes_in_segments(iseg,1:2)
        4692 format(X,I6,2I10)
    enddo
    write(stdout,*) 'segments_in_triangles'
    write(stdout,4693) 'itri','iseg1','iseg2','iseg3'
    4693 format(X,2X,A4,3(5X,A5))
    do itri = 1,ntriangles
        write(stdout,4694) itri,segments_in_triangles(itri,1:3)
        4694 format(X,I6,3I10)
    enddo
    write(stdout,*) 'triangles_in_segments'
    write(stdout,4695) 'iseg','ntri','itri1','itri2'
    4695 format(X,2X,A4,2X,A4,2(5X,A5))
    do iseg = 1,nsegments
        if (triangles_in_segments(iseg,1).eq.1) then
            write(stdout,4696) iseg,triangles_in_segments(iseg,1:2)
            4696 format(X,I6,I6,I10)
        else
            write(stdout,4697) iseg,triangles_in_segments(iseg,1:3)
            4697 format(X,I6,I6,2I10)
        endif
    enddo
endif
if (verbosity.ge.2) then
    write(stdout,*) 'build_segment_list: finished'
endif

return
end subroutine build_segment_list


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!------------------------------------------- OUTPUTS ----------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine write_outline_file()
!----
! Print the outline of the triangle network to a file. Line segments on the outline of the network
! are only connected to one triangle.
!----

use tri_tool_mod, only: outline_file, &
                        nsegments, &
                        nodes_in_segments, &
                        triangles_in_segments, &
                        node_coords

implicit none

! Local variables
integer :: iseg, inode1, inode2


if (.not.allocated(triangles_in_segments)) then
    call usage('write_outline_file: cannot print outline with triangles_in_segments unallocated')
endif

open(unit=61,file=outline_file,status='unknown')

do iseg = 1,nsegments

    ! Segments on the outline have only one triangle attached to them
    if (triangles_in_segments(iseg,1).eq.1) then

        ! Nodes making up the line segment
        inode1 = nodes_in_segments(iseg,1)
        inode2 = nodes_in_segments(iseg,2)

        ! Write the nodal coordinates
        write(61,'(A)') '>'
        write(61,1001) node_coords(inode1,1),node_coords(inode1,2)
        write(61,1001) node_coords(inode2,1),node_coords(inode2,2)
        1001 format(1P2E14.6)

    endif
enddo

close(61)

return
end subroutine write_outline_file

!--------------------------------------------------------------------------------------------------!

subroutine write_neighbor_file()
!----
! Print the number of neighbors and neighbor index (useful for Laplacian smoothing)
!----

use tri_tool_mod, only: neighbor_file, &
                        ntriangles, &
                        segments_in_triangles, &
                        triangles_in_segments

implicit none

! Local variables
integer :: itri, iseg, j, nneighbors
integer :: neighbors(3)


if (.not.allocated(triangles_in_segments)) then
    call usage('write_neighbor_file: cannot print neighbors with triangles_in_segments unallocated')
endif
if (.not.allocated(segments_in_triangles)) then
    call usage('write_neighbor_file: cannot print neighbors with segments_in_triangles unallocated')
endif


open(unit=62,file=neighbor_file,status='unknown')

do itri = 1,ntriangles

    ! Initialize nneighbors and neighbors for the triangle
    neighbors = 0
    nneighbors = 0

    do j = 1,3

        ! The global segment number
        iseg = abs(segments_in_triangles(itri,j))

        ! The triangle itri is obviously connected to this line segment
        ! Add the other connected triangle to the neighbor list

        if (triangles_in_segments(iseg,2).eq.itri) then

            if (triangles_in_segments(iseg,3).ne.0) then
                nneighbors = nneighbors + 1
                neighbors(nneighbors) = triangles_in_segments(iseg,3)
            endif

        elseif (triangles_in_segments(iseg,3).eq.itri) then

            if (triangles_in_segments(iseg,2).ne.0) then
                nneighbors = nneighbors + 1
                neighbors(nneighbors) = triangles_in_segments(iseg,2)
            endif

        endif

    enddo

    ! Got all of the neighbors for this triangle
    ! Print them in the format for fltinv:
    !     tri_index N neighbor1 ... neighborN
    write(62,*) itri,nneighbors,neighbors(1:nneighbors)

enddo

! Finished with file
close(62)

return
end subroutine write_neighbor_file


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!------------------------------------------ COMMAND LINE ------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use io, only: verbosity

use tri_tool_mod, only: np_file, &
                        ele_file, &
                        tri_file, &
                        input_mode, &
                        outline_file, &
                        neighbor_file

implicit none

! Local variables
integer :: i, narg, ierr
character(len=512) :: tag

! Initialize control parameters
np_file = ''
ele_file = ''
tri_file = ''
input_mode = ''
outline_file = ''
neighbor_file = ''
verbosity = 0

! Number of arguments
narg = command_argument_count()
if (narg.eq.0) then
  call usage('')
endif

! Read arguments
i = 1
do while (i.le.narg)

    call get_command_argument(i,tag,status=ierr)

    if (tag.eq.'-np') then
        input_mode = 'triangle'
        i = i + 1
        call get_command_argument(i,np_file,status=ierr)
    elseif (tag.eq.'-ele') then
        input_mode = 'triangle'
        i = i + 1
        call get_command_argument(i,ele_file,status=ierr)

    elseif (tag.eq.'-tri1') then
        input_mode = 'vertices_1line'
        i = i + 1
        call get_command_argument(i,tri_file,status=ierr)

    elseif (tag.eq.'-tri3') then
        input_mode = 'vertices_3line'
        i = i + 1
        call get_command_argument(i,tri_file,status=ierr)

    elseif (tag.eq.'-outline') then
        i = i + 1
        call get_command_argument(i,outline_file,status=ierr)

    elseif (tag.eq.'-neighbor') then
        i = i + 1
        call get_command_argument(i,neighbor_file,status=ierr)

    elseif (tag.eq.'-v') then
        i = i + 1
        call get_command_argument(i,tag,status=ierr)
        read(tag,*) verbosity

    else
        call usage('tri_tool: no option '//trim(tag))
    endif

    if (ierr.ne.0) then
        call usage('tri_tool: error parsing command line arguments')
    endif

    i = i + 1

enddo

return
end subroutine gcmdln

!--------------------------------------------------------------------------------------------------!

subroutine usage(str)

use io, only: stderr

implicit none

! Arguments
character(len=*) :: str

if (str.ne.'') then
    write(stderr,*) trim(str)
    write(stderr,*)
endif

write(stderr,*) 'tri_tool -np NP_FILE -ele ELE_FILE | -tri[1|3] TRI_FILE'
write(stderr,*) '         [-outline OUTLINE_FILE] [-neighbor NEIGHBOR_FILE] [-v LVL]'
write(stderr,*)
write(stderr,*) '-np NP_FILE                 Triangle output node file'
write(stderr,*) '-ele ELE_FILE               Triangle output element file'
write(stderr,*) '-tri1 TRI_FILE              Triangles defined by vertices (on one line)'
write(stderr,*) '-tri3 TRI_FILE              Triangles defined by vertices (on three lines)'
write(stderr,*) '-outline OUTLINE_FILE       Segments that outline network (no neighbors)'
write(stderr,*) '-neighbor NEIGHBOR_FILE     # neighbors and neighbor indices for each triangle'
write(stderr,*) '-v LVL                      Turn on verbose mode'

call error_exit(1)
end subroutine usage
