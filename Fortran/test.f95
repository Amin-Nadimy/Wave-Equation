! PROGRAM sort
! !IMPLICIT NONE
! INTEGER, DIMENSION (:), ALLOCATABLE:: nos
! INTEGER:: i,j,n,temp
! PRINT*, "what's the length og nos?"
! READ*, n
! PRINT*, 'gime me numbers in each row'
! READ*, nos
! ALLOCATE( nos(n))
!
! DO i=1,n-1
!   DO j=i+1,n
!     IF (nos(i)>nos(j)) THEN
!       temp = nos(j)
!       nos(i)=nos(j)
!       nos(i)= temp
!     END IF
!   END DO
! END DO
!
! DO i=1,n
!   PRINT*, "index",i, "values", nos(i)
! END DO
!
! DEALLOCATE(nos)
! END PROGRAM sort

!-------------------------------------------------------------------------------
!!! To run the progrma
! Creates an executable file
! gfortran -o test test.f95

! To run the exe file
! ./test
!-------------------------------------------------------------------------------
! ! give the norm of a vector
! PROGRAM norm
!   IMPLICIT NONE
!   REAL:: x,y,z,a,mynorm
!   PRINT*, "give coordiate"
!   READ*, x,y,z
!   a = mynorm(x,y,z)
!   PRINT*, a
! END PROGRAM norm
!
! FUNCTION mynorm(xx,yy,zz) RESULT(res)
!   IMPLICIT NONE
!   REAL, INTENT(IN):: xx,yy,zz
!   REAL:: res
!   res = SQRT(xx**2 + yy**2 + zz**2)
! END FUNCTION mynorm
!----------------------------- arrays ------------------------------------------
! program test
!   implicit NONE
!   real,allocatable, dimension(:,:)::numbers
!   integer:: nID_size, ndim_size, i,j
!   nID_size = 3
!   ndim_size = 5
!   allocate(numbers(nID_size, ndim_size))
!   do i=1,min(nID_size,ndim_size)
!     numbers(i,i)=1
!   end do
! !  allocate(numbers(size))
!
! ! this prints the numbers in the matrix form
! do i=1,nID_size
!   do j=1,ndim_size
!     !print*, numbers(i,j)
!     write(*,'(f8.1, t3)', advance='no'),numbers(i,j)
!   end do
!   write(*,*)
! end do
!
! !  deallocate(numbers)
!
! end program test

!-------------------------- open, read, allocate ------------------------------
! program test
!   implicit NONE
!   integer:: i,size,j,size2
!   real, allocatable, dimension(:)::ID, x,y,z, ID2, x2,y2,z2
!   size=9
!   size2=size*2-1
!   allocate(ID(size), x(size), y(size), z(size))
!   allocate(ID2(size2), x2(size2), y2(size2), z2(size2))
!
!   print*, 'hi1'
!   ! Reading mesh data
!   open(4, file='node_coo.txt')
!     do i=1,size
!       read(4,*) ID(i), x(i), y(i), z(i)
!     end do
!   close(4)
! print*, 'hi2'
!
! ! inserting new points between initial nodes
!   j=1
!   do i=1,size
!     x2(j) = x(i)
!     y2(j) = y(i)
!     z2(j) = z(i)
!     j=j+2
!   end do
!
! ! adding node IDs
!   ID2(1:size2) =  [(i, i=1,size2)] ! the same as for loop written below that
! ! do i=1,size2
! ! ID2(i) = i
! ! end do
!
! ! calculating midpoint values
!   j=2
!   do while (j<=size2)
!     x2(j)= 0.5*(x2(j-1)+x2(j+1))
!     y2(j)= 0.5*(y2(j-1)+y2(j+1))
!     z2(j)= 0.5*(z2(j-1)+z2(j+1))
!     j=j+2
!   end do
!
! print*, 'hhhhhhhh'
! ! Copy elements of node_coo.txt into node_coo2.txt and double the size
!   open(unit=6, file='node_coo2.txt')
!     do i=1,size2
!       write(6,*) ID2(i), x2(i), y2(i), z2(i)
!     end do
!   close(6)
! print*, 'hi33'
!
!
!   deallocate(ID, x,y,z, ID2, x2,y2,z2)
! end program test

! ---------------------- Dictionary --------------------------------------------
! Program Dictionary
!   implicit none
!   ! real, dimension(:), allocatable:: b
!   real, dimension(:,:), allocatable:: a
!   integer:: i,j,size=2
!   allocate(a(6,size))
!   ! allocate(b(6))
!
!   ! b=(/1,2,3,4,7,8/)
!   do i=1,size
!     a(:,i)=1
!   end do
!   a(1,:)=4
!   a(3,:)=8
!
! ! do i=1,size(b)
! !     if (present(a(i))) then
! !       print*, a(i)
! !     end if
! ! end do
! do i=1,6
!   do j=1,size
!     write(*,'(f8.1, t3)', advance='no') a(i,j)
!   end do
!   write(*,*)
! end do
!   deallocate(a)
! end program Dictionary
! !---------------------- the same as above -----------------------------------
! program writeandread
!   implicit none
!   integer::i, a(5), b(5)
!   open(unit=2, file="writeandread")
!
!   do i=1,5
!     read(2,*) a(i), b(i)
!     a(i)=a(i)*2
!     b(i)=b(i)*3
!   end do
!   close(2)
!
!   open(unit=3, file="writeandread2")
!   do i=1,5
!     write(3,*) a(i), b(i)
!   end do
!   close(3)
!   print*, a,b
!
! end program writeandread
!-------------------------------------------------------------------------------
! ! FINDS THE LARGEST
! PROGRAM X
!   IMPLICIT NONE
!   REAL:: try1, try2, try3,maxim, largest
!   PRINT*, "give 3 numers"
!   READ*, try1, try2,try3
!
!   maxim = largest(try1,try2,try3)
!   PRINT '(1f10.2)', maxim
! END PROGRAM X
!
! FUNCTION largest(a,b,c)
!   IMPLICIT NONE
!   REAL:: a,b,c, res,largest
!   IF (a>b) THEN
!     res = a
!   ELSE
!     res=b
!   END IF
!   IF (c>res) THEN
!     res = c
!   END IF
!   END FUNCTION largest
!
!!!!!!!!!!! error cannot interpret negatives well
!-------------------------------------------------------------------------------
! ! Finding are of a triangle
! real function my_tri(aa,bb,cc) result(res)
!   implicit none
!   real, intent(in):: aa,bb,cc
!   real:: s
!   s=(aa+bb+cc)/2
!   res = sqrt(s*(s-aa)*(s-bb)*(s-cc))
! end function my_tri
!
! program tri_area
! implicit none
! real:: a,b,c,s,area, my_tri
! a=4
! b=4
! c=5
! area = my_tri(a,b,c)
! print*, "area of the triangle is", area
! end program tri_area

! !-------------------------------------------------------------------------------
! !    Copyright (C) 2006 Imperial College London and others.
! !
! !    Please see the AUTHORS file in the main source directory for a full list
! !    of copyright holders.
! !
! !    Prof. C Pain
! !    Applied Modelling and Computation Group
! !    Department of Earth Science and Engineering
! !    Imperial College London
! !
! !    amcgsoftware@imperial.ac.uk
! !
! !    This library is free software; you can redistribute it and/or
! !    modify it under the terms of the GNU Lesser General Public
! !    License as published by the Free Software Foundation,
! !    version 2.1 of the License.
! !
! !    This library is distributed in the hope that it will be useful,
! !    but WITHOUT ANY WARRANTY; without even the implied warranty of
! !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! !    Lesser General Public License for more details.
! !
! !    You should have received a copy of the GNU Lesser General Public
! !    License along with this library; if not, write to the Free Software
! !    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
! !    USA
!
! #include "fdebug.h"
!
! module read_gmsh
!   ! This module reads GMSH files and results in a vector field of
!   ! positions.
!
!   use fldebug
!   use global_parameters, only : OPTION_PATH_LEN
!   use futils
!   use quadrature
!   use elements
!   use spud
!   use parallel_tools
!   use fields
!   use state_module
!   use gmsh_common
!
!   implicit none
!
!   private
!
!   interface read_gmsh_file
!      module procedure read_gmsh_simple
!   end interface
!
!   public :: read_gmsh_file
!
!   integer, parameter:: GMSH_LINE=1, GMSH_TRIANGLE=2, GMSH_QUAD=3, GMSH_TET=4, GMSH_HEX=5, GMSH_NODE=15
!
! contains
!
!   ! -----------------------------------------------------------------
!   ! The main function for reading GMSH files
!
!   function read_gmsh_simple( filename, quad_degree, &
!          quad_ngi, quad_family, mdim ) &
!        result (field)
!     !!< Read a GMSH file into a coordinate field.
!     !!< In parallel the filename must *not* include the process number.
!
!     character(len=*), intent(in) :: filename
!     !! The degree of the quadrature.
!     integer, intent(in), optional, target :: quad_degree
!     !! The degree of the quadrature.
!     integer, intent(in), optional, target :: quad_ngi
!     !! What quadrature family to use
!     integer, intent(in), optional :: quad_family
!     !! Dimension of mesh
!     integer, intent(in), optional :: mdim
!     !! result: a coordinate field
!     type(vector_field) :: field
!
!     type(quadrature_type):: quad
!     type(element_type):: shape
!     type(mesh_type):: mesh
!
!     integer :: fd
!     integer,  pointer, dimension(:) :: sndglno, boundaryIDs, faceOwner
!     integer, ALLOCATABLE, DIMENSION(:,:):: faces2, dummy_faces
!     integer, ALLOCATABLE, DIMENSION(:):: nodes2_ID
!     real, ALLOCATABLE, DIMENSION(:,:):: nodes2_X
!
!     character(len = parallel_filename_len(filename)) :: lfilename
!     integer :: loc, sloc,i,j
!     integer :: numNodes, numElements, numFaces, numFaces2
!     logical :: haveBounds, haveElementOwners, haveRegionIDs
!     integer :: dim, coordinate_dim, gdim
!     integer :: gmshFormat
!     integer :: n, d, e, f, nodeID, check, check2, check3
!
!     type(GMSHnode), pointer :: nodes(:)
!     type(GMSHelement), pointer :: elements(:), faces(:)
!
!
!     ! If running in parallel, add the process number
!     if(isparallel()) then
!        lfilename = trim(parallel_filename(filename)) // ".msh"
!     else
!        lfilename = trim(filename) // ".msh"
!     end if
!
!     fd = free_unit()
!
!     ! Open node file
!     ewrite(2, *) "Opening "//trim(lfilename)//" for reading."
!     open( unit=fd, file=trim(lfilename), err=43, action="read", &
!          access="stream", form="formatted" )
!
!     ! Read in header information, and validate
!     call read_header( fd, lfilename, gmshFormat )
!
!     ! Read in the nodes
!     call read_nodes_coords( fd, lfilename, gmshFormat, nodes )
!
!     ! Read in elements
!     call read_faces_and_elements( fd, lfilename, gmshFormat, &
!          elements, faces, dim)
!
!     call read_node_column_IDs( fd, lfilename, gmshFormat, nodes )
!
!     ! According to fluidity/bin/gmsh2triangle, Fluidity doesn't need
!     ! anything past $EndElements, so we close the file.
!     close( fd )
!
!
!     numNodes = size(nodes)
!     numFaces = size(faces)
!     numElements = size(elements)
!
!     numFaces2=numFaces*2
!
!
!
!
!
!
!     allocate(nodes2_ID(numElements*12), nodes2_x(numElements*12,3), faces2(numFaces*2,7))
!     allocate(dummy_faces(numElements*3,2))
!     forall (i=1:((size(faces2)/7)), j=1:7) faces2(i,j)=-1
!     forall (i=1:numElements*12, j=1:3) nodes2_x(i,j)=-1
!     forall (i=1:numElements*12) nodes2_ID(i)=-1
!     forall (i=1:numElements*3, j=1:2) dummy_faces(i,j)=-1
!
! open(unit=12, file='amin_gmsh.txt')
! write(12,*) '$MeshFormat'
! write(12,*) '2.2 0 8'
! write(12,*) '$EndMeshFormat'
!
! write(12,*) '$Nodes'
!     ! copying old nodes
!     do i=1,numNodes
!       nodes2_ID(i) = nodes(i)%nodeID
!       forall (j=1:3) nodes2_x(i,j) = nodes(i)%x(j)
!       ! write(12,*) nodes2_ID(i), nodes2_x(i,:)
!     end do
!     ! copying new nodes
!     do i=1,numFaces
!       nodes2_ID(i+numNodes) = i+numNodes
!       forall (j=1:3) nodes2_x(i+numNodes,j) = (nodes(faces(i)%nodeIDs(1))%x(j)+ nodes(faces(i)%nodeIDs(2))%x(j))/2.0
!       ! write(12,*) nodes2_ID(i+numNodes), nodes2_X(i,:)
!       faces2(i*2-1,6)=faces(i)%nodeIDs(1)
!       faces2(i*2-1,7)=i+numNodes
!       faces2(i*2-1,2)=faces(i)%type
!       faces2(i*2-1,3)=faces(i)%numTags
!       faces2(i*2-1,4:5)=faces(i)%tags
!
!       faces2(i*2,6)=i+numNodes
!       faces2(i*2,7)=faces(i)%nodeIDs(2)
!       faces2(i*2,2)=faces(i)%type
!       faces2(i*2,3)=faces(i)%numTags
!       faces2(i*2,4:5)=faces(i)%tags
!     end do
!     forall (i=1:(size(faces2)/7)) faces2(i,1)=i
!     numNodes=numNodes+numFaces
!     do i=1,numNodes
!       write(12,*) nodes2_ID(i), nodes2_X(i,:)
!     end do
! write(12,*) '$EndNodes'
!
!
! write(12,*) '$Elements'
! ! write(12,*) 'Old Faces'
! !      do i=1,numFaces
! !        write(12,*) faces(i)%elementID, faces(i)%type, faces(i)%numTags, faces(i)%tags, faces(i)%nodeIDs
! !      end do
!
! write(12,*) 'New Faces'
!      write(12,*) numFaces*2
!      do i=1,((numFaces*2))
!          write(12,*) faces2(i,:)
!      end do
!
!     do i=1,numFaces
!       dummy_faces(i,:)=faces(i)%nodeIDs
!     end do
!
! write(12,*) 'Now element'
!      do i=1, numElements
!        write(12,*) elements(i)%nodeIDs
!        do j=1,(size(dummy_faces)/2)
!          ! checks all combinations of faces if the midpoint is already calculated in the Node section
!          if     ((elements(i)%nodeIDs(1) == dummy_faces(j,1) .and. elements(i)%nodeIDs(2) == dummy_faces(j,2))&
!             .or. (elements(i)%nodeIDs(1) == dummy_faces(j,2) .and. elements(i)%nodeIDs(2) == dummy_faces(j,1))) then
!             check=0
!             exit
!          else
!             check=1
!          end if
!         end do
!         if (check == 1) then
!            ! write(12,*) 0.5*(nodes(elements(i)%nodeIDs(1))%x(1) + nodes(elements(i)%nodeIDs(2))%x(1)),&
!            !             0.5*(nodes(elements(i)%nodeIDs(1))%x(2) + nodes(elements(i)%nodeIDs(2))%x(2)),&
!            !             0.5*(nodes(elements(i)%nodeIDs(1))%x(3) + nodes(elements(i)%nodeIDs(2))%x(3))
!            dummy_faces(numFaces+1,1)=elements(i)%nodeIDs(1)
!            dummy_faces(numFaces+1,2)=elements(i)%nodeIDs(2)
!            numFaces=numFaces+1
!
!            nodes2_X(numNodes+1,1) = 0.5*(nodes(elements(i)%nodeIDs(1))%x(1) + nodes(elements(i)%nodeIDs(2))%x(1))
!            nodes2_X(numNodes+1,2) = 0.5*(nodes(elements(i)%nodeIDs(1))%x(2) + nodes(elements(i)%nodeIDs(2))%x(2))
!            nodes2_X(numNodes+1,3) = 0.5*(nodes(elements(i)%nodeIDs(1))%x(3) + nodes(elements(i)%nodeIDs(2))%x(3))
!            nodes2_ID(numNodes+1) = numNodes+1
!            numNodes = numNodes+1
!          end if
!
!         do j=1,(numFaces)
!           if    ((elements(i)%nodeIDs(1) == dummy_faces(j,1) .and. elements(i)%nodeIDs(3) == dummy_faces(j,2))&
!             .or. (elements(i)%nodeIDs(1) == dummy_faces(j,2) .and. elements(i)%nodeIDs(3) == dummy_faces(j,1))) then
!               check2 =0
!               exit
!           else
!               check2=1
!           end if
!         end do
!         if (check2==1) then
!             dummy_faces(numFaces+1,1)=elements(i)%nodeIDs(1)
!             dummy_faces(numFaces+1,2)=elements(i)%nodeIDs(3)
!             numFaces=numFaces+1
!
!             nodes2_X(numNodes+1,1) = 0.5*(nodes(elements(i)%nodeIDs(1))%x(1) + nodes(elements(i)%nodeIDs(3))%x(1))
!             nodes2_X(numNodes+1,2) = 0.5*(nodes(elements(i)%nodeIDs(1))%x(2) + nodes(elements(i)%nodeIDs(3))%x(2))
!             nodes2_X(numNodes+1,3) = 0.5*(nodes(elements(i)%nodeIDs(1))%x(3) + nodes(elements(i)%nodeIDs(3))%x(3))
!             nodes2_ID(numNodes+1) = numNodes+1
!             numNodes = numNodes+1
!         end if
!
!
!         do j=1,(numFaces)
!         if ((elements(i)%nodeIDs(2) == dummy_faces(j,1) .and. elements(i)%nodeIDs(3) == dummy_faces(j,2) )&
!                 .or. (elements(i)%nodeIDs(2) == dummy_faces(j,2) .and. elements(i)%nodeIDs(3) == dummy_faces(j,1) )) then
!           check3 =0
!           exit
!         else
!           check3=1
!         end if
!         end do
!         if (check3==1) then
!             dummy_faces(numFaces+1,1)=elements(i)%nodeIDs(2)
!             dummy_faces(numFaces+1,2)=elements(i)%nodeIDs(3)
!             numFaces=numFaces+1
!
!             nodes2_X(numNodes+1,1) = 0.5*(nodes(elements(i)%nodeIDs(2))%x(1) + nodes(elements(i)%nodeIDs(3))%x(1))
!             nodes2_X(numNodes+1,2) = 0.5*(nodes(elements(i)%nodeIDs(2))%x(2) + nodes(elements(i)%nodeIDs(3))%x(2))
!             nodes2_X(numNodes+1,3) = 0.5*(nodes(elements(i)%nodeIDs(2))%x(3) + nodes(elements(i)%nodeIDs(3))%x(3))
!             nodes2_ID(numNodes+1) = numNodes+1
!             numNodes = numNodes+1
!         end if
!      end do
! write(12,*) '$EndElements'
!
!      ! do i=1,(size(dummy_faces)/2)
!      !     write(12,*) dummy_faces(i,1), dummy_faces(i,2)
!      !   end do
!      write(12,*) '~nodes=', numNodes, numFaces2,numElements
!
!      do i=1,numNodes
!        write(12,*) nodes2_ID(i), nodes2_X(i,:)
!      end do
!
!     close(12)
!     deallocate(nodes2_ID, nodes2_x, faces2)
!     print*, 'done'
!     stop
    !---------------------- characters  --------------------------------------------
! program character
! implicit none
! character :: a
! character :: b*4 , c*6
! character, DIMENSION(2):: info*6
! info(1)= "amin"
! info(2)="nadimy"
! print*, info
!
! end program character
!-------------------------------------------------------------------------------
! program modi
!   implicit none
!   real :: nominator, denomin, reminder
!
! nominator = 5
! denomin = 3
! reminder = mod(nominator, denomin)
! print*, reminder
!
! end program modi
!------------------------- Factorial -------------------------------------------
! program factorial
! implicit none
! real :: a, fac
! integer:: i, j, n
!
! n= 4
!
! do j=1,n
!   fac = 1
!   do i= 1,j
!     fac = fac*i
!   end do
!   print*, "factorial of", j, "is",fac
! end do
!
! end program factorial
!------------------------ subroutine -------------------------------------------
! program factorial
!   implicit none
!   integer :: n
!   real:: pro
!   n = 4
!   call fact(n,pro)
!   print*, pro
! end program factorial
!
! subroutine fact(n,pro)
!   integer :: n,i
!   real:: pro
!   pro = 0
!   do i=1,n-1
!       pro = pro +n*i
!   end do
! end subroutine fact
!-------------------------------Fibonacci series -------------------------------
! program Fionacci
! integer   :: n
! integer, dimension(:), allocatable::fib,fib_old
! print*, 1
! allocate(fib(2))
! fib(1) = 1
! fib(2) = 1
! fib_old = fib
! print*, fib
! deallocate(fib)
! do n=3,9
!   allocate( fib(n))
!
!   do i= 1,n
!     fib(i) = fib_old(i-1)+fib_old(i)
!   end do
!   fib(n) = 1
!   fib_old = fib
!
!   print*, fib
!   deallocate(fib)
! end do
!
! end program Fionacci
!---------------------------- adding digits of a number ------------------------
! program digit_addition
! implicit none
! integer:: a, sum,b,c
! c=10
! a=125
! sum = 0
! do
!   if (a==0) then
!     exit
!   else
!     b = mod(a,c)
!     !print*,b
!     sum = sum + b
!     a = a/10
!   end if
! end do
! print*,sum
! end program digit_addition
!-------------------------------------------------------------------------------
! program parameter
!   real, parameter:: pi = 3.14
!   write(*,'(f6.4)') pi
! end program parameter

! !---------------------- derived data type -------------------------------------
! program derived_data_type
!   implicit none
!
!   type student
!     integer:: grade
!     character(len=20):: name
!   end type student
!
!   type(student) :: stud1
!   type(student) :: stud2
!
!
!   stud1%grade = 20
!   stud1%name = 'amin'
!
!   print*, stud1
!
!
! end program derived_data_type
!--------------------------------------- pointer -------------------------------
! program pointers
!   implicit none
!   integer, target:: x1, x2
!   integer, pointer:: P1, P2
!   integer, dimension(3,4), target:: M
!   integer, dimension(:,:), pointer:: PM
!
!   x1 = 4
!   x2 = 7
!
!   M=8
!   PM=>M
!
!   pm(1,3)=2
!
!
!   p1=>x1
!   p2=>x2
!   p1 = 5
!   x1 = 10
!
!   print*, M
!   !print*, associated(p1,x1)
!   !nullify(p1)
!   !print*, associated(p1),p1
!
! end program pointers
!-------------------------------------- module ---------------------------------
! module para
!   implicit none
!   real :: pi = 3.14, e = 2.71
! contains
!   subroutine product()
!     print*, 'pi = ',pi
!     print*, 'e = ', e
!   end subroutine product
! end module para
!
! program test
!   use para
!   implicit none
!   real:: produc, aa,bb
!   aa = 2
!   bb = 5
!   produc = pi*bb*e
!   call product()
!   print*, produc
!
! end program test
!-------------------------- subroutine -----------------------------------------
! subroutine amin(a,b,c,d)
!   implicit none
!   real, intent(inout)::a,b,c,d
!   c=a**2 + b**2
!   d=a**2 - b**2
! end subroutine amin
!
! program amin1
!   implicit none
!   real :: x,y,z,w
!   x =2
!   y=3
!   call amin(x,y,z,w)
!   print*, x,y,z,w
! end program amin1
!--------------------------- recursive function --------------------------------
! recursive function amin(n) result(fac)
!   implicit none
!   integer, intent (in):: n
!   integer:: fac
!   select case(n)
!   case(0:1)
!     fac = 1
!     case default
!       fac = n*amin(n-1)
!   end select
! end function amin
!
! program myfact
!   implicit none
!   integer :: a=3,b, amin
!   b = amin(a)
!   print*, "factrorial of",a, "is", b
! end program myfact
!--------------------------- internal subroutine -------------------------------
! program amin
!   implicit NONE
!   real:: x,y
!   x = 5.
!   y = 3.
!   print*, x,y
!   call subamin(x,y)
!   print*, x,y
!
! contains
!   subroutine subamin(a,b)
!     real:: a,b,dummy
!     dummy = a
!     a=b
!     b=dummy
!   end subroutine subamin
! end program amin
!-------------------------------------------------------------------------------
! module operators
!   implicit none
!   !real::x,y,dummy,res
! contains
!   subroutine swap(x,y)
!     real:: x,y,dummy
!     dummy = x
!     x=y
!     y=dummy
!   end subroutine swap
!
!   function times(x) result(res)
!     real:: x,res
!     res = x*10
!   end function times
! end module operators
!
! program amin
!   use operators
!   implicit none
!   real:: a,b,c
!   a=4
!   b=2
!   call swap(a,b)
!   print*, 'a=',a,'b=',b
!   !c = times(a)
!   print*, times(a)
! end program amin
!-------------------------------------------------------------------------------
! program point
!   implicit none
!   real, target:: a=5,t
!   real, pointer:: p
!   p=>a
!   print*, p,a
!   nullify(p)
!   p=>t
!   t=1
!   a=2
!   print'(f5.2)',a,p
! end program point
!------------------------------ interface --------------------------------------
! real function power_real(a,b) result(c)
!   implicit NONE
!   real, intent(in):: a,b
!   c=a**b
! end function power_real
!
! real function power_int(a,b) result(c)
!   implicit NONE
!   integer, intent(in):: a,b
!   c=a**b
! end function power_int
!
! PROGRAM inter
!   implicit NONE
!   real:: a1=2.0, b1=3.0,c1, c2, power
!   integer:: a2=3, b2=2
!
!   interface pow
!     real function power_real(a,b) result(c)
!       implicit NONE
!       real, intent(in):: a,b
!     end function power_real
!
!     real function power_int(a,b) result(c)
!       implicit NONE
!       integer, intent(in):: a,b
!     end function power_int
!
!   end interface pow
!
!   c1=pow(a1,b1)
!   print*, c1
! END PROGRAM inter

!---------------------------- Read and write in external file ------------------
! program read_write
!   implicit NONE
!   real:: num, x,y,z
!   integer:: i
!   open(3, file='test_2')
!   ! do i=1,3
!   !   num = i*2
!   !   write(3,*) num
!   ! end do
!   read(3,*) x,y,z
!   print*, x,y,z
! print *, "done"
!
! end program read_write

!------------------------------- present ---------------------------------------
!The PRESENT function determines whether or not an optional argument is present.
! a is INTENT(IN) and must be an optional dummy argument of the procedure
! in which the PRESENT function appears.
!syntax = present(a)
! function fun1(a,b)
!   real :: fun1,a
!   real,optional :: b
!   if(present(b)) then ! don't reference b unless
!     fun1=a+b          ! it is present
!   else
!     fun1=a
!   end if
! end function
!------------------------------ forall -----------------------------------------
! it is the same as for loop but you can have nested loops in onw line
! syntax: forll(j=1:x, j=1:y)
! you can mask a condition too. meaning that an if statement can be i the same line
! e.g forll(j=1:x, j=1:y, i/=2.andj/=3)
! it menas do the operatin everywhere except when i is not 2 and j is not 3
! note: forall perform simultaneuosly on all components, but when do loop is
! used, the matrix is updated after each i or j operation, so the results
! of do loop and ifall can be different.
! note: if you follow this structure: forall (<variables and consitions>)
!                                              operations
!                                    end forall
! or no need end forall if: forall (<vabls and conditions>) <operations>
!------------------------------------------------------------------------------
!
! subroutine dg_advection_general(vec,c,rhs, totele,nloc,totele_nloc, sngi, ngi, ndim, ndim_navier, nface, max_face_list_no, nc, &
!                   got_shape_funs, n, nlx, nlxx, nlx_lxx, weight, nlx_nod,  face_sn, face_sn2, face_snlx, face_sweigh, npoly,ele_type, & ! shape functions
!                   face_ele, face_list_no,  & ! face info
!                   got_xc, xc, & ! element centres
!                   den,s, u,mu,muvol,x_all, p, & ! fields
!                   den_bc,c_bc,u_bc,mu_bc,mu_bc_on,p_bc,p_bc_on, & ! bcs
!                   n_col_c, coln_c, fin_c, cent_c, h, suf_h, suf_s, & !
!                   k_6h_p_gi, k_4h_p_gi, k_2h_p_gi ) ! pressure stabilization terms
! ! ******************************************************************************
! ! this sub form rhs = A*c -integration of s.**********************************
! ! ******************************************************************************
!       implicit none
!       integer nits_multi_pass
!       real toler, two_thirds
!       parameter(nits_multi_pass=1,toler=1.0e-10,two_thirds=2./3.)
! ! integers representing the length of arrays...
! ! totele=no of elements,nloc=no of nodes per element, totele_nloc=totele*nloc
! ! sngi=no of surface quadrature points of the faces - this is set to the max no of all faces.
! ! ngi=no of surface quadrature points of the faces.
! ! ndim=no of dimensions - including time possibly, nface=no of faces of each elemenet, nc no of fields to solve for.
!       integer, intent(in) :: totele,nloc,totele_nloc, sngi, ngi, ndim, ndim_navier, nface, max_face_list_no, nc
! ! this sub form rhs = A*vec -integration of s.
!       real, intent(in) :: c(nc,totele_nloc), rhs(nc,totele_nloc)
!       real, intent(inout) :: vec(nc,totele_nloc)
! ! shape functions....
! ! if .not.got_shape_funs then get the shape functions else assume we have them already
! ! n, nlx are the volume shape functions and their derivatives in local coordinates.
! ! weight are the weights of the quadrature points.
! ! nlx_nod are the derivatives of the local coordinates at the nods.
! ! nlx_lxx = the 3rd order local derivatives at the nodes.
! ! face info:
! ! face_ele(iface, ele) = given the face no iface and element no return the element next to
! ! the surface or if negative return the negative of the surface element number between element ele and face iface.
! ! face_list_no(iface, ele) returns the possible origantation number which defines the numbering
! ! of the non-zeros of the nabouting element.
!       logical, intent(inout) :: got_shape_funs
!       real, intent(inout) :: n(ngi,nloc), nlx(ngi,ndim,nloc), nlxx(ngi,nloc), nlx_lxx(ngi,ndim,nloc), weight(ngi)
!       real, intent(inout) :: nlx_nod(nloc,ndim,nloc)
!       real, intent(inout) :: face_sn(sngi,nloc,nface), face_sn2(sngi,nloc,max_face_list_no), face_snlx(sngi,ndim,nloc,nface), face_sweigh(sngi,nface)
! ! npoly=order of polynomial in Cartesian space; ele_type=type of element including order of poly.
!       integer, intent(in) :: npoly,ele_type
!       integer, intent(in) :: face_ele( nface, totele), face_list_no( nface, totele) ! face info.
! ! centre of each element infor: got_xc=got the centres else calculate them here, xc(:,ele) is the
! ! centre of the element ele.
!       logical, intent(inout) :: got_xc
!       real, intent(inout) :: xc(ndim,totele)
! ! nodal variables...
! ! den=density used for each field ic; c=solution value for field ic; s=source of field ic.
! ! u=advection velocity; mu=diffusivity at the nodes; x_all=coordinates of the nodes; p=pressure
!       real, intent(in) :: den(nc,totele_nloc),s(nc,totele_nloc)
!       real, intent(in) :: u(ndim,totele_nloc),mu(nc,totele_nloc),muvol(totele_nloc),x_all(ndim,totele_nloc), p(totele_nloc)
!       real, intent(in) :: den_bc(nc,totele_nloc),c_bc(nc,totele_nloc),u_bc(ndim,totele_nloc),mu_bc(nc,totele_nloc),mu_bc_on(nc,totele_nloc)
!       real, intent(in) :: p_bc(totele_nloc),p_bc_on(totele_nloc)
! ! h absorption and surface source surf_s also surface absorption suf_h:
! ! pointers representing the sparse structure of h: n_col_c=no of potentially nonz-zeros in the nc*nc coupling,
! ! coln_c(n_col_c)=colns of this coupling, fin_c(ic) start of the ic row of h matrix, cent_c(ic)=points to the diagonal of row ic.
!       integer, intent(in) :: n_col_c, coln_c(n_col_c), fin_c(nc+1), cent_c(nc)
!       real, intent(in) :: h(n_col_c,totele_nloc), suf_h(nc,totele_nloc), suf_s(nc,totele_nloc)
! ! pressure stabilization terms(6th,4th, and 2nd order stabilization)...
!       real, intent(inout) :: k_6h_p_gi(ngi,ndim_navier,totele), k_4h_p_gi(ngi,ndim_navier,totele), k_2h_p_gi(ngi,ndim_navier,totele)
! ! local variables...
!       integer ele,iloc,jloc,ic,idim,jdim,count_c,jc,gi, iface,ele2,s_list_no, s_gi
!       logical navier_stokes
!       real sarea, norm_toler
! ! allocatable...
!       real, allocatable :: sn(:,:),sn2(:,:),snlx(:,:,:),sweigh(:)
!       real, allocatable :: den_loc(:,:), c_loc(:,:), s_loc(:,:), u_loc(:,:)
!       real, allocatable :: mu_loc(:,:), x_loc(:,:)
!       real, allocatable :: nx(:,:,:),nx_lxx(:,:,:),detwei(:),inv_jac(:,:,:),nx_nod(:,:,:)  ! shape functions
!       real, allocatable :: dengi(:,:), cgi(:,:), sgi(:,:), mugi(:,:), dx_cgi(:,:,:), dx_dengi(:,:,:)
!       real, allocatable :: x_lxx_c(:,:,:), dx_mugi(:,:,:)
!       real, allocatable :: hgi(:,:), ugi(:,:), dx_ugi(:,:), dx_cnod(:,:,:), residual(:,:), dx_pgi(:,:)
!       real, allocatable :: pgi(:), dxx_cgi(:,:,:), divu(:)
! !
!       real, allocatable :: a_star(:,:,:), j_inv_a_star(:,:), k_6h(:,:)
!       real, allocatable :: max_k_6h_pres(:), k_6h_pres(:,:)
!       real, allocatable :: vec_loc(:,:), p_loc(:)
!       real, allocatable :: den_loc2(:,:), c_loc2(:,:), u_loc2(:,:)
!       real, allocatable :: mu_loc2(:,:), mu_on_loc(:,:), p_loc2(:), p_bc_on_loc(:)
!
!       real, allocatable :: xsgi(:,:), usgi(:,:), usgi2(:,:), income(:), snorm(:,:), norm(:), sdetwei(:)
!       real, allocatable :: densgi(:,:), densgi2(:,:), csgi(:,:), csgi2(:,:), musgi12(:,:), s_cont(:)
!       real, allocatable :: psgi(:), psgi2(:), x_suf(:)
!       real, allocatable :: suf_h_loc(:,:), suf_s_loc(:,:), suf_hsgi(:,:), suf_ssgi(:,:)
!
!       real, allocatable :: dx_c_nod(:,:), divu_nod(:), visc_calc_nx_nod(:,:,:)
!       real, allocatable :: ident_dim(:,:), visc_calc_nx_nod_gi(:,:)
!       real, allocatable :: muvol_loc(:), muvolgi(:), loc_vec(:,:), x_lxx_c_norm22(:,:)
! ! for integration...
!       allocate(sn(sngi,nloc),sn2(sngi,nloc),snlx(sngi,ndim,nloc),sweigh(sngi) )
!       allocate(den_loc(nc,nloc), c_loc(nc,nloc), s_loc(nc,nloc), u_loc(ndim,nloc) )
!       allocate(mu_loc(nc,nloc), x_loc(ndim,nloc) )
!       allocate(nx(ngi,ndim,nloc),nx_lxx(ngi,ndim,nloc),detwei(ngi),inv_jac(ngi,ndim,ndim),nx_nod(nloc,ndim,nloc) ) ! shape functions
!       allocate(dengi(ngi,nc), cgi(ngi,nc), sgi(ngi,nc), mugi(ngi,nc), dx_cgi(ngi,ndim,nc), dx_dengi(ngi,ndim,nc) )
!       allocate(x_lxx_c(ngi,ndim,nc), dx_mugi(ngi,ndim,nc) )
!       allocate(hgi(ngi,n_col_c), ugi(ngi,ndim), dx_ugi(ngi,ndim), dx_cnod(ngi,ndim,nc ), residual(ngi,nc), dx_pgi(ngi,ndim) )
!       allocate(pgi(ngi), dxx_cgi(ngi,ndim,nc), divu(ngi) )
! !
!       allocate(a_star(ngi,ndim,nc), j_inv_a_star(ngi,ndim), k_6h(ngi,nc) )
!       allocate(max_k_6h_pres(ngi), k_6h_pres(ngi,ndim_navier) )
!       allocate(vec_loc(nc,nloc), p_loc(nloc) )
!       allocate(den_loc2(nc,nloc), c_loc2(nc,nloc), u_loc2(ndim,nloc) )
!       allocate(mu_loc2(nc,nloc), mu_on_loc(nc,nloc), p_loc2(nloc), p_bc_on_loc(nloc) )
!
!       allocate(xsgi(sngi,ndim), usgi(sngi,ndim), usgi2(sngi,ndim), income(sngi), snorm(sngi,ndim), norm(ndim), sdetwei(sngi) )
!       allocate(densgi(sngi,nc), densgi2(sngi,nc), csgi(sngi,nc), csgi2(sngi,nc), musgi12(sngi,nc), s_cont(nc) )
!       allocate(psgi(sngi), psgi2(sngi), x_suf(ndim) )
!       allocate(suf_h_loc(nc,totele_nloc), suf_s_loc(nc,totele_nloc), suf_hsgi(sngi,nc), suf_ssgi(sngi,nc) )
! !
!       allocate(dx_c_nod(nloc,nloc), divu_nod(nloc), visc_calc_nx_nod(ndim,ndim,nloc), ident_dim(ndim,ndim), visc_calc_nx_nod_gi(nloc,ndim) )
!       allocate(muvol_loc(nloc), muvolgi(ngi), loc_vec(nc,nloc), x_lxx_c_norm22(ngi,nc) )
!
!
!       if(.not.got_shape_funs) then
!          call get_shape_funs_spec(n, nlx, nlx_lxx, nlxx, weight, nlx_nod, &
!              nloc, sngi, ngi, ndim, nface,max_face_list_no, face_sn, face_sn2, face_snlx, face_sweigh, &
!              npoly, ele_type )
!          got_shape_funs=.true.
!       endif
!
!       if(.not.got_xc) then
!          do ele=1,totele
!            xc(idim,ele)= sum(x_all(:,(ele-1)*nloc+1:ele*nloc))/real(nloc)
!          end do
!          got_xc=.true.
!       endif
!
!       navier_stokes = (ndim_navier>0)
!
!       vec = 0.0
!
!       do ele = 1, totele ! VOLUME integral
!           ! Calculate DevFuns%DETWEI,DevFuns%RA,NX,NY,NZ for element ELE
!            den_loc(:,:) = den(:,(ele-1)*nloc+1:ele*nloc)
!            c_loc(:,:) = c(:,(ele-1)*nloc+1:ele*nloc)
!            s_loc(:,:) = s(:,(ele-1)*nloc+1:ele*nloc)
!            u_loc(:,:) = u(:,(ele-1)*nloc+1:ele*nloc)
!            mu_loc(:,:) = mu(:,(ele-1)*nloc+1:ele*nloc)
!            x_loc(:,:) = x_all(:,(ele-1)*nloc+1:ele*nloc)
!            call j_inv_det_nlx( x_loc, n, nlx, nlx_lxx, nx, nlxx, nx_lxx, detwei, weight, inv_jac, ndim,nloc,ngi, nlx_nod, nx_nod )
!
!            dengi=0.0; cgi=0.0; sgi=0.0; mugi=0.0;
!            dx_cgi=0.0; dx_dengi=0.0; x_lxx_c=0.0; dx_mugi=0.0
!            do iloc=1,nloc
!               do ic=1,nc
!                  dengi(:,ic) = dengi(:,ic) + n(:,iloc)*den_loc(ic,iloc)
!                  cgi(:,ic) = cgi(:,ic) + n(:,iloc)*c_loc(ic,iloc)
!                  sgi(:,ic) = sgi(:,ic) + n(:,iloc)*s_loc(ic,iloc)
!                  mugi(:,ic)= mugi(:,ic)+ n(:,iloc)*mu_loc(ic,iloc)
!                  do idim=1,ndim
!                     dx_cgi(:,idim,ic) = dx_cgi(:,idim,ic) + nx(:,idim,iloc)*c_loc(ic,iloc)
!                     dx_dengi(:,idim,ic) = dx_dengi(:,idim,ic) + nx(:,idim,iloc)*den_loc(ic,iloc)
!                     x_lxx_c(:,idim,ic) = x_lxx_c(:,idim,ic) + nx_lxx(:,idim,iloc)*c_loc(ic,iloc)
!                     dx_mugi(:,idim,ic)= dx_mugi(:,idim,ic) + nx(:,idim,iloc)*mu_loc(ic,iloc)
!                  end do
!               end do
!            end do
!            hgi=0.0
!            do iloc=1,nloc
!               do count_c=1,n_col_c
!                  hgi(:,count_c) = hgi(:,count_c) + n(:,iloc)*h(count_c,(ele-1)*nloc+iloc)
!               end do
!            end do
!            do ic=1,nc
!               do count_c=fin_c(ic),fin_c(ic+1)-1
!                  jc=coln_c(count_c)
!                  do iloc=1,nloc
!                     sgi(:,ic) = sgi(:,ic) - n(:,iloc)*hgi(:,count_c)*c_loc(jc,iloc)
!                  end do
!               end do
!            end do
!
!            ugi=0.0; dx_ugi=0.0
!            do iloc=1,nloc
!               do idim=1,ndim
!                  ugi(:,idim) = ugi(:,idim) + n(:,iloc)*u_loc(idim,iloc)
!                  dx_ugi(:,idim) = dx_ugi(:,idim) + nx(:,idim,iloc)*u_loc(idim,iloc)
!               end do
!            end do
! ! calculate second derivative of c: ****************************
!            dx_cnod=0.0 ! unusual variable of length dx_cnod(nloc,ndim,nc)
!            do iloc=1,nloc
!               do ic=1,nc
!                  do idim=1,ndim
!                     dx_cnod(:,idim,ic) = dx_cnod(:,idim,ic) + nx_nod(:,idim,iloc)*c_loc(ic,iloc)
!                  end do
!               end do
!            end do
!
!            dxx_cgi=0.0
!            do iloc=1,nloc
!               do ic=1,nc
!                  do idim=1,ndim
!                     dxx_cgi(:,idim,ic) = dxx_cgi(:,idim,ic) + nx(:,idim,iloc)*dx_cnod(idim,ic,iloc)
!                  end do
!               end do
!            end do
! ! calculate second derivative of c: ****************************
!
!            residual=-sgi ! source
!            do ic=1,nc
!               do idim=1,ndim
!                  residual(:,ic) = residual(:,ic)         &
!                     + mugi(:,ic) * dxx_cgi(:,idim,ic) + dx_mugi(:,idim,ic) * dx_cgi(:,idim,ic) & ! diffusion
!                     + dengi(:,ic)*cgi(:,ic)*dx_ugi(:,idim) &  ! advection
!                     + ugi(:,idim)*( cgi(:,ic)*dx_dengi(:,idim,ic) + dengi(:,ic)*dx_cgi(:,idim,ic) ) ! advection
!               end do
!            end do
!            if(navier_stokes) then ! assume first 3 fields are u,v,w
!               muvol_loc(1:nloc) = muvol( (ele-1)*nloc+1:ele*nloc )
!
!               pgi=0.0; muvolgi=0.0
!               do iloc=1,nloc
!                  pgi(:) = pgi(:) + n(:,iloc)*p((ele-1)*nloc+iloc)
!                  muvolgi(:) = muvolgi(:) + n(:,iloc)*muvol_loc(iloc)
!               end do
!
!               dx_pgi=0.0
!               do iloc=1,nloc
!                  do idim=1,ndim_navier
!                     dx_pgi(:,idim) = dx_pgi(:,idim) + nx(:,idim,iloc)* p((ele-1)*nloc+iloc) ! grad p term.
!                  end do
!               end do
!               residual(:,1:ndim_navier) = residual(:,1:ndim_navier)  + dx_pgi(:,1:ndim_navier)
!               do iloc=1,nloc
!                  do idim=1,ndim_navier
!                     vec(idim,iloc) = vec(idim,iloc) + sum( detwei(:) * nx(:,idim,iloc)*pgi(:) )
!                  end do
!               end do
! ! viscouse terms in stress form...
!               do iloc=1,nloc
!                  divu_nod(iloc) = sum( dx_c_nod(:,iloc) )
!               end do
!
!               visc_calc_nx_nod=0.0
!               do iloc=1,nloc
!               do jloc=1,nloc ! quadrature points
!                  do idim=1,ndim_navier
!                  do jdim=1,ndim_navier
!                     visc_calc_nx_nod(idim,jdim,iloc) = visc_calc_nx_nod(idim,jdim,iloc) &
!   + nx_nod(jloc,jdim,iloc)* ( - (two_thirds*mu_loc(idim,iloc)-muvol_loc(iloc)*ident_dim(idim,jdim)*divu_nod(jloc) ) &
!                        + mu_loc(idim,iloc)*( dx_c_nod(jloc,idim) + dx_c_nod(jloc,jdim) )  )
!                  end do
!                  end do
!               end do
!               end do
!
!               visc_calc_nx_nod_gi=0.0
!               do iloc=1,nloc
!                  do jdim=1,ndim
!                     do idim=1,ndim
!                        visc_calc_nx_nod_gi(iloc,idim) = visc_calc_nx_nod_gi(iloc,idim) + visc_calc_nx_nod(idim,jdim,iloc)
!                     end do
!                  end do
!               end do
!
!               residual(:,:) = residual(:,:) + visc_calc_nx_nod_gi(:,:)
! !
!               do ic=1,ndim_navier ! subtract back the diffsuion  for the velocity
!                  do idim=1,ndim
!                     residual(:,ic) = residual(:,ic)         &
!                     - mugi(:,ic) * dxx_cgi(:,idim,ic) + dx_mugi(:,idim,ic) * dx_cgi(:,idim,ic)  ! diffusion
!                  end do
!               end do
!            endif ! if(navier_stokes) then
!
! ! grad of the Lapalacian - a vector...
!            x_lxx_c = 0.0
!            do iloc=1,nloc
!               do ic=1,nc
!                  do idim=1,ndim
!                     x_lxx_c(:,idim,ic) = x_lxx_c(:,idim,ic) + nx_lxx(:,idim,iloc)*c_loc(ic,iloc)
!                  end do
!               end do
!            end do
!
!            do ic=1,nc
!               do idim=1,ndim
!                  a_star(:,idim,ic) = residual(:,ic) * x_lxx_c(:,idim,ic) / max( x_lxx_c(:,idim,ic), toler)
!               end do
!            end do
!            j_inv_a_star=0.0; x_lxx_c_norm22=0.0
!            do ic=1,nc
!               do idim=1,ndim
!                  do gi=1,ngi
!                     j_inv_a_star(gi,ic) = max( j_inv_a_star(gi,ic), abs( sum(inv_jac(gi,idim,:)*a_star(gi,:,ic)) )    )
!                     x_lxx_c_norm22(gi,ic) = x_lxx_c_norm22(gi,ic) + x_lxx_c(gi,idim,ic)**2
!                  end do
!               end do
!            end do
!            k_6h = residual**2 / max( max( sgi/(sign(cgi,1.0)*max(abs(cgi),toler)), 0.0 ), 4.0*j_inv_a_star*x_lxx_c_norm22, toler ) ! 6th order diffusion coeff
!            if(navier_stokes) then
!
!               max_k_6h_pres=0.0
!               do idim=1,ndim_navier
!                  max_k_6h_pres(:) = max( max_k_6h_pres(:), k_6h(:,idim) )
!               end do
!               k_6h_p_gi=0.0; k_4h_p_gi=0.0; k_2h_p_gi=0.0
!               do idim=1,ndim_navier
!                  k_6h_p_gi(:,idim,ele) =max_k_6h_pres(:)
!               end do
! ! viscouse terms in stress form...
!               divu=0.0
!               do idim=1,ndim_navier
!                  divu(:) = divu(:) + dx_cgi(:,idim,idim)
!               end do
!               do iloc=1,nloc
!                  do idim=1,ndim_navier
!                  do jdim=1,ndim_navier
!                     vec_loc(idim,iloc) = vec_loc(idim,iloc) &
!   + sum( detwei(:) * nx(:,jdim,iloc)* ( - (two_thirds*mugi(idim,:)-muvolgi(:))*ident_dim(idim,jdim)*divu(:) &
!                                         + mugi(idim,:)*( dx_cgi(:,idim,jdim)+dx_cgi(:,jdim,idim) )  )    )
!                  end do
!                  end do
!               end do
!            endif ! if(navier_stokes) then
!
! ! diffusion real and Petrov-Galerkin hyperdiffusion...gghgghjfkkdutytulx;goyurysjvgkkhkylhinkygytjghkyktj  gjygghg
!            do iloc=1,nloc
!               do ic=1,nc
!                  do idim=1,ndim
!                     vec_loc(ic,iloc) = vec_loc(ic,iloc) &
!                                      + sum( detwei(:) * nx(:,idim,iloc)*dengi(:,ic) * ugi(:,idim)* cgi(:,ic)   & ! advection
!                                      -      nx(:,idim,iloc)*detwei(:)*mugi(:,ic)*dx_cgi(:,idim,ic)  & ! diffusion
! !                                       - sum( nx(:,idim,iloc)*k_2h(:,ic)*detwei(:)*dx_c(:,idim,ic) ) & ! diffusion
! !                                       - sum( nlxx(:,idim,iloc)*k_4h(:,ic)*detwei(:)*lxx_c(:,idim,ic) ) & ! hyperdiffusion
!                                      -      nx_lxx(:,idim,iloc)*k_6h(:,ic)*detwei(:)*x_lxx_c(:,idim,ic) )  ! hyperdiffusion
!                  end do
!                  vec_loc(ic,iloc) = vec_loc(ic,iloc) - sum( detwei(:)*n(:,iloc)*sgi(:,ic) )
!               end do
!            end do
!
!       end do !  do ele = 1, totele ! VOLUME integral
!
!
! ! Surface integral...
!
!       do ele = 1, totele ! Surface integral
!           ! for copy local memory copying...
!          den_loc(:,:) = den(:,(ele-1)*nloc+1:ele*nloc)
!          c_loc(:,:) = c(:,(ele-1)*nloc+1:ele*nloc)
!          u_loc(:,:) = u(:,(ele-1)*nloc+1:ele*nloc)
!          mu_loc(:,:) = mu(:,(ele-1)*nloc+1:ele*nloc)
!          x_loc(:,:) = x_all(:,(ele-1)*nloc+1:ele*nloc)
!          p_loc(:) = p((ele-1)*nloc+1:ele*nloc)
!
!          loc_vec = 0.0
!          do iface = 1, nface !  Between_Elements_And_Boundary
!             ele2 = face_ele( iface, ele)
!
!
!               !Surface integral along an element
! !              if(ele>ele2) then ! perform integration along both sides...
!
!                s_list_no = face_list_no( iface, ele)
!                sn = face_sn(:,:,iface); snlx(:,:,:) = face_snlx(:,:,:,iface)
!                sn2 =face_sn2(:,:,s_list_no)
!
! ! nn for element ele
!                if(ele2<0) then ! bcs
!                   den_loc2(:,:) = den_bc(:,(ele-1)*nloc+1:ele*nloc)
!                   c_loc2(:,:) = c_bc(:,(ele-1)*nloc+1:ele*nloc)
!                   u_loc2(:,:) = u_bc(:,(ele-1)*nloc+1:ele*nloc)
!                   mu_loc2(:,:) = mu_bc(:,(ele-1)*nloc+1:ele*nloc)
!                   mu_on_loc(:,:) = mu_bc_on(:,(ele-1)*nloc+1:ele*nloc)
!                   p_loc2(:) = p_bc((ele-1)*nloc+1:ele*nloc)
!                   p_bc_on_loc(:) = p_bc_on((ele-1)*nloc+1:ele*nloc)
!                   suf_h_loc(:,:) = suf_h(:,(ele-1)*nloc+1:ele*nloc)
!                   suf_s_loc(:,:) = suf_s(:,(ele-1)*nloc+1:ele*nloc)
!
!                   x_suf=0.0
!                   do iloc=1,nloc
!                      do idim=1,ndim
!                         x_suf(idim) = x_suf(idim) + sum(sn(:,iloc))*x_loc(idim,iloc)
!                      end do
!                   end do
!                   do idim=1,ndim
!                      x_suf(idim) = x_suf(idim)/sum(sn)
!                      norm(idim) = 2.0*( x_suf(idim) - xc(idim,ele) )
!                   end do
!                else
!                   den_loc2(:,:) = den(:,(ele2-1)*nloc+1:ele2*nloc)
!                   c_loc2(:,:) = c(:,(ele2-1)*nloc+1:ele2*nloc)
!                   u_loc2(:,:) = u(:,(ele2-1)*nloc+1:ele2*nloc)
!                   mu_loc2(:,:) = mu(:,(ele2-1)*nloc+1:ele2*nloc)
!                   mu_on_loc(:,:) = 1.0
!                   p_loc2(:) = p((ele2-1)*nloc+1:ele2*nloc)
!                   p_bc_on_loc(:) = 1.0
!                   suf_h_loc(:,:) = 0.0
!                   suf_s_loc(:,:) = 0.0
!                   do idim=1,ndim
!                      norm(idim) = xc(idim,ele2) - xc(idim,ele)
!                   end do
!                endif
!                norm_toler=max(sqrt( sum( norm**2)),toler)/real(npoly)
!
!                snorm(:,ndim+1)=0.0
!                densgi = 0.0; csgi = 0.0
!                densgi2 = 0.0; csgi2 = 0.0; musgi12 = 0.0
!                suf_hsgi = 0.0; suf_ssgi = 0.0
!                do iloc=1,nloc ! use all of the nodes not just the surface nodes.
!                   do ic=1,nc
!                      densgi(:,ic) = densgi(:,ic) + sn(:,iloc)*den_loc(ic,iloc)
!                      densgi2(:,ic) = densgi2(:,ic) + sn2(:,iloc)*den_loc2(ic,iloc)
!                      csgi(:,ic) = csgi(:,ic) +     sn(:,iloc)*c_loc(ic,iloc)
!                      csgi2(:,ic) = csgi2(:,ic) +     sn2(:,iloc)*c_loc2(ic,iloc)
!                      musgi12(:,ic)= musgi12(:,ic) +0.5*( sn(:,iloc)*mu_loc(ic,iloc) + sn2(:,iloc)*mu_loc2(ic,iloc) )*mu_on_loc(ic,iloc)
!                      suf_hsgi(:,ic) = suf_hsgi(:,ic) + sn(:,iloc)*suf_h_loc(ic,iloc)
!                      suf_ssgi(:,ic) = suf_ssgi(:,ic) + sn(:,iloc)*suf_s_loc(ic,iloc)
!                   end do
!                end do
!                xsgi = 0.0; usgi = 0.0; usgi2 = 0.0
!                do iloc=1,nloc ! use all of the nodes not just the surface nodes.
!                   do idim=1,ndim
!                      xsgi(:,idim)  = xsgi(:,idim)  + sn(:,iloc)*x_loc(idim,iloc)
!                      usgi(:,idim)  = usgi(:,idim)  + sn(:,iloc)*u_loc(idim,iloc)
!                      usgi2(:,idim) = usgi2(:,idim) + sn2(:,iloc)*u_loc2(idim,iloc)
!                   end do
!                end do
! ! start to do the integration...
!                call det_snlx_all( nloc, sngi, ndim-1, ndim, x_loc, sn, snlx, sweigh, sdetwei, sarea, snorm, norm )
! ! Normals flux...
!                do s_gi=1,sngi
!                   income(s_gi)=0.5 + 0.5*sign(1.0, sum(snorm(s_gi,:)*0.5*(usgi(s_gi,:)+usgi2(s_gi,:)))  )
!                end do
! !                 sloc_vec=0.0; sloc_vec2=0.0
!                do iloc=1,nloc
!                    do ic=1,nc
!                       do idim=1,ndim
!                          s_cont(:) = snorm(:,idim)*sdetwei(:) &
!                                     *( (1.-income(:))*densgi(:,ic) * usgi(:,idim)* csgi(:,ic) &
!                                        + income(:)*densgi2(:,ic) * usgi2(:,idim)* csgi2(:,ic) &
!                                        + musgi12(:,ic)*snorm(:,idim)/norm_toler  ) &
!                                        + sdetwei(:)*(suf_hsgi(:,ic)*csgi(:,ic) - suf_ssgi(:,ic))
!                          loc_vec(ic,iloc)  = loc_vec(ic,iloc)  + sum( sn(:,iloc)*s_cont(:) )
! !                           sloc_vec2(ic,iloc) = sloc_vec2(ic,iloc) + sum( sn2(:,iloc)*s_cont(:) )
!                       end do
!                    end do
!                 end do
!                if(navier_stokes) then ! assume first 3 fields are u,v,w
!                   psgi=0.0; psgi2=0.0
!                   do iloc=1,nloc
!                      psgi(:)  = psgi(:)  + sn(:,iloc)*(p_loc(iloc) *p_bc_on_loc(iloc) + (1.-p_bc_on_loc(iloc))*p_loc2(iloc) )
!                      psgi2(:) = psgi2(:) + sn2(:,iloc)*p_loc2(iloc)
!                   end do
!                   do iloc=1,nloc
!                      do idim=1,ndim_navier
!                         s_cont(:) = snorm(:,idim)*sdetwei(:)*0.5*(psgi(:)+psgi2(:))
!                         loc_vec(idim,iloc)  = loc_vec(idim,iloc)  + sum( sn(:,iloc)*s_cont(:) )
! !                          sloc_vec2(idim,iloc) = sloc_vec2(idim,iloc) + sum( sn2(:,iloc)*s_cont(:) ))
!                      end do
!                   end do
! ! viscouse terms in stress form...
!                   divu(:) = (snorm(:,1)*(csgi(:,1)-csgi2(:,1)) + snorm(:,2)*(csgi(:,2)-csgi2(:,2)) &
!                            + snorm(:,3)*(csgi(:,3)-csgi2(:,3)) )/norm_toler
!                   do iloc=1,nloc
!                      do idim=1,ndim_navier
!                      do jdim=1,ndim_navier
!                         loc_vec(idim,iloc) = loc_vec(idim,iloc) &
!   + sum( sdetwei(:) * sn(:,iloc)*snorm(:,jdim)* ( - (two_thirds*mugi(:,idim)-muvolgi(:))*ident_dim(idim,jdim)*divu(:) &
!           + mugi(:,idim)*( snorm(:,idim)*(csgi(:,jdim)-csgi2(:,jdim))+snorm(:,jdim)*(csgi(:,idim)-csgi2(:,idim)) )/norm_toler  )  )
!                      end do
!                      end do
!                   end do
!
!                endif ! if(navier_stokes) then
! ! Put into global vector...
! ! Integrate both sides...
!          end do ! do iface = 1, nface !  Between_Elements_And_Boundary
!          vec(:,(ele-1)*nloc+1:ele*nloc)   = vec(:,(ele-1)*nloc+1:ele*nloc)   + loc_vec(:,:) ! sloc(ic,iloc)
! !           if(ele2>0) vec(:,(ele2-1)*nloc+1:ele2*nloc) = vec(:,(ele2-1)*nloc+1:ele2*nloc) - sloc_vec2(:,:)
!
! !              endif ! IF(ELE > ele2) THEN ELSE
!
!
!       end do ! do ele = 1, totele ! Surface integral
!
!       end subroutine dg_advection_general

! ================================================================================
!       call det_snlx_all( nloc, sngi, ndim-1, ndim, x_loc, sn, snlx, sweigh, sdetwei, sarea, snorm, norm )
! SUBROUTINE det_snlx_all( SNLOC, SNGI, SNDIM, ndim, XSL_ALL, SN, SNLX, SWEIGH, SDETWE, SAREA, NORMXN_ALL, NORMX_ALL )
! !       inv_jac )
!     IMPLICIT NONE
!
!     INTEGER, intent( in ) :: SNLOC, SNGI, SNDIM, ndim
!     REAL, DIMENSION( NDIM, SNLOC ), intent( in ) :: XSL_ALL
!     REAL, DIMENSION( SNGI, SNLOC ), intent( in ) :: SN
!     REAL, DIMENSION( SNGI, SNDIM, SNLOC ), intent( in ) :: SNLX
!     REAL, DIMENSION( SNGI ), intent( in ) :: SWEIGH
!     REAL, DIMENSION( SNGI ), intent( inout ) :: SDETWE
!     REAL, intent( inout ) ::  SAREA
!     REAL, DIMENSION( sngi, NDIM ), intent( inout ) :: NORMXN_ALL
!     REAL, DIMENSION( NDIM ), intent( in ) :: NORMX_ALL
! !    REAL, DIMENSION( NDIM,ndim ), intent( in ) :: inv_jac
!     ! Local variables
!     INTEGER :: GI, SL, IGLX
!     REAL :: DXDLX, DXDLY, DYDLX, DYDLY, DZDLX, DZDLY
!     REAL :: A, B, C, DETJ, RUB3, RUB4
!
!     SAREA=0.
!
!        DO GI=1,SNGI
!
!           DXDLX=0.
!           DXDLY=0.
!           DYDLX=0.
!           DYDLY=0.
!           DZDLX=0.
!           DZDLY=0.
!
!           DO SL=1,SNLOC
!              DXDLX=DXDLX + SNLX(GI,1,SL)*XSL_ALL(1,SL)
!              DXDLY=DXDLY + SNLX(GI,2,SL)*XSL_ALL(1,SL)
!              DYDLX=DYDLX + SNLX(GI,1,SL)*XSL_ALL(2,SL)
!              DYDLY=DYDLY + SNLX(GI,2,SL)*XSL_ALL(2,SL)
!              DZDLX=DZDLX + SNLX(GI,1,SL)*XSL_ALL(3,SL)
!              DZDLY=DZDLY + SNLX(GI,2,SL)*XSL_ALL(3,SL)
!           END DO
!           A = DYDLX*DZDLY - DYDLY*DZDLX
!           B = DXDLX*DZDLY - DXDLY*DZDLX
!           C = DXDLX*DYDLY - DXDLY*DYDLX
!
!           DETJ=SQRT( A**2 + B**2 + C**2)
! !          inv_jac(1,1)=DXDLX; inv_jac(1,2)=DXDLY; inv_jac(1,3)=DXDLZ
! !          inv_jac(2,1)=DyDLX; inv_jac(2,2)=DyDLY; inv_jac(2,3)=DyDLZ
! !          inv_jac(3,1)=DzDLX; inv_jac(3,2)=DzDLY; inv_jac(3,3)=DzDLZ
! !          inv_jac=inv_jac/detj
!           SDETWE(GI)=DETJ*SWEIGH(GI)
!           SAREA=SAREA+SDETWE(GI)
!
!           ! Calculate the normal at the Gauss pts...
!           ! Perform x-product. N=T1 x T2
!             CALL NORMGI(NORMXN_ALL(GI,1),NORMXN_ALL(GI,2),NORMXN_ALL(GI,3), &
!                DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,DZDLY, &
!                NORMX_ALL(1),NORMX_ALL(2),NORMX_ALL(3))
!        END DO
!
!     RETURN
!
!   END SUBROUTINE det_snlx_all
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! program arraycons
!   implicit none
!   integer :: i
!   real :: a(10) = (/(i, i=2,20, 2)/)
!   print *, a
! end program arraycons

! program test
!   implicit none
!   real :: A(3,3), B(3,3)
!   integer::i,j,n=3
!   A = reshape((/1,2,3,4,5,6,7,8,9/),(/3,3/))
!   do i=1,3
!     do j=1,3
!       print*, i,j,A(i,j)
!     end do
!     print*, i,j
!   end do
!   print*, 'now',i,j
!
! end program test
!!!!!!!!!!!!!!!!!!!!!!!!!!! plot !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!
! program plotting
!   implicit none
!   integer:: i
!   real :: dt, dx
!   real :: x(100), y(100), z(100)
!
!   open (unit=48, file='3ddata.dat')
!   do i=-40,40
!     x(i) =i*1.0/16
!     y(i) =i*1.0/16
!     z(i) = cos(x(i))+cos(y(i))
!     write(48,*) x(i), y(i), z(i)
!   end do
!   close(48)
!   call system('gnuplot -p 3ddata.plt')
! end program plotting
 !! save this in a separate file as the plot setting, with the name of 3ddata.plt
!!  set xlabel "x"
!! set ylabel "y"
!! set dgrid3d 29,19
!! set hidden3d
!! set terminal x11 0
!! set nokey
!! set grid
!! set mouse
!! set title "My 3D"
!! splot "3ddata.dat" u 1:2:3 with lines;pause mouse
!! keypress

!! thein the terminal cd to the file the use
!! gfortran 3dplot.f90 -o 3dplot
!! ./3dplot
!! https://www.youtube.com/watch?v=tScjpDeCy1s
