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
!!! To run hte progrma
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

!-------------------------------------------------------------------------------
! program logic
!   implicit none
!   real :: a,b
!   LOGICAL:: m,n
!    m = .true.
!   ! n= .false.
!   a= 5
!   b=4
! PRINT*, "not m is", m
! end program logic
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
!   pm(1,2)=2
!
!
!   p1=>x1
!   p2=>x2
!   p1 = 5
!   x1 = 10
!
!   print*, M
!   print*, associated(p1,x1)
!   nullify(p1)
!   print*, associated(p1),p1
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
program point
  implicit none
  real, target:: a=5,t
  real, pointer:: p
  p=>a
  print*, p,a
  nullify(p)
  p=>t
  t=1
  a=2
  print'(f5.2)',a,p
end program point
