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
!
! end program factorial
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
