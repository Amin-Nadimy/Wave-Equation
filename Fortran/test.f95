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
!-------------------------------------------------------------------------------
program character
implicit none
character :: a
character :: b*4 , c*6
character, DIMENSION(2):: info*6
info(1)= "amin"
info(2)="nadimy"
print*, info

end program character
