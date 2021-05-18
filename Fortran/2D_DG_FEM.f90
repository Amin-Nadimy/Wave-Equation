! 2D DG-FEM wave equation
program wave_equation
implicit none

integer:: nsuf, ng, nt, N_e_r, N_e_c, snloc, i, j, tot_n, tot_e, tot_vol_qp, e
real:: CFL, L, dx, dy, dt
real, dimension(4,2) :: co_ordinates
real, dimension(2):: c
real,dimension(9,2) :: vol_gp
real, dimension(9):: vol_qw
integer, dimension(3):: domain_norm
integer, dimension(4):: sdot, n_hat
real,allocatable,dimension(:,:) :: M, K
real,allocatable,dimension(:) :: U, x, y, x_coo, y_coo, x_dummy, y_dummy



! allocate(M(tot_n, tot_n), K(tot_n, tot_n))
! do i=1,tot_n
!   do j=1,tot_n
!      M(i,j) = 0
!      K(i,j) = 0
!    end do
!  end do
!  forall (i=1:tot_n) U(i)=0

! Costants
! volume quadrature points
vol_gp = transpose(reshape((/-sqrt(0.6), -sqrt(0.6), 0.0, -sqrt(0.6), sqrt(0.6), -sqrt(0.6), -sqrt(0.6),&
                             0.0, 0.0, 0.0, sqrt(0.6), 0.0, -sqrt(0.6),&
                             sqrt(0.6), 0.0, sqrt(0.6), sqrt(0.6), sqrt(0.6)/),(/2,9/)))

! weightings for volume quadrature points
vol_qw = (/(real(25)/81), real(40)/81, real(25)/81, real(40)/81, real(64)/81, &
                        real(40)/81, real(25)/81, real(40)/81, real(25)/81/)

! number local surfaces for each element
nsuf = 4
! number of quadrature points on each surface
ng = 2
CFL = 0.05
!velocity in x-dir(1) and y-dir(2)
c(1) = 0.1
c(2) = 0.1
! length of the domain in each direction
L = 0.5
! number of elements in each row (r) and column (c)
N_e_r = 4
N_e_c= 3
! number of timesteps
nt = 1
! normal to the domain
domain_norm(1) = 0
domain_norm(2) = 0
domain_norm(3) = 1

allocate(x_coo((N_e_r+1)*2), y_coo((N_e_c+1)*2))
allocate(x_dummy((N_e_r+1)*2), y_dummy((N_e_c+1)*2))
allocate(x(N_e_r*2), y(N_e_c*2))

dx = L/(N_e_r)
dy = L/(N_e_c)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! in this block DG x and y coordinates are calculated
! initialisig x-coordinates
x_coo(1)=0
do i=2,(N_e_r+1)
  x_coo(i) = x_coo(i-1)+dx
end do

! initialising y-coordinates
y_coo(1)=0
do i=2,(N_e_c+1)
  y_coo(i) = y_coo(i-1)+dy
end do

! DG x & y-coordinates which have 2 extra entities. they will be deleted in x & y arrays
j=1
do i=1,N_e_r+1
  x_dummy(j) = x_coo(i)
  x_dummy(j+1) = x_coo(i)
  j=j+2
end do

j=1
do i=1,N_e_c+1
  y_dummy(j) = y_coo(i)
  y_dummy(j+1) = y_coo(i)
  j=j+2
end do

! final DG x & y coordinates
x = x_dummy(2:size(x_dummy)-1)
y = y_dummy(2:size(y_dummy)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

dt = CFL/((c(1)/dx)+(c(2)/dy))
! number of local nodes
snloc= 4
! total element
tot_e = N_e_r * N_e_c
! total nodes
tot_n = tot_e * snloc
! total volume quadrature points
tot_vol_qp = size(vol_gp)/2
! array to store dot product of normal and r (vector from center of an element &
! to a point on its edge
sdot=0
! normal to an edge
n_hat=0

! initial condition
allocate(U(tot_n))
do i=1,2
  U(N_e_r*4*i+3:N_e_r*4*i+6)=1
end do

call coordinates(12, N_e_r, dx, dy, co_ordinates)





open(unit=10, file='quad_points.txt')
! do i=1,(size(vol_gp)/2)
!   write(10,*) vol_gp(i,1), vol_gp(i,2), vol_qw(i)
! end do
write(10,*) 'next'
! print in matrix form
! do i=1,3
!   do j=1,3
!   write(10,'(f8.5)', advance='no') M(i,j)
!   end do
!   write(10,*)
! end do


write(10,*) ' x'
do i=1,N_e_r*2
write(10,*) x(i)
end do

write(10,*) ' y'
do i=1,N_e_c*2
write(10,*) y(i)
end do

write(10,*) 'element coordinates'
do i=1,4
  write(10,*) co_ordinates(i,1), co_ordinates(i,2)
end do

write(10,*) 'U='
write(10,*) U
close(10)

deallocate(x_coo,y_coo, x_dummy, y_dummy, x, y)
end program wave_equation














subroutine coordinates(e, N_e_r, dx, dy, co_ordinates)
  implicit none
  integer:: col, row, N_e_r, e
  real:: dx, dy
  real, dimension(4,2):: co_ordinates
  row = ceiling(real(e)/N_e_r)
  col = e-(N_e_r*((ceiling(real(e)/N_e_r))-1))
  co_ordinates(1,1) = (col-1)*dx
  co_ordinates(1,2) = dy*(row-1)
  co_ordinates(2,1) = dx*col
  co_ordinates(2,2) = dy*(row-1)
  co_ordinates(3,1) = dx*(col-1)
  co_ordinates(3,2) = dy*row
  co_ordinates(4,1) = dx*col
  co_ordinates(4,2) = dy*row
end subroutine coordinates
