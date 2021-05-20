! 2D DG-FEM wave equation
program wave_equation
  use cross_product
  implicit none

  integer:: nface, sngi, nt, N_e_r, N_e_c, i, j, tot_n, totele, ngi, e, s_node(4,2)
  integer :: inod, jnod, nloc, snloc, g, iloc, jloc, iface, siloc
  real, dimension(3):: domain_norm
  integer, dimension(4):: sdot, n_hat, glob_no
  real:: CFL, L, dx, dy, dt, xi, eta, det_jac
  real :: sh_func(4),jac(2,2), ddxi_sh(4), ddeta_sh(4), ddx_sh_func, ddy_sh_func
  real, dimension(4,2) :: co_ordinates
  real :: c(2), tangent(3), snormal(3), e_center(3), r(3)
  real :: vol_ngi(9,2), vol_ngw(9), s_ngi(8,2), s_ngw(8)
  real,allocatable,dimension(:,:) :: M, K
  real,allocatable,dimension(:) :: U, x, y, x_coo, y_coo, x_dummy, y_dummy

  ! Costants
  ! volume quadrature points
  vol_ngi = transpose(reshape((/-sqrt(0.6), -sqrt(0.6), 0.0, -sqrt(0.6), sqrt(0.6), -sqrt(0.6), -sqrt(0.6),&
                               0.0, 0.0, 0.0, sqrt(0.6), 0.0, -sqrt(0.6),&
                               sqrt(0.6), 0.0, sqrt(0.6), sqrt(0.6), sqrt(0.6)/),(/2,9/)))

  ! weights of volume quadrature points
  vol_ngw = (/(real(25)/81), real(40)/81, real(25)/81, real(40)/81, real(64)/81, &
                          real(40)/81, real(25)/81, real(40)/81, real(25)/81/)

  ! surface quadrature points
  s_ngi = transpose(reshape((/-3**(-0.5), -1.0, 3**(-0.5), -1.0, 1.0, -3**(-0.5),&
                               1.0, 3**(-0.5), -1.0, -3**(-0.5), -1.0, 3**(-0.5),&
                                -3**(-0.5), 1.0, 3**(-0.5), 1.0/),(/2,8/)))

  ! weights of surface quadrature points
  s_ngw = (/1,1,1,1,1,1,1,1/)
  nloc = 4  ! no of nodes in each element
  snloc = 4 ! no of nodes on each surface
  nface = 4  ! no of faces of each elemenet
  sngi = 8  ! no of surface quadrature points of the faces - this is set to the max no of all faces
  CFL = 0.05

  !velocity in x-dir(1) and y-dir(2)
  c(1) = 0.1
  c(2) = 0

  L = 0.5   ! length of the domain in each direction

  ! number of elements in each row (r) and column (c)
  N_e_r = 4
  N_e_c= 3
  nt = 1  ! number of timesteps

  ! normal to the domain
  domain_norm(1) = 0.0
  domain_norm(2) = 0.0
  domain_norm(3) = 1.0

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
  snloc= 4  ! number of local nodes
  totele = N_e_r * N_e_c    ! total element
  tot_n = totele * nloc    ! total nodes
  ngi = size(vol_ngi)/2      ! total volume quadrature points
  ! array to store dot product of normal and r (vector from center of an element &
  ! to a point on its edge
  sdot=0
  n_hat=0 ! normal to an edge

  allocate(M(tot_n, tot_n), K(tot_n, tot_n))
  do i=1,tot_n
    do j=1,tot_n
       M(i,j) = 0
       K(i,j) = 0
    end do
  end do
  ! forall (i=1:tot_n) U(i)=0

  ! initial condition
  allocate(U(tot_n))
  do i=1,2
    U(N_e_r*4*i+3:N_e_r*4*i+6)=1
  end do

  call sl_global_node(e, s_node, N_e_r)



  do e=1,totele
    ! volume integration
    call global_no(e, N_e_r, glob_no)
    call coordinates(e, N_e_r, dx, dy, co_ordinates, e_center)
    do iloc=1,nloc
      inod = glob_no(iloc)
      do jloc=1,nloc
        jnod = glob_no(jloc)
        det_jac = 0
        do g=1,ngi
          call shape_func(vol_ngi(g,1),vol_ngi(g,2), sh_func)
          call derivatives(vol_ngi(g,1),vol_ngi(g,2), ddxi_sh, ddeta_sh, e, N_e_r, dx, dy, co_ordinates, jac)
          det_jac = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)
          M(inod,jnod) = M(inod,jnod) + vol_ngw(g)*sh_func(iloc)*sh_func(jloc)*det_jac

          ! ddx_sh_func = det_jac**(-1) *(jac(2,2)*ddxi_sh(iloc) - jac(1,2)*ddeta_sh(iloc))
          ! ddy_sh_func = det_jac**(-1) *(jac(1,1)*ddeta_sh(iloc) - jac(2,1)*ddxi_sh(iloc))
          !
          ! K(inod,jnod) = K(inod,jnod) + vol_ngw(g)*dt*det_jac* (c(1)*sh_func(jloc)*ddx_sh_func + c(2)*sh_func(jloc)*ddy_sh_func)
        end do
      end do
    end do
  end do

 ! surface integration
  do e=1,totele
    do iface = 1,nface
      do siloc=1,snloc   ! use all of the nodes not just the surface nodes.
        do g=1,sngi
          call derivatives(s_ngi(g,1),s_ngi(g,2), ddxi_sh, ddeta_sh, e, N_e_r, dx, dy, co_ordinates, jac)
          call coordinates(e, N_e_r, dx, dy, co_ordinates, e_center)

          ! tangent to the iface
          tangent(1) = jac(1,1)
          tangent(2) = jac(1,2)
          tangent(3) = 0
          ! normal to the iface
          snormal = cross(tangent,domain_norm)

          ! vector from the centre of the element to a node on a boundary line
          r(1) = s_ngi(g,1) - e_center(1)
          r(2) = s_ngi(g,2) - e_center(2)
          r(3) = 0
          
        end do
      end do
    end do   ! do iface = 1, nface !  Between_Elements_And_Boundary
  end do   ! do ele = 1, totele ! Surface integral


  open(unit=10, file='quad_points.txt')
  ! print in matrix form
  ! do i=1,1
  !   do j=1,10
  !     write(10,'(f8.5)', advance='no') K(i,j)
  !   end do
  ! end do
  ! do g=1,4
  !   write(10,*) co_ordinates(g,1), co_ordinates(g,2)
  ! end do

  ! write(10,*) ' x'
  ! do i=1,N_e_r*2
  ! write(10,*) x(i)
  ! end do

  ! write(10,*) ' y'
  ! do i=1,N_e_c*2
  ! write(10,*) y(i)
  ! end do

  write(10,*) 'element coordinates'
  ! do i=1,8
  !   write(10,*) s_ngi(i,1), s_ngi(i,2)
  ! end do
  write(10,*) 'U='
  write(10,*)  'done'
  close(10)

  deallocate(U, x, y, x_coo, y_coo, x_dummy, y_dummy, M)!, K)
end program wave_equation





subroutine coordinates(e, N_e_r, dx, dy, co_ordinates, e_center)
  ! this subroutine gets no of elements in a row and ele .no
  ! and gives coordinates of 4 nodes of the ele
  implicit none
  integer:: col, row, N_e_r, e
  real:: dx, dy, e_center(3)
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

  e_center(1) = 0.25* (co_ordinates(1,1) + co_ordinates(2,1) + co_ordinates(3,1) + co_ordinates(4,1))
  e_center(2) = 0.25* (co_ordinates(1,2) + co_ordinates(2,2) + co_ordinates(3,2) + co_ordinates(4,2))
  e_center(3) = 0

end subroutine coordinates



subroutine global_no(e, N_e_r, glob_no)
  ! this subroutine gives global node numbers of an element
  implicit none
  integer :: row, e, N_e_r, glob_no(4)

  row = ceiling(real(e)/N_e_r)
  glob_no(1) = (row-1)*2*N_e_r + 2*(e-1)+1
  glob_no(2) = (row-1)*2*N_e_r + 2*(e-1)+2
  glob_no(3) = (row-1)*2*N_e_r+2*N_e_r+2*(e-1)+1
  glob_no(4) = (row-1)*2*N_e_r+2*N_e_r+2*(e-1)+2

end subroutine global_no



subroutine sl_global_node(e, s_node, N_e_r)
  ! this subroutine give global node numbers of each edge of an element
  implicit none
  integer :: e, s_node(4,2), N_e_r, row

  row = ceiling(real(e)/N_e_r)
  s_node(1,1) = (row-1)*2*N_e_r + 2*(e-1)+1
  s_node(1,2) = (row-1)*2*N_e_r + 2*(e-1)+2
  s_node(2,1) = s_node(1,2)
  s_node(2,2) = (row-1)*2*N_e_r+2*N_e_r+2*(e-1)+2
  s_node(3,1) = s_node(1,1)
  s_node(3,2) = (row-1)*2*N_e_r+2*N_e_r+2*(e-1)+1
  s_node(4,1) = s_node(3,2)
  s_node(4,2) = s_node(2,2)

end subroutine sl_global_node



subroutine shape_func(xi,eta, sh_func)
  ! this subroutine contains shape functions
  implicit none
  real :: xi, eta, sh_func(4)

  sh_func(1) = 0.25*(1-xi)*(1-eta)
  sh_func(2) = 0.25*(1+xi)*(1-eta)
  sh_func(3) = 0.25*(1-xi)*(1+eta)
  sh_func(4) = 0.25*(1+xi)*(1+eta)

end subroutine shape_func




subroutine derivatives(xi, eta, ddxi_sh, ddeta_sh, e, N_e_r, dx, dy, co_ordinates, jac)
  ! this module give derivatives of shape functions, x and y with respect to xi and eta
  implicit none
  real :: xi, eta, dx, dy, co_ordinates(4,2), jac(2,2), ddxi_sh(4), ddeta_sh(4), e_center(3)
  integer :: e, N_e_r, col, row

  ddxi_sh(1) = -0.25*(1-eta)
  ddxi_sh(2) =  0.25*(1-eta)
  ddxi_sh(3) = -0.25*(1+eta)
  ddxi_sh(4) =  0.25*(1+eta)

  ddeta_sh(1) = -0.25*(1-xi)
  ddeta_sh(2) = -0.25*(1+xi)
  ddeta_sh(3) =  0.25*(1-xi)
  ddeta_sh(4) =  0.25*(1+xi)

  call coordinates(e, N_e_r, dx, dy, co_ordinates, e_center)
  ! dx_dxi
  jac(1,1) = 0.25*((eta-1)*co_ordinates(1,1) + (1-eta)*co_ordinates(2,1) - (1+eta)*co_ordinates(3,1) + (1+eta)*co_ordinates(4,1))
  ! dy_dxi
  jac(1,2) = 0.25*((eta-1)*co_ordinates(1,2) + (1-eta)*co_ordinates(2,2) - (1+eta)*co_ordinates(3,2) + (1+eta)*co_ordinates(4,2))
  ! dx_deta
  jac(2,1) = 0.25*((xi-1)*co_ordinates(1,1) - (1+xi)*co_ordinates(2,1) + (1-xi)*co_ordinates(3,1) + (1+xi)*co_ordinates(4,1))
  ! dy_deta
  jac(2,2) = 0.25*((xi-1)*co_ordinates(1,2) - (1+xi)*co_ordinates(2,2) + (1-xi)*co_ordinates(3,2) + (1+xi)*co_ordinates(4,2))

end subroutine derivatives

module cross_product
implicit none
contains

FUNCTION cross(a, b)
  real, DIMENSION(3) :: cross
  real, DIMENSION(3), INTENT(IN) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION cross

end module cross_product
