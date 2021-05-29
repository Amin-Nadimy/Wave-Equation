! 2D DG-FEM wave equation
program wave_equation
  use cross_product
  implicit none

  integer:: nface, sngi, nt, N_e_r, N_e_c, i, j, tot_n, totele, ngi, e, s_node(4,2)
  integer :: inod, jnod, nloc, snloc, g, iloc, jloc, iface, siloc, sinod, ErrorFlag
  integer, dimension(4):: sdot, glob_no
  real, dimension(3):: domain_norm, n_hat
  real:: CFL, L, dx, dy, dt, xi, eta, det_jac, flux(4), K(4)
  real :: sh_func(4),jac(2,2),s_det_jac(4), ddxi_sh(4), ddeta_sh(4), ddx_sh_func, ddy_sh_func
  real, dimension(4,2) :: co_ordinates
  real :: c(2), tangent(4,3), snormal(3), e_center(3), r(3), s_dot
  real :: vol_ngi(9,2), vol_ngw(9), s_ngi(8,2), s_ngw(8)
  real,allocatable,dimension(:,:) :: M, inv_M
  real,allocatable,dimension(:) :: U, BC, x, y, x_coo, y_coo, x_dummy, y_dummy

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
  snloc = 2 ! index of i in phi_i
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
  totele = N_e_r * N_e_c    ! total element
  tot_n = totele * nloc    ! total nodes
  ngi = size(vol_ngi)/2      ! total volume quadrature points
  ! array to store dot product of normal and r (vector from center of an element &
  ! to a point on its edge
  sdot=0
  n_hat=0 ! normal to an edge

  allocate(M(tot_n, tot_n))
  allocate(inv_M(tot_n, tot_n))
  allocate(BC(4*(N_e_r+N_e_c)))
  do i=1,tot_n
    do j=1,tot_n
       M(i,j) = 0
    end do
  end do
  ! forall (i=1:tot_n) U(i)=0

  ! initial condition
  allocate(U(tot_n))
  U = 0
  BC = 0
  ! do i=1,2
  !   U(N_e_r*4*i+3:N_e_r*4*i+6)=1
  ! end do
  U(5:11) = 1

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
          call derivatives(vol_ngi(g,1),vol_ngi(g,2), ddxi_sh, ddeta_sh, e, N_e_r, dx, dy, co_ordinates, jac,&
                           det_jac,s_det_jac, tangent)
          M(inod,jnod) = M(inod,jnod) + vol_ngw(g)*sh_func(iloc)*sh_func(jloc)*det_jac

          ! ddx_sh_func = det_jac**(-1) *(jac(2,2)*ddxi_sh(iloc) - jac(1,2)*ddeta_sh(iloc))
          ! ddy_sh_func = det_jac**(-1) *(jac(1,1)*ddeta_sh(iloc) - jac(2,1)*ddxi_sh(iloc))

          ! K(inod,jnod) = K(inod,jnod) + vol_ngw(g)*dt*det_jac* (c(1)*sh_func(jloc)*ddx_sh_func + c(2)*sh_func(jloc)*ddy_sh_func)
        end do
      end do
    end do
  end do
  call FindInv(M, inv_M, tot_n, ErrorFlag)


open(unit=10, file='quad_points.txt')
 ! surface integration
  do e=1,totele
    call global_no(e, N_e_r, glob_no)
    do iface = 1,nface
      flux(iface) = 0
      K(iface) = 0
      do siloc=1,snloc   ! use all of the nodes not just the surface nodes.
        sinod = glob_no(siloc)
        do g=1,sngi
          call derivatives(s_ngi(g,1),s_ngi(g,2), ddxi_sh, ddeta_sh, e, N_e_r, dx, dy, co_ordinates, jac,&
                           det_jac, s_det_jac, tangent)
          call coordinates(e, N_e_r, dx, dy, co_ordinates, e_center)
          call shape_func(s_ngi(g,1),s_ngi(g,2), sh_func)

          ! normal to the iface
          snormal = cross(tangent(iface,:),domain_norm)

          ! vector from the centre of the element to a node on a boundary line
          r(1) = s_ngi(g,1) - e_center(1)
          r(2) = s_ngi(g,2) - e_center(2)
          r(3) = 0

          ! dot product of snormal and the vector joining e_center and a gaussian point
          s_dot = dot_product(snormal,r)
          if (s_dot <= 0 ) then
            snormal = snormal * (-1)
          end if

          ! unit normal to the face
          n_hat = sign(1.0,snormal)

          ! calculating flux at each quadrature point
          flux(iface) = flux(iface) + s_ngw(g) * s_det_jac(iface) * dt *dot_product(c,n_hat(1:2))&
                                                                   * sh_func(siloc)

          ddx_sh_func = det_jac**(-1) *(jac(2,2)*ddxi_sh(siloc) - jac(1,2)*ddeta_sh(siloc))
          ddy_sh_func = det_jac**(-1) *(jac(1,1)*ddeta_sh(siloc) - jac(2,1)*ddxi_sh(siloc))
          write(10,*) e, iface, siloc, g, s_ngw(g), s_det_jac(iface),dot_product(c,n_hat(1:2)),sh_func(siloc)

          ! calculating K for each corner node (siloc) using all quadrture points
          K(iface) = K(iface) + s_ngw(g)*dt*s_det_jac(iface)* (c(1)*ddx_sh_func + c(2)*ddy_sh_func)

        end do
      end do
    end do   ! do iface = 1, nface !  Between_Elements_And_Boundary
  end do   ! do ele = 1, totele ! Surface integral
! do i=1,tot_n
!   write(10,*) i, vec_K(i)
! end do
close(10)

  ! open(unit=10, file='quad_points.txt')
  ! ! print in matrix form
  ! ! do i=1,1
  ! !   do j=1,10
  ! !     write(10,'(f8.5)', advance='no') M(i,j)
  ! !   end do
  ! ! end do
  ! ! do g=1,4
  ! !   write(10,*) co_ordinates(g,1), co_ordinates(g,2)
  ! ! end do
  !
  ! ! write(10,*) ' x'
  ! ! do i=1,N_e_r*2
  ! ! write(10,*) x(i)
  ! ! end do
  !
  ! ! write(10,*) ' y'
  ! ! do i=1,N_e_c*2
  ! ! write(10,*) y(i)
  ! ! end do
  !
  ! write(10,*) 'element coordinates'
  ! ! do i=1,8
  ! !   write(10,*) s_ngi(i,1), s_ngi(i,2)
  ! ! end do
  ! write(10,*) 'U=', s_det_jac
  ! write(10,*)  'done'
  ! close(10)

  deallocate(U, x, y, x_coo, y_coo, x_dummy, y_dummy, M, BC)
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




subroutine derivatives(xi, eta, ddxi_sh, ddeta_sh, e, N_e_r, dx, dy, co_ordinates, jac, det_jac,s_det_jac, tangent)
  ! this module give derivatives of shape functions, x and y with respect to xi and eta
  implicit none
  real :: xi, eta, dx, dy, det_jac
  real :: s_det_jac(4), jac(2,2), ddxi_sh(4), ddeta_sh(4), e_center(3),tangent(4,3), co_ordinates(4,2)
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

  ! |J| of the element
  det_jac = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)
  ! |J| of the surface
  s_det_jac(1) = sqrt(jac(1,1)**2 + jac(1,2)**2)
  s_det_jac(2) = sqrt(jac(2,1)**2 + jac(2,2)**2)
  s_det_jac(3) = sqrt(jac(2,1)**2 + jac(2,2)**2)
  s_det_jac(4) = sqrt(jac(1,1)**2 + jac(1,2)**2)

  ! tangent to a face
  tangent(1,1) = jac(1,1)
  tangent(1,2) = jac(1,2)
  tangent(1,3) = 0
  tangent(2,1) = jac(2,1)
  tangent(2,2) = jac(2,2)
  tangent(2,3) = 0
  tangent(3,1) = jac(2,1)
  tangent(3,2) = jac(2,2)
  tangent(3,3) = 0
  tangent(4,1) = jac(1,1)
  tangent(4,2) = jac(1,2)
  tangent(4,3) = 0

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



subroutine FINDInv(matrix, inverse, n, errorflag)
  ! Returns the inverse of a matrix calculated by finding the LU
  ! decomposition.  Depends on LAPACK.
  !Subroutine to find the inverse of a square matrix
  !Author : Louisda16th a.k.a Ashwith J. Rego
  !Reference : Algorithm has been well explained in:
  !http://math.uww.edu/~mcfarlat/inverse.htm
  !http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
  !https://www.dreamincode.net/forums/topic/366231-FORTRAN-90%3A-Matrix-Inversion/

  IMPLICIT NONE
  !Declarations
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
  REAL, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
  REAL, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix

  LOGICAL :: FLAG = .TRUE.
  INTEGER :: i, j, k, l
  REAL :: m
  REAL, DIMENSION(n,2*n) :: augmatrix !augmented matrix

  !Augment input matrix with an identity matrix
  DO i = 1, n
      DO j = 1, 2*n
          IF (j <= n ) THEN
              augmatrix(i,j) = matrix(i,j)
          ELSE IF ((i+n) == j) THEN
              augmatrix(i,j) = 1
          Else
              augmatrix(i,j) = 0
          ENDIF
      END DO
  END DO

  !Reduce augmented matrix to upper traingular form
  DO k =1, n-1
      IF (augmatrix(k,k) == 0) THEN
          FLAG = .FALSE.
          DO i = k+1, n
              IF (augmatrix(i,k) /= 0) THEN
                  DO j = 1,2*n
                      augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                  END DO
                  FLAG = .TRUE.
                  EXIT
              ENDIF
              IF (FLAG .EQV. .FALSE.) THEN
                  PRINT*, "Matrix is non - invertible"
                  inverse = 0
                  errorflag = -1
                  return
              ENDIF
          END DO
      ENDIF
      DO j = k+1, n
          m = augmatrix(j,k)/augmatrix(k,k)
          DO i = k, 2*n
              augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
          END DO
      END DO
  END DO

  !Test for invertibility
  DO i = 1, n
      IF (augmatrix(i,i) == 0) THEN
          PRINT*, "Matrix is non - invertible"
          inverse = 0
          errorflag = -1
          return
      ENDIF
  END DO

  !Make diagonal elements as 1
  DO i = 1 , n
      m = augmatrix(i,i)
      DO j = i , (2 * n)
             augmatrix(i,j) = (augmatrix(i,j) / m)
      END DO
  END DO

  !Reduced right side half of augmented matrix to identity matrix
  DO k = n-1, 1, -1
      DO i =1, k
      m = augmatrix(i,k+1)
          DO j = k, (2*n)
              augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
          END DO
      END DO
  END DO

  !store answer
  DO i =1, n
      DO j = 1, n
          inverse(i,j) = augmatrix(i,j+n)
      END DO
  END DO
  errorflag = 0
end subroutine FINDinv
