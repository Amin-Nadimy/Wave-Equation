subroutine coordinates(e, N_e_r, dx, dy, co_ordinates, e_center,row, col)
  ! this subroutine gets no of elements in a row and ele .no
  ! and gives coordinates of 4 nodes of the ele
  implicit none
  integer, intent(in):: N_e_r, e
  integer :: col, row
  real, intent(in):: dx, dy
  real :: e_center(3)
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

module solution_coordinates
  implicit none
  contains

  function U_coordinates(e, N_e_r, N_e_c) result(info)
    ! this function gives row and column number of U or any element
    implicit none
    integer, intent(in) :: e, N_e_r, N_e_c
    integer :: col, row, info(2)
    row = int(ceiling(real(e)/N_e_r))
    col = e-(N_e_r*(row-1))
    info(1) = row
    info(2) = col
  end function U_coordinates
end module solution_coordinates

subroutine global_no(e, N_e_r, glob_no)
  ! this subroutine gives global node numbers of an element
  implicit none
  integer, intent(in) :: e, N_e_r
  integer :: row, glob_no(4)

  row = ceiling(real(e)/N_e_r)
  glob_no(1) = (row-1)*2*N_e_r + 2*(e-1)+1
  glob_no(2) = (row-1)*2*N_e_r + 2*(e-1)+2
  glob_no(3) = (row-1)*2*N_e_r+2*N_e_r+2*(e-1)+1
  glob_no(4) = (row-1)*2*N_e_r+2*N_e_r+2*(e-1)+2

end subroutine global_no



subroutine sl_global_node(e, s_node, N_e_r)
  ! this subroutine give global node numbers of each edge of an element
  implicit none
  integer, intent(in) :: e , N_e_r
  integer :: s_node(4,2), row

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
  real, intent(in) :: xi, eta
  real :: sh_func(4)

  sh_func(1) = 1!1/4*(1-xi)*(1-eta)
  sh_func(2) = 1!1/4*(1+xi)*(1-eta)
  sh_func(3) = 1!1/4*(1-xi)*(1+eta)
  sh_func(4) = 1!1/4*(1+xi)*(1+eta)

end subroutine shape_func

subroutine s_shape_func(iface, s_sh_func)
  ! this subroutine contains shape functions
  implicit none
  integer :: i,j
  integer, intent(in) :: iface
  real, intent(inout):: s_sh_func(4,4)

  s_sh_func=0
  s_sh_func(iface, iface) = 1

end subroutine s_shape_func




subroutine derivatives(xi, eta, ddxi_sh, ddeta_sh, e, N_e_r, dx, dy, co_ordinates, jac, det_jac,s_det_jac, tangent)
  ! this module give derivatives of shape functions, x and y with respect to xi and eta
  implicit none
  real, intent(in) :: xi, eta, dx, dy
  real :: det_jac
  real :: s_det_jac(4), jac(2,2), ddxi_sh(4), ddeta_sh(4), e_center(3),tangent(4,3), co_ordinates(4,2)
  integer, intent(in) :: e, N_e_r
  integer :: col, row

  ddxi_sh(1) = 0!-0.25*(1-eta)
  ddxi_sh(2) = 0!0.25*(1-eta)
  ddxi_sh(3) = 0!-0.25*(1+eta)
  ddxi_sh(4) = 0!0.25*(1+eta)

  ddeta_sh(1) = 0!-0.25*(1-xi)
  ddeta_sh(2) = 0!-0.25*(1+xi)
  ddeta_sh(3) = 0! 0.25*(1-xi)
  ddeta_sh(4) = 0! 0.25*(1+xi)

  call coordinates(e, N_e_r, dx, dy, co_ordinates, e_center, row, col)
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

subroutine n_sign(snormal, n_hat)
  implicit none
  integer :: i
  real, DIMENSION(3) :: n_hat
  real, DIMENSION(3), INTENT(IN) :: snormal

  do i=1,3
    if (snormal(i).eq.0.0) then
      n_hat(i) = 0.0
    else
      n_hat(i) = abs(snormal(i)) / snormal(i)
    end if
  end do
END subroutine n_sign


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
  INTEGER, INTENT(IN) :: n !size of the squre ,matrix
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

!===============================================================================
!===============================================================================
! 2D DG-FEM wave equation
program wave_equation
  use cross_product
  use solution_coordinates
  implicit none

  integer:: nface, sngi, nt, N_e_r, N_e_c, i, j, n, totele, e, s_node(4,2)
  integer :: inod, jnod, nloc, snloc, g, iloc, jloc, iface, siloc, sinod, ErrorFlag
  integer :: row, col, tot_unknowns
  integer, allocatable, dimension(:,:) :: U_coo
  integer, dimension(4):: sdot, glob_no
  real, dimension(3):: domain_norm, n_hat
  real:: CFL, L, dx, dy, dt, xi, eta, det_jac, flux(4), U_hat(4), F
  real :: sh_func(4),jac(2,2),s_det_jac(4), ddxi_sh(4), ddeta_sh(4), ddx_sh_func, ddy_sh_func
  real :: s_sh_func(4,4)
  real, dimension(4,2) :: co_ordinates
  real :: c(2), tangent(4,3), snormal(3), e_center(3), r(3), s_dot
  real :: s_ngi(4,2), s_ngw(size(s_ngi)/2)
  real,allocatable,dimension(:,:) :: M, K, inv_M, U_plot
  real,allocatable,dimension(:) :: U, Un, x, y, BC

  ! Costants
  ! surface quadrature points
  s_ngi = transpose(reshape((/0.0, -1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 1.0/),(/2,4/)))

  ! weights of surface quadrature points
  s_ngw = (/2,2,2,2/)
  nloc = 1  ! no of solution nodes in each element
  snloc = 1 ! index of i in phi_i / no of solution in each surface
  nface = 4  ! no of faces of each elemenet
  sngi = size(s_ngi)/2  ! no of surface quadrature points of the faces - this is set to the max no of all faces
  CFL = 0.05

  !velocity in x-dir(1) and y-dir(2)
  c(1) = 0.1
  c(2) = 0.1

  L = 0.5   ! length of the domain in each direction

  ! number of elements in each row (r) and column (c)
  N_e_r = 20
  N_e_c= 21
  nt = 500 ! number of timesteps

  ! normal to the domain
  domain_norm(1) = 0.0
  domain_norm(2) = 0.0
  domain_norm(3) = 1.0

  dx = L/(N_e_r)
  dy = L/(N_e_c)

  dt = CFL/((c(1)/dx)+(c(2)/dy))
  totele = N_e_r * N_e_c    ! total element
  tot_unknowns = totele

  allocate(x(totele), y(totele))
  allocate(U_coo(tot_unknowns,2))
  allocate(BC(2*(N_e_r+N_e_c)))
  allocate(U(tot_unknowns))
  allocate(U_plot(3,tot_unknowns))

  ! initial condition
  U = 0
  BC = 0
  do i=1,20
    U(N_e_r/5+i*N_e_r:N_e_r/2+i*N_e_r) = 1
  end do

  ! array to store dot product of normal and r (vector from center of an element &
  ! to a point on its edge
  sdot=0
  n_hat=0 ! normal to an edge

  call sl_global_node(e, s_node, N_e_r)
BC
  ! surface integration
  do n=1,nt
    Un = U
    do e=1,totele
      call global_no(e, N_e_r, glob_no)
      do iface = 1,nface
        flux(iface)=0
        do siloc=1,snloc   ! use all of the nodes not just the surface nodes.
          ! sinod = e
          do g=1,sngi
            call derivatives(s_ngi(g,1),s_ngi(g,2), ddxi_sh, ddeta_sh, e, N_e_r, dx, dy, co_ordinates, jac,&
                             det_jac, s_det_jac, tangent)
            call coordinates(e, N_e_r, dx, dy, co_ordinates, e_center, row, col)
            call s_shape_func(iface, s_sh_func)
            ! normal to the iface
            snormal = cross(tangent(iface,:),domain_norm)
            ! vector from the centre of the element to a node on a boundary line
            r(1) = co_ordinates(iface,1) - e_center(1)
            r(2) = co_ordinates(iface,2) - e_center(2)
            r(3) = 0
            ! dot product of snormal and the vector joining e_center and a gaussian point
            s_dot = dot_product(snormal,r)
            if (s_dot <= 0 ) then
              snormal = snormal * (-1)
            end if
            ! unit normal to the face
            call n_sign(snormal, n_hat)
            ! calculating flux at each quadrature point
            flux(iface) = flux(iface) + s_ngw(g) * s_det_jac(iface) * dt *dot_product(c,n_hat(1:2))&
                                                                                  * s_sh_func(iface,g)

          end do
        end do
      end do   ! face loop  Between_Elements_And_Boundary
      ! Upwind values for each surface =========================================
      U_coo(e,:) = U_coordinates(e, N_e_r, N_e_c)

      if (U_coo(e,1).eq.1) then
        U_hat(1) = BC(e)
      else
        U_hat(1) = Un(e-N_e_r)
      end if

      U_hat(2) = Un(e)

      if (U_coo(e,2).eq.1) then
        U_hat(3) = BC(2*N_e_r+U_coo(e,1))
      else
        U_hat(3) = Un(e-1)
      end if

      U_hat(4) = Un(e)
      ! ========================================================================
      F = flux(1)*U_hat(1) + flux(2)*U_hat(2) + flux(3)*U_hat(3) + flux(4)*U_hat(4)   ! calculating total flux of element e
      U(e) = Un(e) - 1/(dx*dy) * F
    end do   ! element loop
    if (n.eq.1) then
      U_plot(1,:) = U
    elseif (n.eq.nt/2) then
      U_plot(2,:) = U
    elseif (n.eq.nt) then
      U_plot(3,:) = U
    end if
  end do ! time loop
!!!!!!!!!!!!!!!!!!!!!!!! axis info !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! in this block DG x and y coordinates are calculated
  do e=1,totele
  call coordinates(e, N_e_r, dx, dy, co_ordinates, e_center, row, col)
  x(e) = e_center(1)
  y(e) = e_center(2)
  end do
  !!!!!!!!!!!!!!!!!!!!!!!!save plot info !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! save the 1st timestep info
  open (unit=40, file='n1.dat')
  do i=1,tot_unknowns
    write(40,*) x(i), y(i), U_plot(1,i)
  end do
  close(40)

  ! save the mid-timestep info
  open (unit=41, file='n2.dat')
  do i=1,tot_unknowns
    write(41,*) x(i), y(i), U_plot(2,i)
  end do
  close(41)

  ! save the last timestep info
  open (unit=42, file='n3.dat')
  do i=1,tot_unknowns
    write(42,*) x(i), y(i), U_plot(3,i)
  end do
  close(42)

!!!!!!!!!!!!!!!!!!!!!!!!! plotting commands !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! gnuplot
  ! splot "n1.dat" title "first timestep" w l, "n2.dat" title "mid-timestep" w l, "n3.dat" title "last timestep" w l
  ! set xlabel "x"
  ! set ylabel "y"
  ! set zlabel "U"
  ! set zrange[0:1]

  deallocate(U, U_plot, BC, x, y)
end program wave_equation
