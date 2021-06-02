subroutine coordinates(ele, no_ele_per_row, dx, dy, co_ordinates, ele_center,row, col)
  ! this subroutine gets no of elements in a row and ele .no
  ! and gives coordinates of 4 nodes of the ele
  implicit none
  integer, intent(in):: no_ele_per_row, ele
  integer :: col, row
  real, intent(in):: dx, dy
  real :: ele_center(3)
  real, dimension(4,2):: co_ordinates

  row = ceiling(real(ele)/no_ele_per_row)
  col = ele-(no_ele_per_row*((ceiling(real(ele)/no_ele_per_row))-1))
  co_ordinates(1,1) = (col-1)*dx
  co_ordinates(1,2) = dy*(row-1)
  co_ordinates(2,1) = dx*col
  co_ordinates(2,2) = dy*(row-1)
  co_ordinates(3,1) = dx*(col-1)
  co_ordinates(3,2) = dy*row
  co_ordinates(4,1) = dx*col
  co_ordinates(4,2) = dy*row

  ele_center(1) = 0.25* (co_ordinates(1,1) + co_ordinates(2,1) + co_ordinates(3,1) + co_ordinates(4,1))
  ele_center(2) = 0.25* (co_ordinates(1,2) + co_ordinates(2,2) + co_ordinates(3,2) + co_ordinates(4,2))
  ele_center(3) = 0

end subroutine coordinates

module solution_coordinates
  implicit none
  contains

  function U_position(ele, no_ele_per_row, no_ele_per_col) result(info)
    ! this function gives row and column number of U or any element
    implicit none
    integer, intent(in) :: ele, no_ele_per_row, no_ele_per_col
    integer :: col, row, info(2)
    row = int(ceiling(real(ele)/no_ele_per_row))
    col = ele-(no_ele_per_row*(row-1))
    info(1) = row
    info(2) = col
  end function U_position
end module solution_coordinates

subroutine global_no(ele, no_ele_per_row, glob_no)
  ! this subroutine gives global node numbers of an element
  implicit none
  integer, intent(in) :: ele, no_ele_per_row
  integer :: row, glob_no(4)

  row = ceiling(real(ele)/no_ele_per_row)
  glob_no(1) = (row-1)*2*no_ele_per_row + 2*(ele-1)+1
  glob_no(2) = (row-1)*2*no_ele_per_row + 2*(ele-1)+2
  glob_no(3) = (row-1)*2*no_ele_per_row+2*no_ele_per_row+2*(ele-1)+1
  glob_no(4) = (row-1)*2*no_ele_per_row+2*no_ele_per_row+2*(ele-1)+2

end subroutine global_no



subroutine sl_global_node(ele, s_node, no_ele_per_row)
  ! this subroutine give global node numbers of each edge of an element
  implicit none
  integer, intent(in) :: ele , no_ele_per_row
  integer :: s_node(4,2), row

  row = ceiling(real(ele)/no_ele_per_row)
  s_node(1,1) = (row-1)*2*no_ele_per_row + 2*(ele-1)+1
  s_node(1,2) = (row-1)*2*no_ele_per_row + 2*(ele-1)+2
  s_node(2,1) = s_node(1,2)
  s_node(2,2) = (row-1)*2*no_ele_per_row+2*no_ele_per_row+2*(ele-1)+2
  s_node(3,1) = s_node(1,1)
  s_node(3,2) = (row-1)*2*no_ele_per_row+2*no_ele_per_row+2*(ele-1)+1
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




subroutine derivatives(xi, eta, ddxi_sh, ddeta_sh, ele, no_ele_per_row, dx, dy, co_ordinates, jac, det_jac,surf_det_jac, tangent)
  ! this module give derivatives of shape functions, x and y with respect to xi and eta
  implicit none
  real, intent(in) :: xi, eta, dx, dy
  real :: det_jac
  real :: surf_det_jac(4), jac(2,2), ddxi_sh(4), ddeta_sh(4), ele_center(3),tangent(4,3), co_ordinates(4,2)
  integer, intent(in) :: ele, no_ele_per_row
  integer :: col, row

  ddxi_sh(1) = 0!-0.25*(1-eta)
  ddxi_sh(2) = 0!0.25*(1-eta)
  ddxi_sh(3) = 0!-0.25*(1+eta)
  ddxi_sh(4) = 0!0.25*(1+eta)

  ddeta_sh(1) = 0!-0.25*(1-xi)
  ddeta_sh(2) = 0!-0.25*(1+xi)
  ddeta_sh(3) = 0! 0.25*(1-xi)
  ddeta_sh(4) = 0! 0.25*(1+xi)

  call coordinates(ele, no_ele_per_row, dx, dy, co_ordinates, ele_center, row, col)
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
  surf_det_jac(1) = sqrt(jac(1,1)**2 + jac(1,2)**2)
  surf_det_jac(2) = sqrt(jac(2,1)**2 + jac(2,2)**2)
  surf_det_jac(3) = sqrt(jac(2,1)**2 + jac(2,2)**2)
  surf_det_jac(4) = sqrt(jac(1,1)**2 + jac(1,2)**2)

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

subroutine n_sign(surf_normal, unit_norm)
  implicit none
  integer :: i
  real, DIMENSION(3) :: unit_norm
  real, DIMENSION(3), INTENT(IN) :: surf_normal

  do i=1,3
    if (surf_normal(i).eq.0.0) then
      unit_norm(i) = 0.0
    else
      unit_norm(i) = abs(surf_normal(i)) / surf_normal(i)
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

  integer:: nface, sngi, no_timesteps, no_ele_per_row, no_ele_per_col, i, j, n, totele, ele, s_node(4,2)
  integer :: inod, jnod, no_sol_per_ele, no_surf_point, g, iloc, jloc, iface, siloc, sinod, ErrorFlag
  integer :: row, col, tot_unknowns
  integer, allocatable, dimension(:,:) :: U_pos
  integer, dimension(4):: glob_no
  real, dimension(3):: domain_norm, unit_norm
  real:: CFL, L, dx, dy, dt, xi, eta, det_jac, loc_flux(4), U_hat(4), tot_flux
  real :: sh_func(4),jac(2,2),surf_det_jac(4), ddxi_sh(4), ddeta_sh(4), ddx_sh_func, ddy_sh_func
  real :: s_sh_func(4,4)
  real, dimension(4,2) :: co_ordinates
  real :: velocity(2), tangent(4,3), surf_normal(3), ele_center(3), r(3), surf_dot
  real :: surf_gaus_point(4,2), surf_gaus_weight(size(surf_gaus_point)/2)
  real,allocatable,dimension(:,:) :: M, K, inv_M, U_plot
  real,allocatable,dimension(:) :: U, Un, x, y, BC

  ! Costants
  ! surface quadrature points
  surf_gaus_point = transpose(reshape((/0.0, -1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 1.0/),(/2,4/)))

  ! weights of surface quadrature points
  surf_gaus_weight = (/2,2,2,2/)
  no_sol_per_ele = 1  ! no of solution nodes in each element
  no_surf_point = 1 ! index of i in phi_i / no of solution in each surface
  nface = 4  ! no of faces of each elemenet
  sngi = size(surf_gaus_point)/2  ! no of surface quadrature points of the faces - this is set to the max no of all faces
  CFL = 0.05

  !velocity in x-dir(1) and y-dir(2)
  velocity(1) = 0.1
  velocity(2) = 0

  L = 0.5   ! length of the domain in each direction

  ! number of elements in each row (r) and column (c)
  no_ele_per_row = 15
  no_ele_per_col= 1
  no_timesteps = 10 ! number of timesteps

  ! normal to the domain
  domain_norm(1) = 0.0
  domain_norm(2) = 0.0
  domain_norm(3) = 1.0

  dx = L/(no_ele_per_row)
  dy = L/(no_ele_per_col)

  dt = CFL/((velocity(1)/dx)+(velocity(2)/dy))
  totele = no_ele_per_row * no_ele_per_col    ! total elemeno_timesteps
  tot_unknowns = totele

  allocate(x(totele), y(totele))
  allocate(U_pos(tot_unknowns,2))
  allocate(BC(2*(no_ele_per_row+no_ele_per_col)))
  allocate(U(tot_unknowns))
  allocate(U_plot(3,tot_unknowns))

  ! initial condition
  U = 0
  BC = 0
  ! do i=1,1
  !   U(no_ele_per_row/5+i*no_ele_per_row:no_ele_per_row/2+i*no_ele_per_row) = 1
  ! end do
  U(5:10)=1

  unit_norm=0 ! normal to an edge

  call sl_global_node(ele, s_node, no_ele_per_row)

  ! surface integration
  do n=1,no_timesteps
    Un = U
    do ele=1,totele
      call global_no(ele, no_ele_per_row, glob_no)
      do iface = 1,nface
        loc_flux(iface)=0
        do siloc=1,no_surf_point   ! use all of the nodes not just the surface nodes.
          ! sinod = e
          do g=1,sngi
            call derivatives(surf_gaus_point(g,1),surf_gaus_point(g,2), ddxi_sh, ddeta_sh, &
                             ele, no_ele_per_row, dx, dy, co_ordinates, jac,&
                             det_jac, surf_det_jac, tangent)
            call coordinates(ele, no_ele_per_row, dx, dy, co_ordinates, ele_center, row, col)
            call s_shape_func(iface, s_sh_func)
            ! normal to the iface
            surf_normal = cross(tangent(iface,:),domain_norm)
            ! vector from the centre of the element to a node on a boundary line
            r(1) = co_ordinates(iface,1) - ele_center(1)
            r(2) = co_ordinates(iface,2) - ele_center(2)
            r(3) = 0
            ! dot product of surf_normal and the vector joining e_center and a gaussian point
            surf_dot = dot_product(surf_normal,r)
            if (surf_dot <= 0 ) then
              surf_normal = surf_normal * (-1)
            end if
            ! unit normal to the face
            call n_sign(surf_normal, unit_norm)
            ! calculating flux at each quadrature point
            loc_flux(iface) = loc_flux(iface) + surf_gaus_weight(g) * surf_det_jac(iface) &
                                                * dt *dot_product(velocity,unit_norm(1:2))&
                                                                         * s_sh_func(iface,g)
          end do
        end do
      end do   ! face loop  Between_Elements_And_Boundary
      ! Upwind values for each surface =========================================
      U_pos(ele,:) = U_position(ele, no_ele_per_row, no_ele_per_col)

      if (U_pos(ele,1).eq.1) then
        U_hat(1) = BC(ele)
      else
        U_hat(1) = Un(ele-no_ele_per_row)
      end if

      U_hat(2) = Un(ele)

      if (U_pos(ele,2).eq.1) then
        U_hat(3) = BC(2*no_ele_per_row+U_pos(ele,1))
      else
        U_hat(3) = Un(ele-1)
      end if

      U_hat(4) = Un(ele)
      ! ========================================================================
      tot_flux = loc_flux(1)*U_hat(1) + loc_flux(2)*U_hat(2) + loc_flux(3)*U_hat(3) + loc_flux(4)*U_hat(4)   ! calculating total flux of element e
      U(ele) = Un(ele) - 1/(dx*dy) * tot_flux
    end do   ! element loop
    if (n.eq.1) then
      U_plot(1,:) = U
    elseif (n.eq.no_timesteps/2) then
      U_plot(2,:) = U
    elseif (n.eq.no_timesteps) then
      U_plot(3,:) = U
    end if
  end do ! time loop
!!!!!!!!!!!!!!!!!!!!!!!! axis info !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! in this block DG x and y coordinates are calculated
  do ele=1,totele
  call coordinates(ele, no_ele_per_row, dx, dy, co_ordinates, ele_center, row, col)
  x(ele) = ele_center(1)
  y(ele) = ele_center(2)
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
