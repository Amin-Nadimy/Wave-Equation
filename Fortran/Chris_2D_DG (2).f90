! 2D DG-FEM wave equation
program wave_equation
  use cross_product
  implicit none

  integer:: nface, sngi, nt, N_e_r, N_e_c, i, j, tot_n, totele, ngi, e, s_node(4,2)
  integer :: inod, jnod, nloc, snloc, g, iloc, jloc, iface, siloc, sinod
  integer :: n, row, col, U_hat(4), ndim, nonods
  integer, dimension(4):: sdot, glob_no
  integer, allocatable, dimension(:,:) :: U_pos
  real, dimension(3):: domain_norm, n_hat
  real:: CFL, L, dx, dy, dt, xi, eta, det_jac, flux(4), K_loc(4), M_sum(4)
  real :: sh_func(4),jac(2,2),s_det_jac(4), ddxi_sh(4), ddeta_sh(4), ddx_sh_func, ddy_sh_func
  real, dimension(:,:), allocatable :: loc_coordinate
  real :: c(2), tangent(4,3), snormal(3), e_center(3), r(3), s_dot, F, K
  real :: ngi_3p(9,2), ngi_3p_wei(9), ngi_2p_wei(4), s_ngi(4,2), ngi_2p(4,2)
  real :: s_ngw(size(s_ngi)/2), s_sh_func(4,2)
  real,allocatable,dimension(:,:) :: M, inv_M, x_loc, x_all
  real,allocatable,dimension(:) :: U, Un, BC, x, y, x_coo, y_coo, x_dummy, y_dummy

  ! Costants
  ! volume quadrature points 3 points
  ngi_3p = transpose(reshape((/-sqrt(0.6), -sqrt(0.6), 0.0, -sqrt(0.6), sqrt(0.6), -sqrt(0.6), -sqrt(0.6),&
                               0.0, 0.0, 0.0, sqrt(0.6), 0.0, -sqrt(0.6),&
                               sqrt(0.6), 0.0, sqrt(0.6), sqrt(0.6), sqrt(0.6)/),(/2,9/)))

  ngi_2p = transpose(reshape((/-sqrt(3)**-1, -sqrt(3)**-1, sqrt(3)**-1, -sqrt(3)**-1, &
                               -sqrt(3)**-1, sqrt(3)**-1, sqrt(3)**-1, sqrt(3)**-1/),(/2,4/)))

  ! weights of volume quadrature points
  ngi_3p_wei = (/(real(25)/81), real(40)/81, real(25)/81, real(40)/81, real(64)/81, &
                          real(40)/81, real(25)/81, real(40)/81, real(25)/81/)

  ngi_2p_wei = (/1,1,1,1/)

  ! surface quadrature points
  s_ngi = transpose(reshape((/0.0, -1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 1.0/),(/2,4/)))

  ! weights of surface quadrature points
  s_ngw = (/2,2,2,2/)
  ndim = 2 ! no of spatial dimensions e.g. 2
  nloc = 4  ! no of nodes in each element
  snloc = 2 ! no of nodes on each boundary line
  nface = 4  ! no of faces of each elemenet
  sngi = 2 ! no of surface quadrature points of the faces - this is set to the max no of all faces
  CFL = 0.05

  !velocity in x-dir(1) and y-dir(2)
  c(1) = 0.1
  c(2) = 0

  L = 0.5   ! length of the domain in each direction

  ! number of elements in each row (r) and column (c)
  N_e_r = 20
  N_e_c= 1
  nt = 10  ! number of timesteps

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
  ngi = size(ngi_3p)/2      ! total volume quadrature points
  ! array to store dot product of normal and r (vector from center of an element &
  ! to a point on its edge
  sdot=0
  n_hat=0 ! normal to an edge

  allocate(M(tot_n, tot_n))
  allocate(inv_M(tot_n, tot_n))
  allocate(BC(4*(N_e_r+N_e_c)))

  nonods=nloc*totele ! no of dg nodes...

!AMIN it isn't consistant with x_loc(1:ndim,:) in line 152
!AMIN also not consistant with line 152 x_all(ele, ,1:ndim,:)
  allocate(x_all(nonods,ndim))
  allocate(loc_coordinate(ndim,nloc))
  allocate(u_new(nonods),u_old(nonods))

  ! forall (i=1:tot_n) U(i)=0

  ! initial condition
  allocate(U(tot_n))
! initialize t_new
  U = 0
  BC = 0
  ! do i=1,2
  !   U(N_e_r*4*i+3:N_e_r*4*i+6)=1
  ! end do
  U(N_e_r*2/5:N_e_r) = 1
  U(N_e_r*2+N_e_r*2/5:N_e_r*2+N_e_r) = 1

  call sl_global_node(e, s_node, N_e_r)

  ngi=4
  allocate(l1(ngi), l2(ngi), l3(ngi), l4(ngi) )

  allocate(weight_space(ngi_space), n_space(ngi_space,nloc_space), nlx_space(ngi_space,ndim_space,nloc_space))

  allocate(weight_t(ngi_t), n_t(ngi_t,nloc_t), nlx_t(ngi_t,ndim_t,nloc_t) )l2(ngi), l3(ngi), l4(ngi) )

  allocate(x_loc(ndim,nloc), x_all(1,ndim,nloc))
!AMIN what does it do?
!  call SHATRInew(L1, L2, L3, L4, weight, &
!                 nloc,ngi,ndim,  n,nlx)
!
   mloc=1
   allocate(rdum(10))
   mdum(1:10)=0.0
   LOWQUA=.false.
   allocate(WEIGHT(ngi))
   call RE2DN4(LOWQUA,NGI,NLOC,MLOC,   M,WEIGHT,N,NLX(:,1,:),NLX(:,2,:),  SNGI,SNLOC,SWEIGH,SN,SNLX)

  do itime=1,ntime
    t_old =t_new

    do ele=1,totele
      ! volume integration
!AMIN! subroutine coordinates(e, N_e_r, dx, dy, loc_coordinate, e_center, row, col, ndim, nloc)
      x_loc(1:ndim,:)=x_all(ele,1:ndim,:) ! x_all contains the coordinates of the corner nodes
      call det_nlx( x_loc, n, nlx, nx, detwei, weight, ndim,nloc,ngi )
      volume=sum(detwei)

!AMIN what is U_ele?
      u_loc(:,:) = u_ele(:,:,ele) ! this is the
!AMIN changed to t_old
      t_loc(:) =  t_old(:,ele) ! this is the FEM element - 2nd index is the local node number

      ! calculates vel at gi sum = phi_1(gi)*c1 +phi_2(gi)*c2 +phi_3(gi)*c3 +phi_4(gi)*c4
      do idim=1,ndim
        do gi=1,ngi
           !ugi_x(gi,idim)=sum(nx(gi,idim,:)*u_loc(idim,:))
           ugi(gi,idim)=sum(n(gi,:)*u_loc(idim,:))
        end do
      end do

      ! calculates t at gi sum = phi_1(gi)*t1 +phi_2(gi)*t2 +phi_3(gi)*t3 +phi_4(gi)*t4
      do gi=1,ngi
       !ugi_x(gi,idim)=sum(nx(gi,idim,:)*u_loc(idim,:))
!AMIN n(gi,:) not t(gi,:)
        tgi(gi)=sum(n(gi,:)*t_loc(:))
      end do

      rhs_loc=0.0 ! the local to element rhs vector.
      mass_ele=0.0 ! initialize mass matric to be zero.
      do iloc=1,nloc
        inod = glob_no(iloc,ele)
        do jloc=1,nloc
          jnod = glob_no(ele,jloc)

          do gi=1,ngi
            !M(inod,jnod) = M(inod,jnod) + ngi_3p_wei(g)*sh_func(iloc)*sh_func(jloc)*det_jac
            Mass_ele(iloc,jloc) = Mass_ele(iloc,jloc) + n(gi,iloc)*n(gi,jloc)*detwei(gi)
            do idim=1,ndim
               rhs_loc(iloc) = rhs_loc(iloc) + nx(gi,idim,iloc)*ugi(gi,idim)*tgi(gi)*detwei(gi)
            end do
          end do ! quadrature
        end do ! ijloc
      end do ! iloc

      ! Include the surface integral here:
      do iface = 1,nface
        ele2 = face_ele( iface, ele)
!AMIN dimension isn't consistant with line 235 t_new(ele,iloc)
        t_loc2(:)=t_new(:,ele2)
        u_loc2(:,:)=u_new(:,:,ele2)

        !Surface integral along an element
        !if(ele>ele2) then ! perform integration along both sides...

        s_list_no = face_list_no( iface, ele)
        sn = face_sn(:,:,iface); snlx(:,:,:) = face_snlx(:,:,:,iface)
        sn2 =face_sn2(:,:,s_list_no)

        usgi=0.0; usgi2=0.0; xsgi=0.0; tsgi=0.0; tsgi2=0.0 
        do iloc=1,nloc ! use all of the nodes not just the surface nodes.
          do idim=1,ndim
!AMIN
            ! xsgi(:,idim)  = xsgi(:,idim)  + sn(:,iloc)*x_loc(idim,iloc)
            usgi(:,idim)  = usgi(:,idim)  + sn(:,iloc)*u_loc(idim,iloc)
            usgi2(:,idim) = usgi2(:,idim) + sn2(:,iloc)*u_loc2(idim,iloc)
            xsgi(:,idim)  = xsgi(:,idim)  + sn(:,iloc)*x_loc(idim,iloc)

          end do
          tsgi(:)  = tsgi(:)  + sn(:,iloc)*t_loc(iloc)
          tsgi2(:) = tsgi2(:) + sn2(:,iloc)*t_loc2(iloc)
        end do
! this is the approximate normal direction...
        do idim=1,ndim
           norm(idim) = sum(xsgi(:,idim))/real(sngi) - sum(x_loc(:,idim)/real(nloc))
        end do
        ! start to do the integration
! AMIN I can usit only once an dsave 4 values and use in later on
        call det_snlx_all( nloc, sngi, ndim-1, ndim, x_loc, sn, snlx, sweigh, sdetwei, sarea, snorm, norm )

        do s_gi=1,sngi
          income(s_gi)=0.5 + 0.5*sign(1.0, sum(snorm(s_gi,:)*0.5*(usgi(s_gi,:)+usgi2(s_gi,:)))  )
        end do
        ! sloc_vec=0.0; sloc_vec2=0.0
        do iloc=1,nloc
          do idim=1,ndim
            s_cont(:) = snorm(:,idim)*sdetwei(:) &
!AMIN got rid of (:)
                      *( (1.-income(:))* usgi(:,idim)*tsgi(:) + income(:)*usgi2(:,idim)*tsgi2(:) )
!AMIN- rhs_loc(iloc) not (ic,iloc)
            rhs_loc(iloc)  = rhs_loc(ic,iloc)  + sum( sn(:,iloc)*s_cont(:) )
            !sloc_vec2(ic,iloc) = sloc_vec2(ic,iloc) + sum( sn2(:,iloc)*s_cont(:) )
          end do
        end do
      end do ! iface

      mass_ele_inv=invert(mass_ele) ! inverse of the mass matric (nloc,nloc)
      do iloc=1,nloc
         t_new(ele,iloc)=t_old(ele,iloc) + dt*sum( mass_ele_inv(iloc,:) * rhs_loc(:) )
      end do

    end do ! do ele=1,totele
  end do ! do itime=1,ntime

  deallocate(U, x, y, x_coo, y_coo, x_dummy, y_dummy, M, BC, inv_M, U_pos&
              x_loc, x_all, loc_coordinate)
end program wave_equation





subroutine coordinates(e, N_e_r, dx, dy, loc_coordinate, e_center, row, col, ndim, nloc)
  ! this subroutine gets no of elements in a row and ele .no
  ! and gives coordinates of 4 nodes of the ele
  implicit none
  integer:: col, row, N_e_r, e, ndim, nloc
  real:: dx, dy, e_center(3)
  real,allocatable, dimension(:,:) :: loc_coordinate

  allocate(loc_coordinate(ndim, nloc))

  row = ceiling(real(e)/N_e_r)
  col = e-(N_e_r*((ceiling(real(e)/N_e_r))-1))
  loc_coordinate(1,1) = dx*(col-1)
  loc_coordinate(2,1) = dy*(row-1)
  loc_coordinate(1,2) = dx*col
  loc_coordinate(2,2) = dy*(row-1)
  loc_coordinate(1,3) = dx*(col-1)
  loc_coordinate(2,3) = dy*row
  loc_coordinate(1,4) = dx*col
  loc_coordinate(2,4) = dy*row

  e_center(1) = 0.25* (loc_coordinate(1,1) + loc_coordinate(1,2) + loc_coordinate(1,3) + loc_coordinate(1,4))
  e_center(2) = 0.25* (loc_coordinate(2,1) + loc_coordinate(2,2) + loc_coordinate(2,3) + loc_coordinate(2,4))
  e_center(3) = 0

end subroutine coordinates


function U_position(e, N_e_r, N_e_c) result(info)
  ! this function gives row and column number of U or any element
  implicit none
  integer, intent(in) :: e, N_e_r, N_e_c
  integer :: col, row, info(2)
  row = int(ceiling(real(e)/N_e_r))
  col = e-(N_e_r*(row-1))
  info(1) = row
  info(2) = col
end function U_position


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
  integer :: nloc, gi
  real :: xi, eta, sh_func(4)
  real, dimension(:,:), allocatable :: sh_func

  allocate(sh_func(gi,nloc))

  sh_func(1) = 0.25*(1-xi)*(1-eta)
  sh_func(2) = 0.25*(1+xi)*(1-eta)
  sh_func(3) = 0.25*(1-xi)*(1+eta)
  sh_func(4) = 0.25*(1+xi)*(1+eta)

end subroutine shape_func


subroutine s_shape_func(iface, xi, eta, s_sh_func)
  ! this subroutine contains shape functions
  implicit none
  integer :: i,j
  integer, intent(in) :: iface
  real, intent(inout):: s_sh_func(4,2),xi, eta

  s_sh_func = 0
  s_sh_func(1,1) = 0.25*(1-xi)*(1-eta)
  s_sh_func(1,2) = 0.25*(1+xi)*(1-eta)
  s_sh_func(2,1) = 0.25*(1+xi)*(1-eta)
  s_sh_func(2,2) = 0.25*(1+xi)*(1+eta)
  s_sh_func(3,1) = 0.25*(1-xi)*(1-eta)
  s_sh_func(3,2) = 0.25*(1-xi)*(1+eta)
  s_sh_func(4,1) = 0.25*(1-xi)*(1+eta)
  s_sh_func(4,2) = 0.25*(1+xi)*(1+eta)

end subroutine s_shape_func


subroutine derivatives(xi, eta, ddxi_sh, ddeta_sh, e, N_e_r, dx, dy, loc_coordinate, jac, det_jac,s_det_jac, tangent)
  ! this module give derivatives of shape functions, x and y with respect to xi and eta
  implicit none
  real :: xi, eta, dx, dy, det_jac
  real :: s_det_jac(4), jac(2,2), ddxi_sh(4), ddeta_sh(4), e_center(3),tangent(4,3), loc_coordinate(2,4)
  integer :: e, N_e_r, col, row

  ddxi_sh(1) = -0.25*(1-eta)
  ddxi_sh(2) =  0.25*(1-eta)
  ddxi_sh(3) = -0.25*(1+eta)
  ddxi_sh(4) =  0.25*(1+eta)

  ddeta_sh(1) = -0.25*(1-xi)
  ddeta_sh(2) = -0.25*(1+xi)
  ddeta_sh(3) =  0.25*(1-xi)
  ddeta_sh(4) =  0.25*(1+xi)

  call coordinates(e, N_e_r, dx, dy, loc_coordinate, e_center, row, col, ndim, nloc)
  ! dx_dxi
  jac(1,1) = 0.25*((eta-1)*loc_coordinate(1,1) + (1-eta)*loc_coordinate(1,2) - (1+eta)*loc_coordinate(1,3) + (1+eta)*loc_coordinate(1,4))
  ! dy_dxi
  jac(1,2) = 0.25*((eta-1)*loc_coordinate(2,1) + (1-eta)*loc_coordinate(2,2) - (1+eta)*loc_coordinate(2,3) + (1+eta)*loc_coordinate(2,4))
  ! dx_deta
  jac(2,1) = 0.25*((xi-1)*loc_coordinate(1,1) - (1+xi)*loc_coordinate(1,2) + (1-xi)*loc_coordinate(1,3) + (1+xi)*loc_coordinate(1,4))
  ! dy_deta
  jac(2,2) = 0.25*((xi-1)*loc_coordinate(2,1) - (1+xi)*loc_coordinate(2,2) + (1-xi)*loc_coordinate(2,3) + (1+xi)*loc_coordinate(2,4))

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



subroutine det_nlx( x_all, n, nlx, nx, detwei, weight, ndim,nloc,ngi )
! ****************************************************
! This sub form the derivatives of the shape functions
! ****************************************************
! x_all: spatial nodes.
! n, nlx, nlx_lxx: shape function and local derivatives of the shape functions (nlx_lxx is local grad of the local laplacian)- defined by shape functional library.
! nx, nx_lxx: derivatives of the shape functions.
! detwei, inv_jac: determinant at the quadrature pots and inverse of Jacobian at quadrature pts.
! ndim,nloc,ngi: no of dimensions, no of local nodes within an element, no of quadrature points.
! nlx_nod, nx_nod: same as nlx and nx but formed at the nodes not quadrature points.
  implicit none
  integer, intent( in ) :: ndim,nloc,ngi

  REAL, DIMENSION( ndim,nloc ), intent( in ) :: x_all
  REAL, DIMENSION( ngi, nloc ), intent( in ) :: N
  REAL, DIMENSION( ngi, ndim, nloc ), intent( in ) :: nlx
  REAL, DIMENSION( ngi ), intent( in ) :: WEIGHT
  REAL, DIMENSION( ngi ), intent( inout ) :: DETWEI
  REAL, DIMENSION( ngi, ndim, nloc ), intent( inout ) :: nx
  ! Local variables
  REAL :: AGI, BGI, CGI, DGI, EGI, FGI, GGI, HGI, KGI, A11, A12, A13, A21, &
          A22, A23, A31, A32, A33, DETJ
  INTEGER :: GI, L, IGLX, ii

  ! conventional:
  do  GI=1,NGI! Was loop 331

    AGI=0.
    BGI=0.
    CGI=0.

    DGI=0.
    EGI=0.
    FGI=0.

    GGI=0.
    HGI=0.
    KGI=0.

    do  L=1,NLOC! Was loop 79
      IGLX=L
      !ewrite(3,*)'xndgln, x, nl:', &
      !     iglx, l, x(iglx), y(iglx), z(iglx), NLX(L,GI), NLY(L,GI), NLZ(L,GI)
      ! NB R0 does not appear here although the z-coord might be Z+R0.
      AGI=AGI+NLX(GI,1,L)*X_ALL(1,IGLX)
      BGI=BGI+NLX(GI,1,L)*X_ALL(2,IGLX)
      CGI=CGI+NLX(GI,1,L)*X_ALL(3,IGLX)

      DGI=DGI+NLX(GI,2,L)*X_ALL(1,IGLX)
      EGI=EGI+NLX(GI,2,L)*X_ALL(2,IGLX)
      FGI=FGI+NLX(GI,2,L)*X_ALL(3,IGLX)

      GGI=GGI+NLX(GI,3,L)*X_ALL(1,IGLX)
      HGI=HGI+NLX(GI,3,L)*X_ALL(2,IGLX)
      KGI=KGI+NLX(GI,3,L)*X_ALL(3,IGLX)
    end do ! Was loop 79

    DETJ=AGI*(EGI*KGI-FGI*HGI)&
        -BGI*(DGI*KGI-FGI*GGI)&
        +CGI*(DGI*HGI-EGI*GGI)
    DETWEI(GI)=ABS(DETJ)*WEIGHT(GI)
    !ewrite(3,*)'gi, detj, weight(gi)', gi, detj, weight(gi)
    !rsum = rsum + detj
    !rsumabs = rsumabs + abs( detj )
    ! For coefficient in the inverse mat of the jacobian.
    A11= (EGI*KGI-FGI*HGI) /DETJ
    A21=-(DGI*KGI-FGI*GGI) /DETJ
    A31= (DGI*HGI-EGI*GGI) /DETJ

    A12=-(BGI*KGI-CGI*HGI) /DETJ
    A22= (AGI*KGI-CGI*GGI) /DETJ
    A32=-(AGI*HGI-BGI*GGI) /DETJ

    A13= (BGI*FGI-CGI*EGI) /DETJ
    A23=-(AGI*FGI-CGI*DGI) /DETJ
    A33= (AGI*EGI-BGI*DGI) /DETJ

    do  L=1,NLOC! Was loop 373
      NX(GI,1,L)= A11*NLX(GI,1,L)+A12*NLX(2,L,GI)+A13*NLX(GI,3,L)
      NX(GI,2,L)= A21*NLX(GI,1,L)+A22*NLX(2,L,GI)+A23*NLX(GI,3,L)
      NX(GI,3,L)= A31*NLX(GI,1,L)+A32*NLX(2,L,GI)+A33*NLX(GI,3,L)
    end do ! Was loop 373

  end do ! GI Was loop 331
end subroutine det_nlx



!
       SUBROUTINE RE2DN4(LOWQUA,NGI,NLOC,MLOC,&
     &      M,WEIGHT,N,NLX,NLY, &
     &      SNGI,SNLOC,SWEIGH,SN,SNLX)
!       use FLDebug
       IMPLICIT NONE
! NB might have to define surface elements for p and (u,v,w) 
! in here as well. 
!      This subroutine defines the shape functions M and N and their
!      derivatives at the Gauss points
!       REAL M(1,NGI),WEIGHT(NGI),N(4,NGI),NLX(4,NGI),NLY(4,NGI)
       INTEGER NGI,NLOC,MLOC
       REAL M(MLOC,NGI),WEIGHT(NGI)
       REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI)
       REAL POSI,TLY
       REAL LX(16),LY(16),LXP(4),LYP(4)
       REAL WEIT(16)
       INTEGER SNGI,SNLOC
       REAL SWEIGH(SNGI)
       REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI)
       INTEGER P,Q,CORN,GPOI,ILOC,JLOC,NDGI
       LOGICAL LOWQUA,GETNDP
       INTEGER I
! NB LXP(I) AND LYP(I) ARE THE LOCAL X AND Y COORDS OF NODAL POINT I
       
!       ewrite(3,*)'inside re2dn4, nloc,mloc,ngi',&
!     &               nloc,mloc,ngi
!
       LXP(1)=-1
       LYP(1)=-1
!
       LXP(2)=1
       LYP(2)=-1
!
!       LXP(3)=1
!       LYP(3)=1
!
!       LXP(4)=-1
!       LYP(4)=1
!
       LXP(3)=-1
       LYP(3)=1
!
       LXP(4)=1
       LYP(4)=1
!
       IF(NGI.EQ.4) THEN
          POSI=1.0/SQRT(3.0)
          LX(1)=-POSI
          LY(1)=-POSI
          LX(2)= POSI
          LY(2)= POSI
!
      do  Q=1,2! Was loop 23
      do  P=1,2! Was loop 24
      do  CORN=1,4! Was loop 25
                   GPOI=(Q-1)*2 + P
!
                   IF(MLOC.EQ.1)  M(1,GPOI)=1.
                   WEIGHT(GPOI)=1.
!
                   N(CORN,GPOI)=0.25*(1.+LXP(CORN)*LX(P))&
     &                   *(1.+LYP(CORN)*LY(Q))
                   NLX(CORN,GPOI)=0.25*LXP(CORN)*(1.+LYP(CORN)*LY(Q))
                   NLY(CORN,GPOI)=0.25*LYP(CORN)*(1.+LXP(CORN)*LX(P))
      end do ! Was loop 25
!
      end do ! Was loop 24
      end do ! Was loop 23
!       ewrite(3,*) 'here 1'
!       ewrite(3,*) 'N:',N
!       ewrite(3,*) 'NLX:',NLX
!       ewrite(3,*) 'NLY:',NLY
! Surface shape functions
          IF((SNGI.GT.1).AND.(SNLOC.GT.1)) THEN
!       ewrite(3,*) '***************** SNGI=',SNGI
      do  P=1,2! Was loop 27
      do  CORN=1,2! Was loop 27
                   GPOI=P
                   SN(CORN,GPOI)=0.5*(1.+LXP(CORN)*LX(P))
                   SNLX(CORN,GPOI)=0.5*LXP(CORN)
                   SWEIGH(GPOI)=1.
      end do ! Was loop 27
      end do ! Was loop 27
          ENDIF
! IF(NGI.EQ.4) THEN ...
       ELSE
          NDGI =INT(SQRT(NGI+0.1) +0.1)
!          ewrite(3,*) 'ndgi,ngi,sngi:',ndgi,ngi,sngi
!
          GETNDP=.FALSE.
          CALL LAGROT(WEIT,LX,NDGI,GETNDP)
          LY(1:NDGI) = LX(1:NDGI)
!          ewrite(3,*) 'weit:',weit
!          ewrite(3,*) 'lx:',lx
!
      do  Q=1,NDGI! Was loop 323
      do  P=1,NDGI! Was loop 324
      do  CORN=1,4! Was loop 325
!            ewrite(3,*) 'q,p,corn:',q,p,corn
                   GPOI=(Q-1)*NDGI + P
!
                   IF(MLOC.EQ.1)  M(1,GPOI)=1.
                   WEIGHT(GPOI)=WEIT(P)*WEIT(Q)
!             ewrite(3,*) 'here1'
!
                   N(CORN,GPOI)=0.25*(1.+LXP(CORN)*LX(P))&
     &                   *(1.+LYP(CORN)*LY(Q))
!             ewrite(3,*) 'here2'
                   NLX(CORN,GPOI)=0.25*LXP(CORN)*(1.+LYP(CORN)*LY(Q))
                   NLY(CORN,GPOI)=0.25*LYP(CORN)*(1.+LXP(CORN)*LX(P))
!             ewrite(3,*) 'here3'
      end do ! Was loop 325
!
      end do ! Was loop 324
      end do ! Was loop 323
!       ewrite(3,*) 'here 1'
!       ewrite(3,*) 'N:',N
!       ewrite(3,*) 'NLX:',NLX
!       ewrite(3,*) 'NLY:',NLY
! Surface shape functions
!       ewrite(3,*) '***************** SNGI=',SNGI
          IF(SNGI.GT.0) THEN
             GETNDP=.FALSE.
             CALL LAGROT(WEIT,LX,SNGI,GETNDP)
      do  P=1,SNGI! Was loop 327
      do  CORN=1,2! Was loop 327
                   GPOI=P
                   SN(CORN,GPOI)=0.5*(1.+LXP(CORN)*LX(P))
                   SNLX(CORN,GPOI)=0.5*LXP(CORN)
                   SWEIGH(GPOI)=WEIT(P)
      end do ! Was loop 327
      end do ! Was loop 327
! ENDOF IF(SNGI.GT.0) THEN...
          ENDIF
!
! END OF IF(NGI.EQ.4) THEN ELSE ...
       ENDIF
!
       IF(MLOC.EQ.NLOC) THEN
      do  I=1,4! Was loop 2545
      do  CORN=1,4! Was loop 2545
                M(CORN,I)=N(CORN,I)
      end do ! Was loop 2545
      end do ! Was loop 2545
       ENDIF
!       ewrite(3,*) 'in re2dn4.f here 2 ngi,sngi',ngi,sngi
!       ewrite(3,*) 'N:',N
!       ewrite(3,*) 'NLX:',NLX
!       ewrite(3,*) 'NLY:',NLY
       END





      SUBROUTINE LAGROT(WEIT,QUAPOS,NDGI,GETNDP)
!        use RGPTWE_module
      IMPLICIT NONE
!     This computes the weight and points for standard Gaussian quadrature.
!     IF(GETNDP) then get the POSITION OF THE NODES 
!     AND DONT BOTHER WITH THE WEITS.
      INTEGER NDGI
      REAL WEIT(NDGI),QUAPOS(NDGI)
      LOGICAL GETNDP
      LOGICAL WEIGHT
      INTEGER IG
! real function...
      real RGPTWE
!     
      IF(.NOT.GETNDP) THEN
         WEIGHT=.TRUE.
         do IG=1,NDGI
            WEIT(IG)=RGPTWE(IG,NDGI,WEIGHT)
         END DO
!     
         WEIGHT=.FALSE.
         do IG=1,NDGI
            QUAPOS(IG)=RGPTWE(IG,NDGI,WEIGHT)
         END DO
      ELSE
         IF(NDGI.EQ.1) THEN
            QUAPOS(1)=0.
         ELSE
            do IG=1,NDGI
               QUAPOS(IG)= -1+2.*REAL(IG-1)/REAL(NDGI-1)
            END DO
         ENDIF
      ENDIF
      END SUBROUTINE LAGROT





  SUBROUTINE det_snlx_all( SNLOC, SNGI, SNDIM, ndim, XSL_ALL, SN, SNLX, SWEIGH, SDETWE, SAREA, NORMXN_ALL, NORMX_ALL )
!       inv_jac )
    IMPLICIT NONE

    INTEGER, intent( in ) :: SNLOC, SNGI, SNDIM, ndim
    REAL, DIMENSION( NDIM, SNLOC ), intent( in ) :: XSL_ALL
    REAL, DIMENSION( SNGI, SNLOC ), intent( in ) :: SN
    REAL, DIMENSION( SNGI, SNDIM, SNLOC ), intent( in ) :: SNLX
    REAL, DIMENSION( SNGI ), intent( in ) :: SWEIGH
    REAL, DIMENSION( SNGI ), intent( inout ) :: SDETWE
    REAL, intent( inout ) ::  SAREA
    REAL, DIMENSION( sngi, NDIM ), intent( inout ) :: NORMXN_ALL
    REAL, DIMENSION( NDIM ), intent( in ) :: NORMX_ALL
!    REAL, DIMENSION( NDIM,ndim ), intent( in ) :: inv_jac
    ! Local variables
    INTEGER :: GI, SL, IGLX
    REAL :: DXDLX, DXDLY, DYDLX, DYDLY, DZDLX, DZDLY
    REAL :: A, B, C, DETJ, RUB3, RUB4

    SAREA=0.

       DO GI=1,SNGI

          DXDLX=0.
          DXDLY=0.
          DYDLX=0.
          DYDLY=0.
          DZDLX=0.
          DZDLY=0.

          DO SL=1,SNLOC
             DXDLX=DXDLX + SNLX(GI,1,SL)*XSL_ALL(1,SL)
             DXDLY=DXDLY + SNLX(GI,2,SL)*XSL_ALL(1,SL)
             DYDLX=DYDLX + SNLX(GI,1,SL)*XSL_ALL(2,SL)
             DYDLY=DYDLY + SNLX(GI,2,SL)*XSL_ALL(2,SL)
             DZDLX=DZDLX + SNLX(GI,1,SL)*XSL_ALL(3,SL)
             DZDLY=DZDLY + SNLX(GI,2,SL)*XSL_ALL(3,SL)
          END DO
          A = DYDLX*DZDLY - DYDLY*DZDLX
          B = DXDLX*DZDLY - DXDLY*DZDLX
          C = DXDLX*DYDLY - DXDLY*DYDLX

          DETJ=SQRT( A**2 + B**2 + C**2)
!          inv_jac(1,1)=DXDLX; inv_jac(1,2)=DXDLY; inv_jac(1,3)=DXDLZ
!          inv_jac(2,1)=DyDLX; inv_jac(2,2)=DyDLY; inv_jac(2,3)=DyDLZ
!          inv_jac(3,1)=DzDLX; inv_jac(3,2)=DzDLY; inv_jac(3,3)=DzDLZ
!          inv_jac=inv_jac/detj 
          SDETWE(GI)=DETJ*SWEIGH(GI)
          SAREA=SAREA+SDETWE(GI)

          ! Calculate the normal at the Gauss pts...
          ! Perform x-product. N=T1 x T2
          CALL NORMGI(NORMXN_ALL(GI,1),NORMXN_ALL(GI,2),NORMXN_ALL(GI,3), &
               DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,DZDLY, &
               NORMX_ALL(1),NORMX_ALL(2),NORMX_ALL(3))
       END DO

    RETURN

  END SUBROUTINE det_snlx_all



  SUBROUTINE NORMGI( NORMXN, NORMYN, NORMZN, &
       DXDLX, DYDLX, DZDLX, DXDLY, DYDLY, DZDLY, &
       NORMX, NORMY, NORMZ)
    ! Calculate the normal at the Gauss pts
    ! Perform x-product. N=T1 x T2
    implicit none
    REAL, intent( inout ) :: NORMXN, NORMYN, NORMZN
    REAL, intent( in )    :: DXDLX, DYDLX, DZDLX, DXDLY, DYDLY, DZDLY
    REAL, intent( in )    :: NORMX, NORMY, NORMZ
    ! Local variables
    REAL :: RN, SIRN

    CALL XPROD1( NORMXN, NORMYN, NORMZN, &
         DXDLX, DYDLX, DZDLX, &
         DXDLY, DYDLY, DZDLY )

    RN = SQRT( NORMXN**2 + NORMYN**2 + NORMZN**2 )

    SIRN = SIGN( 1.0 / RN, NORMXN * NORMX + NORMYN * NORMY + NORMZN * NORMZ )

    NORMXN = SIRN * NORMXN
    NORMYN = SIRN * NORMYN
    NORMZN = SIRN * NORMZN

    RETURN

  END SUBROUTINE NORMGI



  SUBROUTINE XPROD1( AX, AY, AZ, &
       BX, BY, BZ, &
       CX, CY, CZ )
    implicit none
    REAL, intent( inout ) :: AX, AY, AZ
    REAL, intent( in )    :: BX, BY, BZ, CX, CY, CZ

    ! Perform x-product. a=b x c
    AX =    BY * CZ - BZ * CY
    AY = -( BX * CZ - BZ * CX )
    AZ =    BX * CY - BY * CX

    RETURN
  END subroutine XPROD1


