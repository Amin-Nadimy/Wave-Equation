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
  real, dimension(2,4) :: loc_coordinate
  real :: c(2), tangent(4,3), snormal(3), e_center(3), r(3), s_dot, F, K
  real :: vol_ngi(9,2), vol_ngw(9), s_ngi(4,2), s_ngw(size(s_ngi)/2), s_sh_func(4,2)
  real,allocatable,dimension(:,:) :: M, inv_M, x_loc, x_all
  real,allocatable,dimension(:) :: U, Un, BC, x, y, x_coo, y_coo, x_dummy, y_dummy

  ! Costants
  ! volume quadrature points
  vol_ngi = transpose(reshape((/-sqrt(0.6), -sqrt(0.6), 0.0, -sqrt(0.6), sqrt(0.6), -sqrt(0.6), -sqrt(0.6),&
                               0.0, 0.0, 0.0, sqrt(0.6), 0.0, -sqrt(0.6),&
                               sqrt(0.6), 0.0, sqrt(0.6), sqrt(0.6), sqrt(0.6)/),(/2,9/)))

  ! weights of volume quadrature points
  vol_ngw = (/(real(25)/81), real(40)/81, real(25)/81, real(40)/81, real(64)/81, &
                          real(40)/81, real(25)/81, real(40)/81, real(25)/81/)

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
  ngi = size(vol_ngi)/2      ! total volume quadrature points
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

  call SHATRInew(L1, L2, L3, L4, weight, &
                 nloc,ngi,ndim,  n,nlx)

  do itime=1,ntime
    t_old =t_new

    do ele=1,totele
      ! volume integration
!AMIN is x_all 3D? ( , , )
      x_loc(1:ndim,:)=x_all(ele,1:ndim,:) ! x_all contains the coordinates of the corner nodes
      call det_nlx( x_loc, n, nlx, nx, detwei, weight, ndim,nloc,ngi )
      volume=sum(detwei)

      u_loc(:,:) = u_ele(:,:,ele) ! this is the
      t_loc(:) =  t_ele(:,ele) ! this is the FEM element - 2nd index is the local node number

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
        tgi(gi)=sum(t(gi,:)*t_loc(:))
      end do

      rhs_loc=0.0 ! the local to element rhs vector.
      mass_ele=0.0 ! initialize mass matric to be zero.
      do iloc=1,nloc
        inod = glob_no(iloc,ele)
        do jloc=1,nloc
          jnod = glob_no(ele,jloc)

          do gi=1,ngi
            !M(inod,jnod) = M(inod,jnod) + vol_ngw(g)*sh_func(iloc)*sh_func(jloc)*det_jac
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

        do iloc=1,nloc ! use all of the nodes not just the surface nodes.
          do idim=1,ndim
            xsgi(:,idim)  = xsgi(:,idim)  + sn(:,iloc)*x_loc(idim,iloc)
            usgi(:,idim)  = usgi(:,idim)  + sn(:,iloc)*u_loc(idim,iloc)
            usgi2(:,idim) = usgi2(:,idim) + sn2(:,iloc)*u_loc2(idim,iloc)

          end do
          tsgi(:)  = tsgi(:)  + sn(:,iloc)*t_loc(iloc)
          tsgi2(:) = tsgi2(:) + sn2(:,iloc)*t_loc2(iloc)
        end do
        ! start to do the integration...
        call det_snlx_all( nloc, sngi, ndim-1, ndim, x_loc, sn, snlx, sweigh, sdetwei, sarea, snorm, norm )

        do s_gi=1,sngi
          income(s_gi)=0.5 + 0.5*sign(1.0, sum(snorm(s_gi,:)*0.5*(usgi(s_gi,:)+usgi2(s_gi,:)))  )
        end do
        ! sloc_vec=0.0; sloc_vec2=0.0
        do iloc=1,nloc
          do idim=1,ndim
            s_cont(:) = snorm(:,idim)*sdetwei(:) &
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
              x_loc, x_all)
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
  real :: xi, eta, sh_func(4)

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
