! 2D DG-FEM wave equation
! face_ele(iface, ele) = given the face no iface and element no return the element next to
! the surface or if negative return the negative of the surface element number between element ele and face iface.

program wave_equation
  implicit none

  logical :: LOWQUA
  integer :: nloc, gi, ngi, sngi, iface, nface, totele, ele, ele2, s_list_no, s_gi, iloc, jloc
  integer:: itime, ntime, idim, ndim, nonodes, snloc, mloc
  integer :: no_ele_col, no_ele_row, errorflag, max_face_list_no
  integer :: face_ele!(1,1)
  integer :: face_list_no!(:,:)

  real :: sarea, volume, dt, CFL, L, dx, dy

  real :: n(:,:), nlx(:,:,:), nx(:,:,:)
  real :: weight(:), detwei(:), sdetwei(:)
  real :: sn(:,:),sn2(:,:),snlx(:,:,:),sweigh(:), s_cont(:)
  real :: t_loc(:), t_loc2(:), t_old(:,:), t_new(:,:), tgi(:), tsgi(:), tsgi2(:)
  real :: face_sn(:,:,:), face_sn2(:,:,:), face_snlx(:,:,:,:), face_sweigh(:,:)
  real :: u_loc(:,:), u_ele(:,:,:), u_loc2(:,:), x_loc(:,:), ugi(:,:), u_new(:),u_old(:)
  real :: x_all(:,:,:)
  real :: xsgi(:,:), usgi(:,:), usgi2(:,:), income(:), snorm(:,:), norm(:)
  real :: mass_ele(:,:), mass_ele_inv(:,:), rhs_loc(:)
  ! real :: rdum(:)

  CFL = 0.05
  L = 0.5   ! length of the domain in each direction
  no_ele_row = 10
  no_ele_col = 1
  totele = no_ele_row * no_ele_col
  nface = totele * 4
  dx = L/(no_ele_row)
  dy = L/(no_ele_col)
  ndim = 2
  snloc = 2 ! no of nodes on each boundary line
  ngi = 4
  sngi = 2 ! no of surface quadrature points of the faces - this is set to the max no of all faces
  nloc = 4
  ntime = 10

  allocate(n( ngi, nloc ), nx( ngi, ndim, nloc ), nlx( ngi, ndim, nloc ))
  allocate(weight(ngi), detwei(ngi), sdetwei(sngi))
  !allocate(face_list_no( nface, totele))
!AMIN l1, l2, l3, l4 aren't used anywhere
  ! allocate(l1(ngi), l2(ngi), l3(ngi), l4(ngi))
  allocate(sn(sngi,nloc),sn2(sngi,nloc),snlx(sngi,ndim,nloc),sweigh(sngi))
  allocate(u_loc(ndim,nloc), u_loc2(ndim,nloc), ugi(ngi,ndim), u_ele(ndim,nloc,totele))
  allocate(xsgi(sngi,ndim), usgi(sngi,ndim), usgi2(sngi,ndim), income(sngi), snorm(sngi,ndim), norm(ndim))
  allocate(t_loc(nloc), t_loc2(nloc), t_old(nloc,totele), t_new(nloc,totele), tgi(ngi), tsgi(sngi), tsgi2(sngi))
  allocate(face_sn(sngi,nloc,nface), face_sn2(sngi,nloc,max_face_list_no), face_snlx(sngi,ndim,nloc,nface), face_sweigh(sngi,nface))
  allocate(x_all(ndim,nloc,totele), x_loc(ndim,nloc))
  allocate(mass_ele(nloc,nloc), mass_ele_inv(nloc,nloc), rhs_loc(nloc))
!AMIN, rdum isn't used anywhere
  ! allocate(rdum(10))


  mloc=1
!AMIN mdum isn't used anywhere
  ! mdum(1:10)=0.0
  LOWQUA=.false.
  ! initial conditions
  t_new(:,2:5) = 1
  U_ele(:,:,:) = 0.1
  dt = CFL/((u_ele(1,1,1)/dx)+(u_ele(2,1,1)/dy))

  call RE2DN4(LOWQUA,NGI,NLOC,MLOC,M,WEIGHT,N,NLX(:,1,:),NLX(:,2,:),  SNGI,SNLOC,SWEIGH,SN,SNLX)

  do itime=1,ntime
    t_old =t_new

    do ele=1,totele
      ! volume integration
      !x_loc(1:ndim,:)=x_all(ele,1:ndim,:) ! x_all contains the coordinates of the corner nodes
      x_loc(1:ndim,:)=x_all(1:ndim,:,ele) ! x_all contains the coordinates of the corner nodes
      ! call det_nlx( x_loc, n, nlx, nx, detwei, weight, ndim,nloc,ngi )
      !volume=sum(detwei)

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
        do idim=1,ndim
           ugi(gi,idim)=sum(n(gi,:)*u_loc(idim,:))
        end do
        tgi(gi)=sum(n(gi,:)*t_loc(:))
      end do

      rhs_loc=0.0 ! the local to element rhs vector.
      mass_ele=0.0 ! initialize mass matric to be zero.
      do iloc=1,nloc
        !inod = glob_no(iloc,ele)
        do jloc=1,nloc
          !jnod = glob_no(ele,jloc)

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
        call ele_neighbour(iface, totele, nface, ele, face_ele, no_ele_row, row, row2, s_list_no)
        ele2 = face_ele!( iface, ele) ! calculate this somehwere!
        t_loc2(:)=t_new(:,ele2)
        u_loc2(:,:)=u_new(:,:,ele2)

        !Surface integral along an element
        !if(ele>ele2) then ! perform integration along both sides...

!AMIN s_list_no is local face number - need to work on this also
        !s_list_no = face_list_no( iface, ele)
        sn = face_sn(:,:,iface); snlx(:,:,:) = face_snlx(:,:,:,iface)
        sn2 =face_sn2(:,:,s_list_no)

        usgi=0.0; usgi2=0.0; xsgi=0.0; tsgi=0.0; tsgi2=0.0
        do iloc=1,nloc ! use all of the nodes not just the surface nodes.
          do idim=1,ndim

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
        ! call det_snlx_all( nloc, sngi, ndim-1, ndim, x_loc, sn, snlx, sweigh, sdetwei, sarea, snorm, norm )

        do s_gi=1,sngi
          income(s_gi)=0.5 + 0.5*sign(1.0, sum(snorm(s_gi,:)*0.5*(usgi(s_gi,:)+usgi2(s_gi,:)))  )
        end do
        ! sloc_vec=0.0; sloc_vec2=0.0
        do iloc=1,nloc
          do idim=1,ndim
            s_cont(:) = snorm(:,idim)*sdetwei(:) &
!AMIN got rid of (:)
                      *( (1.-income(:))* usgi(:,idim)*tsgi(:) + income(:)*usgi2(:,idim)*tsgi2(:) )
            rhs_loc(iloc)  = rhs_loc(iloc)  + sum( sn(:,iloc)*s_cont(:) )
            !sloc_vec2(ic,iloc) = sloc_vec2(ic,iloc) + sum( sn2(:,iloc)*s_cont(:) )
          end do
        end do
      end do ! iface

      ! call FINDInv(mass_ele, mass_ele_inv, size(mass_ele(1,:)), errorflag)
      !mass_ele_inv=invert(mass_ele) ! inverse of the mass matric (nloc,nloc)
      do iloc=1,nloc
         t_new(iloc,ele)=t_old(iloc,ele) + dt*sum( mass_ele_inv(iloc,:) * rhs_loc(:) )
      end do
    end do ! do ele=1,totele
  end do ! do itime=1,ntime

  deallocate(x_loc, x_all, n, nx, nlx, weight, detwei, sdetwei,&
              sn, sn2, snlx, sweigh, u_loc, u_loc2, ugi,u_ele, xsgi, usgi,usgi2,income,&
              snorm,norm,t_loc,t_old,t_new,tgi,u_new, u_old,mass_ele,rhs_loc)
end program wave_equation


subroutine coordinates(ndim, totele, x_all, no_ele_row, no_ele_col, row, col, nloc, dx, dy)
  implicit none
  integer :: ele, row, col
  integer, intent(in) :: totele,no_ele_row, no_ele_col, nloc, ndim
  real, intent(in) :: dx, dy
  real :: x_all(ndim,nloc,totele)

  do ele = 1,totele
    row = ceiling(real(ele)/no_ele_row)
    col = ele-(no_ele_row*(row-1))

    x_all(1,1,ele) = dx*(col-1); x_all(2,1,ele) = dy*(row-1)
    x_all(1,2,ele) = dx*col    ; x_all(2,2,ele) = dy*(row-1)
    x_all(1,3,ele) = dx*(col-1); x_all(2,3,ele) = dy*row
    x_all(1,4,ele) = dx*col    ; x_all(2,4,ele) = dy*row
  end do
end subroutine coordinates



subroutine ele_neighbour(iface, totele, nface, ele, face_ele, no_ele_row, row, row2, s_list_no)
  ! ordering the face numbers: bottom face=1, right=1, left=3 and top=4
  ! row and row2 are row number associated with ele and ele2
  ! no_ele_row is total number of element in each row
  ! face_ele(iface, ele) = given the face no iface and element no return  &
  ! the element next to the surface or if negative return 0 incase o ele=1 &
  ! or the negative of the surface element number between element ele and face iface.
  implicit none
  integer, intent(in) :: ele, iface, no_ele_row, totele, nface
  integer, intent(inout) :: row, row2, face_ele, s_list_no

  row = ceiling(real(ele)/no_ele_row)
  if (iface==1) then
    face_ele = ele - no_ele_row
    s_list_no = 4
    row2 = ceiling(real(ele)/no_ele_row)
    if (row2.EQ.1) s_list_no = -1*s_list_no

  elseif ( iface==2 ) then
    face_ele = ele + 1
    s_list_no = 3
    row2 = ceiling(real(face_ele)/no_ele_row)
    if (row2.NE.row) then
      face_ele = -1*face_ele  !It is a boundary element located at the right side of the domain
      s_list_no = -1*s_list_no
    end if

  elseif ( iface==3 ) then
    face_ele = ele -1   !It is a boundary element
    s_list_no = 2
    row2 = ceiling(real(face_ele)/no_ele_row)
    if (row2.NE.row) then
      face_ele = -1*face_ele  !It is a boundary element located at the lest side of the domain
      s_list_no = -1*s_list_no
    end if

  elseif ( iface==4 ) then
    face_ele = ele + no_ele_row
    s_list_no = 1
    if (face_ele.GT.totele) then
      face_ele = -1*face_ele  !It is a boundary element located at the top of the domain
      s_list_no = -1*s_list_no
    end if
  end if
end subroutine ele_neighbour



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
    ! ewrite(3,*)'gi, detj, weight(gi)', gi, detj, weight(gi)
    ! rsum = rsum + detj
    ! rsumabs = rsumabs + abs( detj )
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



SUBROUTINE RE2DN4(LOWQUA,NGI,NLOC,MLOC,M,WEIGHT,N,NLX,NLY,SNGI,SNLOC,SWEIGH,SN,SNLX)
  ! use FLDebug
  IMPLICIT NONE
  ! NB might have to define surface elements for p and (u,v,w)
  ! in here as well.
  ! This subroutine defines the shape functions M and N and their
  ! derivatives at the Gauss points
  ! REAL M(1,NGI),WEIGHT(NGI),N(4,NGI),NLX(4,NGI),NLY(4,NGI)
  INTEGER, intent(in):: NGI,NLOC,MLOC
  REAL:: M(MLOC,NGI),WEIGHT(NGI)
  REAL:: N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI)
  REAL:: POSI,TLY
  REAL:: LX(16),LY(16),LXP(4),LYP(4)
  REAL:: WEIT(16)
  INTEGER:: SNGI,SNLOC
  REAL ::SWEIGH(SNGI)
  REAL:: SN(SNLOC,SNGI),SNLX(SNLOC,SNGI)
  INTEGER:: P,Q,CORN,GPOI,ILOC,JLOC,NDGI
  LOGICAL:: LOWQUA,GETNDP
  INTEGER:: I
  ! NB LXP(I) AND LYP(I) ARE THE LOCAL X AND Y COORDS OF NODAL POINT I

  ! ewrite(3,*)'inside re2dn4, nloc,mloc,ngi',&
                ! nloc,mloc,ngi

  LXP(1)=-1
  LYP(1)=-1

  LXP(2)=1
  LYP(2)=-1

  ! LXP(3)=1
  ! LYP(3)=1

  ! LXP(4)=-1
  ! LYP(4)=1

  LXP(3)=-1
  LYP(3)=1

  LXP(4)=1
  LYP(4)=1

  IF(NGI.EQ.4) THEN
    POSI=1.0/SQRT(3.0)
    LX(1)=-POSI
    LY(1)=-POSI
    LX(2)= POSI
    LY(2)= POSI

    do  Q=1,2! Was loop 23
      do  P=1,2! Was loop 24
        do  CORN=1,4! Was loop 25
          GPOI=(Q-1)*2 + P

          IF(MLOC.EQ.1)  M(1,GPOI)=1.
            WEIGHT(GPOI)=1.

            N(CORN,GPOI)=0.25*(1.+LXP(CORN)*LX(P))&
                        *(1.+LYP(CORN)*LY(Q))
            NLX(CORN,GPOI)=0.25*LXP(CORN)*(1.+LYP(CORN)*LY(Q))
            NLY(CORN,GPOI)=0.25*LYP(CORN)*(1.+LXP(CORN)*LX(P))
        end do ! Was loop 25
      end do ! Was loop 24
    end do ! Was loop 23
    ! ewrite(3,*) 'here 1'
    ! ewrite(3,*) 'N:',N
    ! ewrite(3,*) 'NLX:',NLX
    ! ewrite(3,*) 'NLY:',NLY
    ! Surface shape functions
    IF((SNGI.GT.1).AND.(SNLOC.GT.1)) THEN
       ! ewrite(3,*) '***************** SNGI=',SNGI
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
    ! ewrite(3,*) 'ndgi,ngi,sngi:',ndgi,ngi,sngi

    GETNDP=.FALSE.
    CALL LAGROT(WEIT,LX,NDGI,GETNDP)
    LY(1:NDGI) = LX(1:NDGI)
    ! ewrite(3,*) 'weit:',weit
    ! ewrite(3,*) 'lx:',lx

    do  Q=1,NDGI! Was loop 323
      do  P=1,NDGI! Was loop 324
        do  CORN=1,4! Was loop 325
          ! ewrite(3,*) 'q,p,corn:',q,p,corn
          GPOI=(Q-1)*NDGI + P
          IF(MLOC.EQ.1)  M(1,GPOI)=1.
          WEIGHT(GPOI)=WEIT(P)*WEIT(Q)
          ! ewrite(3,*) 'here1'
          N(CORN,GPOI)=0.25*(1.+LXP(CORN)*LX(P))&
                           *(1.+LYP(CORN)*LY(Q))
          ! ewrite(3,*) 'here2'
          NLX(CORN,GPOI)=0.25*LXP(CORN)*(1.+LYP(CORN)*LY(Q))
          NLY(CORN,GPOI)=0.25*LYP(CORN)*(1.+LXP(CORN)*LX(P))
          ! ewrite(3,*) 'here3'
        end do ! Was loop 325
      end do ! Was loop 324
    end do ! Was loop 323
    ! ewrite(3,*) 'here 1'
    ! ewrite(3,*) 'N:',N
    ! ewrite(3,*) 'NLX:',NLX
    ! ewrite(3,*) 'NLY:',NLY
    ! Surface shape functions
    ! ewrite(3,*) '***************** SNGI=',SNGI
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
  ! END OF IF(NGI.EQ.4) THEN ELSE ...
  ENDIF

  IF(MLOC.EQ.NLOC) THEN
    do  I=1,4! Was loop 2545
      do  CORN=1,4! Was loop 2545
        M(CORN,I)=N(CORN,I)
      end do ! Was loop 2545
    end do ! Was loop 2545
  ENDIF
  ! ewrite(3,*) 'in re2dn4.f here 2 ngi,sngi',ngi,sngi
  ! ewrite(3,*) 'N:',N
  ! ewrite(3,*) 'NLX:',NLX
  ! ewrite(3,*) 'NLY:',NLY
  ! END
END SUBROUTINE RE2DN4





SUBROUTINE LAGROT(WEIT,QUAPOS,NDGI,GETNDP)
  ! use RGPTWE_module
  IMPLICIT NONE
  ! This computes the weight and points for standard Gaussian quadrature.
  ! IF(GETNDP) then get the POSITION OF THE NODES
  ! AND DONT BOTHER WITH THE WEITS.
  INTEGER:: NDGI
  REAL:: WEIT(NDGI),QUAPOS(NDGI)
  LOGICAL:: GETNDP
  LOGICAL:: WEIGHT
  INTEGER ::IG
  !real function...
  real:: RGPTWE

  IF(.NOT.GETNDP) THEN
    WEIGHT=.TRUE.
    do IG=1,NDGI
      WEIT(IG)=RGPTWE(IG,NDGI,WEIGHT)
    END DO

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
  !inv_jac )
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
  !REAL, DIMENSION( NDIM,ndim ), intent( in ) :: inv_jac
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
    ! inv_jac(1,1)=DXDLX; inv_jac(1,2)=DXDLY; inv_jac(1,3)=DXDLZ
    ! inv_jac(2,1)=DyDLX; inv_jac(2,2)=DyDLY; inv_jac(2,3)=DyDLZ
    ! inv_jac(3,1)=DzDLX; inv_jac(3,2)=DzDLY; inv_jac(3,3)=DzDLZ
    ! inv_jac=inv_jac/detj
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
