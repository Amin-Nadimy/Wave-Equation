! 2D DG-FEM wave equation
! face_ele(iface, ele) = given the face no iface and element no return the element next to
! the surface or if negative return the negative of the surface element number between element ele and face iface.

program wave_equation
  implicit none

  logical :: LOWQUA
  integer :: nloc, gi, ngi, sngi, iface, nface, totele, ele, ele2, ele22
  integer :: s_list_no, s_gi, iloc, jloc
  integer:: itime, ntime, idim, ndim, sndim, nonodes, snloc, mloc, col, n_s_list_no, row, row2
  integer :: no_ele_col, no_ele_row, errorflag, max_face_list_no, i_got_boundary, jac_its, njac_its
  integer, allocatable :: face_ele(:,:), face_list_no(:,:)

  real :: sarea, volume, dt, CFL, L, dx, dy, r_got_boundary
  logical :: direct_solver

  real, allocatable :: n(:,:), nlx(:,:,:), nx(:,:,:), M(:,:)
  real, allocatable :: weight(:), detwei(:), sdetwei(:)
  real, allocatable :: sn(:,:),sn2(:,:),snlx(:,:,:),sweigh(:), s_cont(:,:)
  real, allocatable :: t_loc(:), t_loc2(:), t_old(:,:), t_new(:,:), t_new_nonlin(:,:)
  real, allocatable :: t_bc(:,:), u_bc(:,:,:)
  real, allocatable :: tgi(:), tsgi(:), tsgi2(:)
  real, allocatable :: face_sn(:,:,:), face_sn2(:,:,:), face_snlx(:,:,:,:), face_sweigh(:,:)
  real, allocatable :: u_loc(:,:), u_ele(:,:,:), u_loc2(:,:), x_loc(:,:), ugi(:,:),u_old(:)
  real, allocatable :: x_all(:,:,:), jac(:)
  real, allocatable :: xsgi(:,:), usgi(:,:), usgi2(:,:), income(:), snorm(:,:), norm(:)
  real, allocatable :: mass_ele(:,:), mass_ele_inv(:,:), rhs_loc(:), ml_ele(:), rhs_jac(:), mass_t_new(:)
  real, allocatable :: SN_orig(:,:),SNLX_orig(:,:,:)
  ! real :: rdum(:)

  CFL = 0.01
  L = 1.0   ! length of the domain in each direction
  no_ele_row = 100
  no_ele_col = 1
  totele = no_ele_row * no_ele_col
  nface =  4
  !dx = L/(no_ele_row)
  !dy = L/(no_ele_col)
  dx = 1.0
  dy = 1.0
  ndim = 2
  sndim = ndim-1 ! dimensions of the surface element coordinates.
  snloc = 2 ! no of nodes on each boundary line
  ngi = 4
  sngi = 2 ! no of surface quadrature points of the faces - this is set to the max no of all faces
  nloc = 4
  mloc=1
  ntime = 3000
  n_s_list_no = 4
  njac_its=10 ! no of Jacobi iterations
  direct_solver=.false. ! use directi solver?
  ! max_face_list_no = 4

  allocate(n( ngi, nloc ), nx( ngi, ndim, nloc ), nlx( ngi, ndim, nloc ), M(ngi,nloc))
  allocate(weight(ngi), detwei(ngi), sdetwei(sngi))
  allocate(face_list_no( nface, totele), face_ele(nface,totele))
!AMIN l1, l2, l3, l4 aren't used anywhere
  ! allocate(l1(ngi), l2(ngi), l3(ngi), l4(ngi))
  allocate(sn(sngi,nloc),sn2(sngi,nloc),snlx(sngi,sndim,nloc),sweigh(sngi))
  allocate(SN_orig(sngi,snloc),SNLX_orig(sngi,sndim,snloc), s_cont(sngi,ndim) )
  allocate(u_loc(ndim,nloc), u_loc2(ndim,nloc), ugi(ngi,ndim), u_ele(ndim,nloc,totele))
  allocate(xsgi(sngi,ndim), usgi(sngi,ndim), usgi2(sngi,ndim), income(sngi), snorm(sngi,ndim), norm(ndim))
  allocate(t_loc(nloc), t_loc2(nloc), t_old(nloc,totele), t_new(nloc,totele), t_new_nonlin(nloc,totele))
  allocate(t_bc(nloc,totele), u_bc(ndim,nloc,totele))
  allocate(tgi(ngi), tsgi(sngi), tsgi2(sngi))
  allocate(face_sn(sngi,nloc,nface), face_sn2(sngi,nloc,n_s_list_no), face_snlx(sngi,sndim,nloc,nface) )
  allocate(face_sweigh(sngi,nface))
  allocate(x_all(ndim,nloc,totele), x_loc(ndim,nloc), jac(ndim*ndim))
  allocate(mass_ele(nloc,nloc), mass_ele_inv(nloc,nloc), rhs_loc(nloc), ml_ele(nloc))
  allocate(rhs_jac(nloc), mass_t_new(nloc) )
!AMIN, rdum isn't used anywhere
  ! allocate(rdum(10))

!AMIN mdum isn't used anywhere
  ! mdum(1:10)=0.0
  t_bc(:,:)=0.0 ! this contains the boundary conditions just outside the domain
  u_bc(:,:,:)=0.0 ! this contains the boundary conditions on velocity just outside the domain
  ! initial conditions
  t_new(:,:) = 0.0
  t_new(:,no_ele_row/5:no_ele_row/2) = 1.0 ! a bit off - is this correct -only correct for 1D?
  U_ele(:,:,:) = 0
  U_ele(1,:,1:totele) = 1.0 ! suggest this should be unity
  u_bc(:,:,:)=u_ele(:,:,:) ! this contains the boundary conditions on velocity just outside the domain
!  dt = CFL/((u_ele(1,1,1)/dx)+(u_ele(2,1,1)/dy))
  dt = CFL*dx
  print *,'dt=',dt

  LOWQUA=.false.
  call RE2DN4(LOWQUA,NGI,NLOC,MLOC,M,WEIGHT,N,NLX(:,1,:),NLX(:,2,:),  SNGI,SNLOC,SWEIGH,SN_orig,SNLX_orig)
  call ele_info(totele, nface, face_ele, no_ele_row, row, row2, &
                x_all, dx, dy, ndim, nloc, no_ele_col, col)
  call surface_pointers_sn(nface, sngi, snloc, nloc, ndim, sndim, totele, n_s_list_no, no_ele_row, no_ele_col, &
               sn_orig, snlx_orig, face_sn, face_snlx, face_sn2, face_list_no)
! print*, face_snlx(:,:,4,1) ! it gives results for iloc 1,3 and 4 but not 2.

  do itime=1,ntime
    t_old =t_new ! prepare for next time step
    t_new_nonlin = t_new ! for non-linearity
    do ele=1,totele
      ! volume integration
      x_loc(:,:)=x_all(:,:,ele) ! x_all contains the coordinates of the corner nodes
      call det_nlx( x_loc, n, nlx, nx, detwei, weight, ndim, nloc, ngi )
      !volume=sum(detwei)
!      do gi=1,ngi
!         print *,'gi=',gi
!         print *,'nx(gi,1,:):',nx(gi,1,:)
!         print *,'nx(gi,2,:):',nx(gi,2,:)
!      end do
!      stop 221

      u_loc(:,:) = u_ele(:,:,ele) ! this is the
!AMIN changed to t_old
      t_loc(:) =  t_old(:,ele) ! this is the FEM element - 2nd index is the local node number

      ! calculates vel at gi sum = phi_1(gi)*c1 +phi_2(gi)*c2 +phi_3(gi)*c3 +phi_4(gi)*c4
      ! do idim=1,ndim
      !   do gi=1,ngi
      !      !ugi_x(gi,idim)=sum(nx(gi,idim,:)*u_loc(idim,:))
      !      ugi(gi,idim)=sum(n(gi,:)*u_loc(idim,:))
      !   end do
      ! end do
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
          end do ! quadrature
        end do ! ijloc
!
        ml_ele(iloc)=sum(n(:,iloc)*detwei(:)) ! lumped mass matrix in element (use later for solver).
        do gi=1,ngi
          do idim=1,ndim
            rhs_loc(iloc) = rhs_loc(iloc) + nx(gi,idim,iloc)*ugi(gi,idim)*tgi(gi)*detwei(gi)
          end do
        end do ! quadrature
      end do ! iloc
!       print *,'ele,tgi(:):',ele,tgi(:)
!       print *,'t_loc(:):',t_loc(:)
!       print *,'ugi(:,:):',ugi(:,:)

      !Include the surface integral here:
      do iface = 1,nface
        ele22 = face_ele( iface, ele) ! if ele22 is negative assume on boundary
        i_got_boundary = (sign(1, -ele22) + 1 )/2
        r_got_boundary = real(i_got_boundary)
        ele2 = ele*i_got_boundary + ele22*(1-i_got_boundary)
! r_got_boundary=1.0 if we want to use the boundary conditions and have incomming velocity.
! r_got_boundary=0.0 not on boundary

        t_loc2(:)=t_new(:,ele2) * (1.0-r_got_boundary)     + t_bc(:,ele)  * r_got_boundary
        u_loc2(:,:)=u_ele(:,:,ele2)* (1.0-r_got_boundary)  + u_bc(:,:,ele)* r_got_boundary

        !Surface integral along an element
        ! need to work on this also
        s_list_no = face_list_no( iface, ele) ! correct
        sn(:,:)     = face_sn(:,:,iface)  ! correct
        snlx(:,:,:) = face_snlx(:,:,:,iface)
        sn2(:,:)    = face_sn2(:,:,s_list_no) ! correct

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
           norm(idim) = sum(xsgi(:,idim))/real(sngi) - sum(x_loc(idim,:))/real(nloc)
        end do
        ! start to do the integration
        call det_snlx_all( nloc, sngi, sndim, ndim, x_loc, sn, snlx, sweigh, sdetwei, sarea, snorm, norm )
!      print *,'iface,sum(sdetwei(:)):',iface,sum(sdetwei(:))

! income=1 if info is comming from neighbouring element.
        do s_gi=1,sngi
          income(s_gi)=0.5 + 0.5*sign(1.0, -sum(snorm(s_gi,:)*0.5*(usgi(s_gi,:)+usgi2(s_gi,:)))  )
        end do
!        print *,'iface,income:',iface,income
!        print *,'snorm(1,:):',snorm(1,:)
!        print *,'snorm(2,:):',snorm(2,:)
!        print *,'ele22, i_got_boundary, r_got_boundary, ele, ele2:',ele22, i_got_boundary, r_got_boundary, ele, ele2

        ! sloc_vec=0.0; sloc_vec2=0.0
        do idim=1,ndim
          s_cont(:,idim) = snorm(:,idim)*sdetwei(:) &
                      *( (1.-income(:))* usgi(:,idim)*tsgi(:) + income(:)*usgi2(:,idim)*tsgi2(:) )
        end do

        do iloc=1,nloc
          do idim=1,ndim
            rhs_loc(iloc)  = rhs_loc(iloc)  - sum( sn(:,iloc)*s_cont(:,idim) )
          end do
        end do
      end do ! iface

      if(direct_solver) then
      ! inverse of the mass matric (nloc,nloc)
         call FINDInv(mass_ele, mass_ele_inv, nloc, errorflag)
         do iloc=1,nloc
            t_new_nonlin(iloc,ele)=t_old(iloc,ele) + dt*sum( mass_ele_inv(iloc,:) * rhs_loc(:) )
         end do
      else ! iterative solver
         do iloc=1,nloc
            rhs_jac(iloc)= sum( mass_ele(iloc,:) * t_old(:,ele) ) + dt*rhs_loc(iloc)
         end do

         do jac_its=1,njac_its ! jacobi iterations...
            do iloc=1,nloc
               mass_t_new(iloc)= sum( mass_ele(iloc,:) * t_new_nonlin(:,ele) )
            end do
            do iloc=1,nloc
               t_new_nonlin(iloc,ele)=  (ml_ele(iloc)*t_new_nonlin(iloc,ele) - mass_t_new(iloc) + rhs_jac(iloc) ) / ml_ele(iloc)
            end do
         end do
      endif ! endof if(direct_solver) then else

      ! print*, ele, t_new(:,ele)
    end do ! do ele=1,totele
    t_new=t_new_nonlin
  end do ! do itime=1,ntime

!       print *,'t_new:',t_new
!       print *,'t_old:',t_old
   print *,' '
   do ele=1,totele
      print *,dx*real(ele-1),t_new(1,ele)
      print *,dx*real(ele),  t_new(2,ele)
   end do
   print *,' '
   do ele=1,totele
!      print *,dx*real(ele-1),t_new(1,ele)
!      print *,dx*real(ele),  t_new(2,ele)
      print *,dx*real(ele-1),t_old(1,ele)
      print *,dx*real(ele),  t_old(2,ele)
   end do

  deallocate(n, nx, nlx, M, weight, detwei, sdetwei, face_list_no, face_ele, &
             sn,sn2,snlx,sweigh, SN_orig,SNLX_orig, s_cont, u_loc, u_loc2, ugi, u_ele, &
             xsgi, usgi, usgi2, income, snorm, norm, t_loc, t_loc2, t_old, t_new, &
             tgi, tsgi, tsgi2, face_sn, face_sn2, face_snlx, x_all, x_loc, &
             mass_ele, mass_ele_inv, rhs_loc)

end program wave_equation




subroutine surface_pointers_sn(nface, sngi, snloc, nloc, ndim, sndim, totele, n_s_list_no, nele_x, nele_y, &
               sn_orig, snlx_orig, face_sn, face_snlx, face_sn2, face_list_no)
  ! ********************************************************************************************************
  ! calculate face_list_no( iface, ele), face_sn(:,:,iface); face_snlx(:,:,:,iface), face_sn2(:,:,s_list_no)
  ! ********************************************************************************************************
  ! ****integers:
  ! nface=no of faces of each elemenet.
  ! sngi=no of surface quadrature points of the faces - this is set to the max no of all faces.
  ! snloc=no of nodes on a surface element.
  ! nloc=no of nodes within a volume element.
  ! ndim=no of dimensions - including time possibly.
  ! sndim=ndim-1 dimensions of the surface elements.
  ! totele=no of elements.
  ! n_s_list_no= no of different oriantations for the surface element numbering.
  ! nele_x = no of elements across in the x-direction.
  ! nele_y = no of elements across in the y-direction.
  !
  ! ****original surface element shape functions:
  ! sn_orig(sgi,siloc) = shape function at quadrature point sgi and surface node siloc
  ! snlx_orig(sgi,1,siloc) = shape function derivative (only one in 2D) at quadrature point sgi and surface node siloc
  !
  ! ****new surface element shape functions designed to be highly efficient:
  ! face_sn(sgi,iloc,iface) = shape function at quadrature point sgi and volume node iloc for face iface of element.
  ! face_snlx(sgi,1,iloc,iface) = shape function derivative at quadrature point sgi and volume node iloc for face iface of element.
  ! face_sn2(sgi,iloc,n_s_list_no) = shape function at quadrature point sgi but on the otherside of face and volume node iloc for face iface of element.
  ! This works from: sn2 =face_sn2(:,:,s_list_no) in which  s_list_no = face_list_no( iface, ele) .
  implicit none
  integer, intent(in) :: nface, sngi, snloc, nloc, ndim, sndim, totele, n_s_list_no, nele_x, nele_y
  real, intent(in) :: sn_orig(sngi,snloc), snlx_orig(sngi,sndim,snloc)
  real, intent(out) :: face_sn(sngi,nloc,nface), face_snlx(sngi,sndim,nloc,nface), face_sn2(sngi,nloc,n_s_list_no)
  integer, intent(out) :: face_list_no(nface,totele)
  ! local variables...
  integer iface,lnod1,lnod2,i_s_list_no

  face_sn(:,:,:)=0.0
  face_snlx(:,:,:,:)=0.0
  face_sn2(:,:,:)=0.0

  ! local node numbers:
  !  3   4
  !  1   2
  !
  ! face numbers:
  !    4
  !  2   3
  !    1

  iface=1
  lnod1=2
  lnod2=1
  face_sn(:,lnod1,iface)=sn_orig(:,1)
  face_sn(:,lnod2,iface)=sn_orig(:,2)
  face_snlx(:,1,lnod1,iface)=snlx_orig(:,1,1)
  face_snlx(:,1,lnod2,iface)=snlx_orig(:,1,2)
  i_s_list_no = iface
  lnod1=4
  lnod2=3
  face_sn2(:,lnod1,i_s_list_no)=sn_orig(:,1)
  face_sn2(:,lnod2,i_s_list_no)=sn_orig(:,2)
  face_list_no(iface,:) = i_s_list_no

  iface=2
  lnod1=1
  lnod2=3
  face_sn(:,lnod1,iface)=sn_orig(:,1)
  face_sn(:,lnod2,iface)=sn_orig(:,2)
  face_snlx(:,1,lnod1,iface)=snlx_orig(:,1,1)
  face_snlx(:,1,lnod2,iface)=snlx_orig(:,1,2)
  i_s_list_no = iface
  lnod1=2
  lnod2=4
  face_sn2(:,lnod1,i_s_list_no)=sn_orig(:,1)
  face_sn2(:,lnod2,i_s_list_no)=sn_orig(:,2)
  face_list_no(iface,:) = i_s_list_no

  iface=3
  lnod1=4
  lnod2=2
  face_sn(:,lnod1,iface)=sn_orig(:,1)
  face_sn(:,lnod2,iface)=sn_orig(:,2)
  face_snlx(:,1,lnod1,iface)=snlx_orig(:,1,1)
  face_snlx(:,1,lnod2,iface)=snlx_orig(:,1,2)
  i_s_list_no = iface
  lnod1=3
  lnod2=1
  face_sn2(:,lnod1,i_s_list_no)=sn_orig(:,1)
  face_sn2(:,lnod2,i_s_list_no)=sn_orig(:,2)
  face_list_no(iface,:) = i_s_list_no

  iface=4
  lnod1=3
  lnod2=4
  face_sn(:,lnod1,iface)=sn_orig(:,1)
  face_sn(:,lnod2,iface)=sn_orig(:,2)
  face_snlx(:,1,lnod1,iface)=snlx_orig(:,1,1)
  face_snlx(:,1,lnod2,iface)=snlx_orig(:,1,2)
  i_s_list_no = iface
  lnod1=1
  lnod2=2
  face_sn2(:,lnod1,i_s_list_no)=sn_orig(:,1)
  face_sn2(:,lnod2,i_s_list_no)=sn_orig(:,2)
  face_list_no(iface,:) = i_s_list_no
  !
  ! calculate face_list_no(nface,totele) =
  !  do jele=1,nele_y
  !    do iele=1,nele_x
  !       ele=(jele-1)*nele_x + iele
  !    end do
  !  end do

end subroutine surface_pointers_sn



subroutine ele_info(totele, nface, face_ele, no_ele_row, row, row2, &
                         x_all, dx, dy, ndim, nloc, no_ele_col, col)
  ! ordering the face numbers: bottom face=1, right=1, left=3 and top=4
  ! row and row2 are row number associated with ele and ele2
  ! no_ele_row is total number of element in each row
  ! no_ele_col is total number of element in each column
  ! ndim is no of dimensions
  ! nloc is no of corner nodes
  ! face_ele(iface, ele) = given the face no iface and element no return the element next to
  ! the surface or if negative return the negative of the surface element number between element ele and face iface.

  implicit none
  integer, intent(in) :: no_ele_row, totele, no_ele_col, nloc, ndim, nface
  integer, intent(inout) :: row, row2 , col ! why are these inout????
  integer:: face_ele(nface,totele), face_list_no(nface,totele), ele, iface

  real, intent(in) :: dx, dy
  real, intent(out) :: x_all(ndim,nloc,totele)

  do ele=1,totele
    row = ceiling(real(ele)/no_ele_row)
    col = ele-(no_ele_row*(row-1))

    ! corner node coordinates
    x_all(1,1,ele) = dx*(col-1); x_all(2,1,ele) = dy*(row-1)
    x_all(1,2,ele) = dx*col    ; x_all(2,2,ele) = dy*(row-1)
    x_all(1,3,ele) = dx*(col-1); x_all(2,3,ele) = dy*row
    x_all(1,4,ele) = dx*col    ; x_all(2,4,ele) = dy*row

    ! findin neighbiuring ele and face numbers
    do iface=1,nface
      if (iface==1) then
        face_ele(iface,ele) = ele - no_ele_row
        ! face_list_no(iface,ele) = 4
        ! row2 = ceiling(real(ele)/no_ele_row)
        ! if (row2.EQ.1) face_list_no(iface,ele) = -1*face_list_no(iface,ele)

      elseif ( iface==2 ) then
        face_ele(iface,ele) = ele - 1
        ! face_list_no(iface,ele) = 3
        row2 = ceiling(real(face_ele(iface,ele))/no_ele_row)
        if (row2.NE.row) then
          face_ele(iface,ele) = -1*face_ele(iface,ele)  !It is a boundary element located at the right side of the domain
          ! face_list_no(iface,ele) = -1*face_list_no(iface,ele)
        end if

      elseif ( iface==3 ) then
        face_ele(iface,ele) = ele +1   !It is a boundary element
        ! face_list_no(iface,ele) = 2
        row2 = ceiling(real(face_ele(iface,ele))/no_ele_row)
        if (row2.NE.row) then
          face_ele(iface,ele) = -1*face_ele(iface,ele)  !It is a boundary element located at the lest side of the domain
          ! face_list_no(iface,ele) = -1*face_list_no(iface,ele)
        end if

      elseif ( iface==4 ) then
        face_ele(iface,ele) = ele + no_ele_row
        ! face_list_no(iface,ele) = 1
        if (face_ele(iface,ele).GT.totele) then
          face_ele(iface,ele) = -1*face_ele(iface,ele)  !It is a boundary element located at the top of the domain
          ! face_list_no(iface,ele) = -1*face_list_no(iface,ele)
        end if
      end if
    end do
  end do
end subroutine ele_info



!subroutine det_nlx( x_loc, n, nlx, nx, detwei, weight, ndim, nloc, ngi, jac )
subroutine det_nlx( x_loc, n, nlx, nx, detwei, weight, ndim, nloc, ngi )
  ! ****************************************************
  ! This sub form the derivatives of the shape functions
  ! ****************************************************
  ! x_loc: spatial nodes.
  ! n, nlx, nlx_lxx: shape function and local derivatives of the shape functions (nlx_lxx is local grad of the local laplacian)- defined by shape functional library.
  ! nx, nx_lxx: derivatives of the shape functions.
  ! detwei, inv_jac: determinant at the quadrature pots and inverse of Jacobian at quadrature pts.
  ! ndim,nloc,ngi: no of dimensions, no of local nodes within an element, no of quadrature points.
  ! nlx_nod, nx_nod: same as nlx and nx but formed at the nodes not quadrature points.
  implicit none
  integer, intent( in ) :: ndim,nloc,ngi

  REAL, DIMENSION( ndim,nloc ), intent( in ) :: x_loc
  REAL, DIMENSION( ngi, nloc ), intent( in ) :: N
  REAL, DIMENSION( ngi, ndim, nloc ), intent( in ) :: nlx
  REAL, DIMENSION( ngi ), intent( in ) :: WEIGHT
  REAL, DIMENSION( ngi ), intent( inout ) :: DETWEI
  REAL, DIMENSION( ngi, ndim, nloc ), intent( inout ) :: nx
!  real, DIMENSION( ndim*ndim ) :: jac
  ! Local variables
  REAL :: AGI, BGI, CGI, DGI, EGI, FGI, GGI, HGI, KGI, A11, A12, A13, A21, &
          A22, A23, A31, A32, A33, DETJ
  INTEGER :: GI, L, IGLX, ii

  if (ndim==2) then
   ! conventional:
    do  GI=1,NGI! Was loop 331

      AGI=0.
      BGI=0.

      CGI=0.
      DGI=0.

      do  L=1,NLOC! Was loop 79
        IGLX=L

        AGI=AGI+NLX(GI,1,L)*x_loc(1,L)
        BGI=BGI+NLX(GI,1,L)*x_loc(2,L)

        CGI=CGI+NLX(GI,2,L)*x_loc(1,L)
        DGI=DGI+NLX(GI,2,L)*x_loc(2,L)

      end do ! Was loop 79

      DETJ= AGI*DGI-BGI*CGI
      DETWEI(GI)=ABS(DETJ)*WEIGHT(GI)
      ! For coefficient in the inverse mat of the jacobian.
      A11= DGI /DETJ
      A21=-BGI /DETJ

      A12=-CGI /DETJ
      A22= AGI /DETJ

      do  L=1,NLOC! Was loop 373
        NX(GI,1,L)= A11*NLX(GI,1,L)+A12*NLX(GI,2,L)
        NX(GI,2,L)= A21*NLX(GI,1,L)+A22*NLX(GI,2,L)
      end do ! Was loop 373
    end do ! GI Was loop 331

!    jac(1) = AGI; jac(2) = DGI ; jac(3) = BGI ; jac(4) = EGI

  elseif ( ndim.eq.3 ) then
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
        AGI=AGI+NLX(GI,1,L)*x_loc(1,IGLX)
        BGI=BGI+NLX(GI,1,L)*x_loc(2,IGLX)
        CGI=CGI+NLX(GI,1,L)*x_loc(3,IGLX)

        DGI=DGI+NLX(GI,2,L)*x_loc(1,IGLX)
        EGI=EGI+NLX(GI,2,L)*x_loc(2,IGLX)
        FGI=FGI+NLX(GI,2,L)*x_loc(3,IGLX)

        GGI=GGI+NLX(GI,3,L)*x_loc(1,IGLX)
        HGI=HGI+NLX(GI,3,L)*x_loc(2,IGLX)
        KGI=KGI+NLX(GI,3,L)*x_loc(3,IGLX)
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
        NX(GI,1,L)= A11*NLX(GI,1,L)+A12*NLX(GI,2,L)+A13*NLX(GI,3,L)
        NX(GI,2,L)= A21*NLX(GI,1,L)+A22*NLX(GI,2,L)+A23*NLX(GI,3,L)
        NX(GI,3,L)= A31*NLX(GI,1,L)+A32*NLX(GI,2,L)+A33*NLX(GI,3,L)
      end do ! Was loop 373
    end do ! GI Was loop 331
  end if
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
  REAL:: M(NGI,MLOC),WEIGHT(NGI)
  REAL:: N(NGI,NLOC),NLX(NGI,NLOC),NLY(NGI,NLOC)
  REAL:: POSI,TLY
  REAL:: LX(16),LY(16),LXP(4),LYP(4)
  REAL:: WEIT(16)
  INTEGER:: SNGI,SNLOC
  REAL ::SWEIGH(SNGI)
  REAL:: SN(SNGI,SNLOC),SNLX(SNGI,SNLOC)
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

          IF(MLOC.EQ.1)  M(GPOI,1)=1.
            WEIGHT(GPOI)=1.

            N(GPOI,CORN)=0.25*(1.+LXP(CORN)*LX(P))&
                        *(1.+LYP(CORN)*LY(Q))
            NLX(GPOI,CORN)=0.25*LXP(CORN)*(1.+LYP(CORN)*LY(Q))
            NLY(GPOI,CORN)=0.25*LYP(CORN)*(1.+LXP(CORN)*LX(P))
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
                     SN(GPOI,CORN)=0.5*(1.+LXP(CORN)*LX(P))
                     SNLX(GPOI,CORN)=0.5*LXP(CORN)
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
          IF(MLOC.EQ.1)  M(GPOI,1)=1.
          WEIGHT(GPOI)=WEIT(P)*WEIT(Q)
          ! ewrite(3,*) 'here1'
          N(GPOI,CORN)=0.25*(1.+LXP(CORN)*LX(P))&
                           *(1.+LYP(CORN)*LY(Q))
          ! ewrite(3,*) 'here2'
          NLX(GPOI,CORN)=0.25*LXP(CORN)*(1.+LYP(CORN)*LY(Q))
          NLY(GPOI,CORN)=0.25*LYP(CORN)*(1.+LXP(CORN)*LX(P))
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
          SN(GPOI,CORN)=0.5*(1.+LXP(CORN)*LX(P))
          SNLX(GPOI,CORN)=0.5*LXP(CORN)
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
        M(I,CORN)=N(I,CORN)
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
  real :: RGPTWE
  !real, allocatable:: RGPTWE(:,:,:)
  !allocate(RGPTWE(ndgi,1,1))

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
  REAL, DIMENSION( sngi, ndim ), intent( inout ) :: NORMXN_ALL
  REAL, DIMENSION( ndim ), intent( in ) :: NORMX_ALL
  !REAL, DIMENSION( NDIM,ndim ), intent( in ) :: inv_jac
  ! Local variables
  INTEGER :: GI, SL, IGLX
  REAL :: DXDLX, DXDLY, DYDLX, DYDLY, DZDLX, DZDLY
  REAL :: A, B, C, DETJ, RUB3=0.0!, RUB4

  SAREA=0.

  if (ndim==2) then
    DO GI=1,SNGI

      DXDLX=0.
      DXDLY=0.
      DYDLX=0.
      DYDLY=0.
      DZDLX=0.
      DZDLY=0.

      DO SL=1,SNLOC
        DXDLX=DXDLX + SNLX(GI,1,SL)*XSL_ALL(1,SL)
        ! DXDLY=DXDLY + SNLX(GI,2,SL)*XSL_ALL(1,SL)
        DYDLX=DYDLX + SNLX(GI,1,SL)*XSL_ALL(2,SL)
        ! DYDLY=DYDLY + SNLX(GI,2,SL)*XSL_ALL(2,SL)
        ! DZDLX=DZDLX + SNLX(GI,1,SL)*XSL_ALL(3,SL)
        ! DZDLY=DZDLY + SNLX(GI,2,SL)*XSL_ALL(3,SL)
      END DO

      A = DYDLX!*DZDLY - DYDLY*DZDLX
      B = DXDLX!*DZDLY - DXDLY*DZDLX
      ! C = DXDLX*DYDLY - DXDLY*DYDLX

      DETJ=SQRT( A**2 + B**2)! + C**2)
      ! inv_jac(1,1)=DXDLX; inv_jac(1,2)=DXDLY; inv_jac(1,3)=DXDLZ
      ! inv_jac(2,1)=DyDLX; inv_jac(2,2)=DyDLY; inv_jac(2,3)=DyDLZ
      ! inv_jac(3,1)=DzDLX; inv_jac(3,2)=DzDLY; inv_jac(3,3)=DzDLZ
      ! inv_jac=inv_jac/detj
      SDETWE(GI)=DETJ*SWEIGH(GI)
      SAREA=SAREA+SDETWE(GI)

      ! Calculate the normal at the Gauss pts...
      ! Perform x-product. N=T1 x T2
      CALL NORMGI(NORMXN_ALL(GI,1),NORMXN_ALL(GI,2),rub3, &
                  DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,1.0, &
                  NORMX_ALL(1),NORMX_ALL(2),1.0)
    END DO

  elseif ( ndim==3 ) then
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
  end if

  RETURN

END SUBROUTINE det_snlx_all



REAL FUNCTION RGPTWE(IG,ND,WEIGHT)
    IMPLICIT NONE
    !     NB If WEIGHT is TRUE in function RGPTWE then return the Gauss-pt weight
    !     else return the Gauss-pt.
    !     NB there are ND Gauss points we are looking for either the
    !     weight or the x-coord of the IG'th Gauss point.
    INTEGER IG,ND
    LOGICAL WEIGHT

      IF(WEIGHT) THEN
         GO TO (10,20,30,40,50,60,70,80,90,100) ND
         !     +++++++++++++++++++++++++++++++
         !     For N=1 +++++++++++++++++++++++
  10     CONTINUE
         RGPTWE=2.0
         GO TO 1000
         !     For N=2 +++++++++++++++++++++++
  20     CONTINUE
         RGPTWE=1.0
         GO TO 1000
         ! For N=3 +++++++++++++++++++++++
  30     CONTINUE
         GO TO (11,12,11) IG
  11     RGPTWE= 0.555555555555556
         GO TO 1000
  12     RGPTWE= 0.888888888888889
         GO TO 1000
         ! For N=4 +++++++++++++++++++++++
  40     CONTINUE
         GO TO (21,22,22,21) IG
  21     RGPTWE= 0.347854845137454
         GO TO 1000
  22     RGPTWE= 0.652145154862546
         GO TO 1000
         ! For N=5 +++++++++++++++++++++++
  50     CONTINUE
         GO TO (31,32,33,32,31) IG
  31     RGPTWE= 0.236926885056189
         GO TO 1000
  32     RGPTWE= 0.478628670499366
         GO TO 1000
  33     RGPTWE= 0.568888888888889
         GO TO 1000
         ! For N=6 +++++++++++++++++++++++
  60     CONTINUE
         GO TO (41,42,43,43,42,41) IG
  41     RGPTWE= 0.171324492379170
         GO TO 1000
  42     RGPTWE= 0.360761573048139
         GO TO 1000
  43     RGPTWE= 0.467913934572691
         GO TO 1000
         ! For N=7 +++++++++++++++++++++++
  70     CONTINUE
         GO TO (51,52,53,54,53,52,51) IG
  51     RGPTWE= 0.129484966168870
         GO TO 1000
  52     RGPTWE= 0.279705391489277
         GO TO 1000
  53     RGPTWE= 0.381830050505119
         GO TO 1000
  54     RGPTWE= 0.417959183673469
         GO TO 1000
         ! For N=8 +++++++++++++++++++++++
  80     CONTINUE
         GO TO (61,62,63,64,64,63,62,61) IG
  61     RGPTWE= 0.101228536290376
         GO TO 1000
  62     RGPTWE= 0.222381034453374
         GO TO 1000
  63     RGPTWE= 0.313706645877877
         GO TO 1000
  64     RGPTWE= 0.362683783378362
         GO TO 1000
         ! For N=9 +++++++++++++++++++++++
  90     CONTINUE
         GO TO (71,72,73,74,75,74,73,72,71) IG
  71     RGPTWE= 0.081274388361574
         GO TO 1000
  72     RGPTWE= 0.180648160694857
         GO TO 1000
  73     RGPTWE= 0.260610696402935
         GO TO 1000
  74     RGPTWE= 0.312347077040003
         GO TO 1000
  75     RGPTWE= 0.330239355001260
         GO TO 1000
         ! For N=10 +++++++++++++++++++++++
  100    CONTINUE
         GO TO (81,82,83,84,85,85,84,83,82,81) IG
  81     RGPTWE= 0.066671344308688
         GO TO 1000
  82     RGPTWE= 0.149451349150581
         GO TO 1000
  83     RGPTWE= 0.219086362515982
         GO TO 1000
  84     RGPTWE= 0.269266719309996
         GO TO 1000
  85     RGPTWE= 0.295524224714753
         !
  1000   CONTINUE
      ELSE
         GO TO (210,220,230,240,250,260,270,280,290,200) ND
         ! +++++++++++++++++++++++++++++++
         ! For N=1 +++++++++++++++++++++++ THE GAUSS POINTS...
  210    CONTINUE
         RGPTWE=0.0
         GO TO 2000
         ! For N=2 +++++++++++++++++++++++
  220    CONTINUE
         RGPTWE= 0.577350269189626
         GO TO 2000
         ! For N=3 +++++++++++++++++++++++
  230    CONTINUE
         GO TO (211,212,211) IG
  211    RGPTWE= 0.774596669241483
         GO TO 2000
  212    RGPTWE= 0.0
         GO TO 2000
         ! For N=4 +++++++++++++++++++++++
  240    CONTINUE
         GO TO (221,222,222,221) IG
  221    RGPTWE= 0.861136311594953
         GO TO 2000
  222    RGPTWE= 0.339981043584856
         GO TO 2000
         ! For N=5 +++++++++++++++++++++++
  250    CONTINUE
         GO TO (231,232,233,232,231) IG
  231    RGPTWE= 0.906179845938664
         GO TO 2000
  232    RGPTWE= 0.538469310105683
         GO TO 2000
  233    RGPTWE= 0.0
         GO TO 2000
         ! For N=6 +++++++++++++++++++++++
  260    CONTINUE
         GO TO (241,242,243,243,242,241) IG
  241    RGPTWE= 0.932469514203152
         GO TO 2000
  242    RGPTWE= 0.661209386466265
         GO TO 2000
  243    RGPTWE= 0.238619186083197
         GO TO 2000
         ! For N=7 +++++++++++++++++++++++
  270    CONTINUE
         GO TO (251,252,253,254,253,252,251) IG
  251    RGPTWE= 0.949107912342759
         GO TO 2000
  252    RGPTWE= 0.741531185599394
         GO TO 2000
  253    RGPTWE= 0.405845151377397
         GO TO 2000
  254    RGPTWE= 0.0
         GO TO 2000
         ! For N=8 +++++++++++++++++++++++
  280    CONTINUE
         GO TO (261,262,263,264,264,263,262,261) IG
  261    RGPTWE= 0.960289856497536
         GO TO 2000
  262    RGPTWE= 0.796666477413627
         GO TO 2000
  263    RGPTWE= 0.525532409916329
         GO TO 2000
  264    RGPTWE= 0.183434642495650
         GO TO 2000
         ! For N=9 +++++++++++++++++++++++
  290    CONTINUE
         GO TO (271,272,273,274,275,274,273,272,271) IG
  271    RGPTWE= 0.968160239507626
         GO TO 2000
  272    RGPTWE= 0.836031107326636
         GO TO 2000
  273    RGPTWE= 0.613371432700590
         GO TO 2000
  274    RGPTWE= 0.324253423403809
         GO TO 2000
  275    RGPTWE= 0.0
         GO TO 2000
         ! For N=10 +++++++++++++++++++++++
  200    CONTINUE
         GO TO (281,282,283,284,285,285,284,283,282,281) IG
  281    RGPTWE= 0.973906528517172
         GO TO 2000
  282    RGPTWE= 0.865063366688985
         GO TO 2000
  283    RGPTWE= 0.679409568299024
         GO TO 2000
  284    RGPTWE= 0.433395394129247
         GO TO 2000
  285    RGPTWE= 0.148874338981631
         !
  2000   CONTINUE
         IF(IG.LE.INT((ND/2)+0.1)) RGPTWE=-RGPTWE
      ENDIF
END FUNCTION RGPTWE



! normal at GI
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


! cross product
SUBROUTINE XPROD1( AX, AY, AZ, &
                   BX, BY, BZ, &
                   CX, CY, CZ )
  implicit none
  REAL, intent( inout ) :: AX, AY, AZ
  REAL, intent( in )    :: BX, BY, BZ, CX, CY, CZ

  ! Perform cross product(x-product). a=b x c
  AX =    BY * CZ - BZ * CY
  AY = -( BX * CZ - BZ * CX )
  AZ =    BX * CY - BY * CX

  RETURN
END subroutine XPROD1


! Finds inverse of a matrix
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
