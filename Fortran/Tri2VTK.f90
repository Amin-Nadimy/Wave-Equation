module Tri2VTK
	use Lib_VTK_IO
	use structures
	use TriangleOps, only: printMatrix, getUpTriangles, Lex2TriCord_SD, getVoronoiCenter
	use strings
	use generic
	use MeshOps, only:getCollideVertexes, getEqvCoord
    implicit none

	private
	public :: Tri2VTK_SD,Tri2VTK_CC, Tri2VTK_MeshCC, Tri2VTK_MeshSD, Tri2VTK_Mix
	public :: Tri2VTK_MeshMix, Tri2VTK_Sys,Tri2VTK_Sys_Mesh

contains

	!This subroutine creates an VTK archive to be read by ParaView
	!Only works for Mix nodes
	!Tri(in) :: type(triangle)
	!(Optional) var(in) :: Integer. Select the variable to be plotted
	!								1 -> U(default)
	!								2 -> F
	!								3 -> Residual
	!								4 -> Auxiliar
	!(Optional) relief(in) :: Logical. Choose whether plot in 3d(default) or 2d
	!(Optional) filename(in) :: Character(*). Name of the archive to save into
	!(Optional) T(in) :: Integer. Number to add to the filename. Useful to time dependant problems
	!					or to see diferent iterations. As paraview permits to move along time.
	!(Optional) title(in) :: Character(*). Title to show in paraview
	!(Optional) Sys (in) :: integer. If 1 then 2 x 2 system of scalar
	subroutine Tri2VTK_Sys(Tri, var, relief, filename, T, title, Sys)
		!Global variables
		type(Triangle), intent(in) :: Tri
		logical, optional, intent(in) :: relief
		integer, optional, intent(in) :: var, T, Sys
		character (len = *), optional, intent(in) :: filename, title
		!Local variables
		integer :: sis, var2
		logical :: relieve
		character (len = 200) :: filename2, title2
		!Prepare data
		sis = size(Tri%UpUF,3)/2

		if (present(Sys)) Sis = Sys
		!*****************Optional area*******************************
		filename2 = "Triangle_Sys"
		if (present(filename)) filename2 = filename
		title2 = "Semi-structured triangle"
		if (present(title)) title2 = title
		relieve = .false.
		if (present(relief)) relieve = relief
		var2 = 1
		if (present(var)) var2 = var
		!Call methods
		select case (sis)
		    case (1)
		        if (present(t)) then
		            call Tri2VTK_CC(Tri,var=var2,relief=relieve,filename=trim(filename2),title=title2,T=t,Sys=1)
		        else
		        	call Tri2VTK_CC(Tri,var=var2,relief=relieve,filename=trim(filename2),title=title2,Sys=1)
		        end if
		    case (2)
		        if (present(t)) then
		            call Tri2VTK_Mix(Tri, var2, relieve, trim(filename2), T, title2)
		        else
		        	call Tri2VTK_Mix(Tri,var=var2,relief=relieve,filename=trim(filename2),title=title2)
		        end if
		end select

	end subroutine Tri2VTK_Sys


	!This subroutine creates an VTK archive to be read by ParaView
	!Only works for side centered nodes
	!meshL(in) :: type(mesh)
	!(Optional) var(in) :: Integer. Select the variable to be plotted
	!								1 -> U(default)
	!								2 -> F
	!								3 -> Residual
	!								4 -> Auxiliar
	!(Optional) relief(in) :: Logical. Choose whether plot in 3d(default) or 2d
	!(Optional) filename(in) :: Character(*). Name of the archive to save into
	!(Optional) T(in) :: Integer. Number to add to the filename. Useful to time dependant problems
	!					or to see diferent iterations. As paraview permits to move along time.
	!(Optional) title(in) :: Character(*). Title to show in paraview
	!(Optional) Coarse :: logical. If true then only prints the coarse mesh
	!(Optional) Sys (in) :: integer. If 1 then 2 x 2 system of scalar
	subroutine Tri2VTK_Sys_Mesh(meshL, var, relief, filename, T, title, coarse, Sys)
		Implicit none
		!Global variables
		type(mesh), intent(inout), dimension(:) :: meshL
		logical, optional, intent(in) :: relief, coarse
		integer, optional, intent(in) :: var, T, Sys
		character (len = *), optional, intent(in) :: filename, title
		!Local variables
		integer :: sis, var2
		logical :: relieve, corase
		character (len = 200) :: filename2, title2
		!Prepare data
		sis = size(meshL(1)%Tri(1)%UpUF,3)/2

		if (present(Sys)) Sis = Sys
		!*****************Optional area*******************************
		filename2 = "Mesh_Sys.vtk"
		if (present(filename)) filename2 = filename//".vtk"
		title2 = "Semi-structured triangle"
		if (present(title)) title2 = title
		relieve = .false.
		if (present(relief)) relieve = relief
		var2 = 1
		if (present(var)) var2 = var
		corase = .false.
		if (present(coarse)) corase = coarse
		!Call methods
		select case (sis)
		    case (1)
		        if (present(t)) then
		        	call Tri2VTK_MeshCC(meshL,var=var2,relief=relieve,filename=trim(filename2),title=title2,T=t,Sys=Sis)
		        else
		        	call Tri2VTK_MeshCC(meshL,var=var2,relief=relieve,filename=trim(filename2),title=title2,Sys=Sis)
		        end if
		    case (2)
		        if (present(t)) then
		            call Tri2VTK_meshMix(meshL, var2, relieve, trim(filename2), T, title2,coarse=corase)
		        else
		        	call Tri2VTK_meshMix(meshL,var=var2,relief=relieve,filename=trim(filename2),title=title2,coarse=corase)
		        end if
		end select

	end subroutine


	!This subroutine creates an VTK archive to be read by ParaView
	!Only works for Mix nodes
	!Tri(in) :: type(triangle)
	!(Optional) var(in) :: Integer. Select the variable to be plotted
	!								1 -> U(default)
	!								2 -> F
	!								3 -> Residual
	!								4 -> Auxiliar
	!(Optional) relief(in) :: Logical. Choose whether plot in 3d(default) or 2d
	!(Optional) filename(in) :: Character(*). Name of the archive to save into
	!(Optional) T(in) :: Integer. Number to add to the filename. Useful to time dependant problems
	!					or to see diferent iterations. As paraview permits to move along time.
	!(Optional) title(in) :: Character(*). Title to show in paraview
	subroutine Tri2VTK_Mix(Tri, var, relief, filename, T, title)
		Implicit none
		!Global variables
		type(Triangle), intent(in) :: Tri
		logical, optional, intent(in) :: relief
		integer, optional, intent(in) :: var, T
		character (len = *), optional, intent(in) :: filename, title
		!Local variables
		integer :: i,j,k, E_IO, Nn, Nf, pos, var2, a
		logical :: relieve
		double precision, dimension(2) :: v1,v2, aux, point
		real, allocatable, dimension(:,:) :: XYZ
		integer, allocatable, dimension(:) :: connect, type
		double precision, allocatable, dimension(:) :: values
		double precision, allocatable, dimension(:,:,:) :: Aux2, Aux2Down
		character (len = 200) :: filename2, title2, cadena, path
		logical :: bol
		!*****************Optional area*******************************
		filename2 = "Triangle_Mix.vtk"
		if (present(filename)) filename2 = filename//".vtk"
		title2 = "Semi-structured triangle"
		if (present(title)) title2 = title
		relieve = .true.
		if (present(relief)) relieve = relief
		var2 = 1
		if (present(var)) var2 = var
		a = size(Tri%UpUF,3)
		!*************************************************************
		!Number of nodes
		Nn = TotalNodes(Tri,.false.)
		!Number of figures = Up + down triangles
		Nf = Tri%SizeUp**2
		!Allocate arrays
		allocate(XYZ(3,Nn))
		allocate(type(Nf))
		allocate(connect(Nf*7))
		allocate(values(Nn))
		allocate(Aux2(Tri%SizeUp+1,Tri%SizeUp+1,4))
		allocate(Aux2Down(0:Tri%SizeUp+1,0:Tri%SizeUp+1,a))

		!**Select data to plot**
		select case (var2)
		    case (2)!F
		    	do k = 1, 4
				    do i = 2, Tri%SizeUp+1
				        do j = 1, i-1
				            Aux2(j,i,k) = Tri%UpUF(i,j,k)
				        end do
				    end do
			   	end do
			    do i =2, Tri%SizeDown+1
			    	do j = 1, i-1
			    	    Aux2Down(j,i,a) = Tri%DownUF(i,j,a)
			    	end do
			    end do
		    case (3)!Residual
				do k = 1, 4
				    do i = 2, Tri%SizeUp+1
				        do j = 1, i-1
				            Aux2(j,i,k) = Tri%UpResAux(j,i,k)
				        end do
				    end do
			    end do
			    do i =2, Tri%SizeDown+1
			    	do j = 1, i-1
			    	    Aux2Down(j,i,a) = Tri%DownResAux(j,i,a)
			    	end do
			    end do
		    case (4)!Aux
		    	do k = 1, 4
				    do i = 2, Tri%SizeUp+1
				        do j = 1, i-1
				            Aux2(j,i,k) = Tri%UpResAux(i,j,k)
				        end do
				    end do
				end do
			    do i =2, Tri%SizeDown+1
			    	do j = 1, i-1
			    	    Aux2Down(j,i,a) = Tri%DownResAux(i,j,a)
			    	end do
			    end do
		    case default!U
		        !Aux2Up = Tri%UpUF
		        !Aux2Down = Tri%DownUF
				do k = 1, 4
				    do i = 2, Tri%SizeUp+1
				        do j = 1, i-1
				            Aux2(j,i,k) = Tri%UpUF(j,i,k)
				        end do
				    end do
			    end do
			    do i =2, Tri%SizeDown+1
			    	do j = 1, i-1
			    	    Aux2Down(j,i,a) = Tri%DownUF(j,i,a)
			    	end do
			    end do
		end select


		!**Store type**
		!As we are with side nodes only the type is 22, cuadratic triangle
		Type = 22
		!Create file		!output_format = 'ASCII'
		if (present(T)) then
			i = len_trim(filename2) - 3
			write(cadena,*) T
			call removebksl(cadena)
			call insertstr(filename2,cadena,i)
!****************************ONLY WORKS WITH LINUX AND INTEL FORTRAN************************!
!			inquire( DIRECTORY='VTK_Images', exist=bol)
!			if (.not. bol) call system ('mkdir VTK_Images')
!****************************DISABLE IF NECESSARY, AS IT IS NOT CRITICAL*********************
			filename2 = "VTK_Images/"//filename2

		end if
		!Get absolute path
		call getCWD (path)
		!To avoid using debug as working directory
		call delsubstr(path,"/Debug")
		path = trim(path)//"/"//trim(filename2)
        E_IO = VTK_INI(output_format = 'ASCII',filename = path,&
        title = title2, mesh_topology = 'UNSTRUCTURED_GRID')

		!**Position of the nodes**
		!Initial quantity of nodes in the L2 side
		i = 1 + Tri%SizeUP * 2
		!Store coordinates
		v1 = (Tri%Xp2-Tri%Xp1) / (2d0*Tri%SizeUp)
		v2 = (Tri%Xp3-Tri%Xp1) / (2d0*Tri%SizeUp)
		point = Tri%Xp1
		k = 1
		do i = i, 1, -1
			aux = point
			do j = 1, i
				XYZ(1,k) = point(1)
				XYZ(2,k) = point(2)
				!2D
				XYZ(3,k) = 0d0
				k = k + 1
				point = point + v1
			end do
			point = aux + v2
		end do

		!**Connectivity**
		!For up triangles
		pos = 1
		j = 3
		k = Tri%SizeUp * 2 + 1
		do i = 1, getUpTriangles(Tri)*7, 7 !7 points is a triangle. Therefore, one triangle per iteration
		    connect(i) = 6
			connect(i+1) = pos - 1
			connect(i+2) = pos + 1
			connect(i+3) = pos + 2 * ( k - 1)
			connect(i+4) = pos
			connect(i+5) = pos + k
			connect(i+6) = pos + k - 1

			if (j >= k) then
				k = k - 2
				j = 3
				pos = pos + k + 4
			else
				j = j + 2
				pos = pos + 2
			end if
		end do
		!For Down triangles
		pos = 3
		j = 3
		k = Tri%SizeUp * 2 + 1
		do i = getUpTriangles(Tri)*7 + 1, size(connect,1), 7
			connect(i) = 6
			connect(i+1) = pos - 1
			connect(i+2) = pos + 2 * (k - 1)
			connect(i+3) = pos + 2 * (k - 2)
			connect(i+4) = pos - 1 + k
			connect(i+5) = pos + 2 * (k - 1) - 1
			connect(i+6) = pos - 2 + k

			if (j >= k-2) then
				k = k - 2
				j = 3
				pos = pos + k + 6
			else
				j = j + 2
				pos = pos + 2
			end if
		end do
		!**********************The places with !Aux2(col,row,node)...*******************
		!**** are the boundaries, that now are stablished manually as zero!*************
		!**Store data**
		pos = 1
		k = Tri%SizeUp * 2 + 1
		do j = 1, Tri%SizeUp
		    do i = j+1, Tri%SizeUp + 1
		        if (i == j + 1) then
		            values(pos) = 0d0!Aux2(j,i,2)
		            values(pos+1) = Aux2(j,i,2)
		            values(pos+k) = Aux2(j,i,1)
		            values(pos+k+1) = Aux2(j,i,3)
		        else if (i == Tri%SizeUp + 1) then
		            values(pos) = (Aux2(j,i,2) + Aux2(j,i-1,2)) / 2d0
		            values(pos+1) = Aux2(j,i,2)
		            values(pos+k) = Aux2(j,i,1)
		            values(pos+k+1) = Aux2(j,i,3)
		            values(pos+2) = 0d0!Aux2(j,i,2)
		        else
		            values(pos) = (Aux2(j,i,2) + Aux2(j,i-1,2)) / 2d0
		            values(pos+1) = Aux2(j,i,2)
		            values(pos+k) = Aux2(j,i,1)
		            values(pos+k+1) = Aux2(j,i,3)

		        end if
			pos = pos + 2
		    end do
		pos = pos + k
		k = k - 2
		end do
		!manually 2 nodes from Xp3
		values(size(values,1)-3) = 0d0!Aux2(Tri%SizeUp,Tri%SizeUp+1,2)
		values(size(values,1)) = 0d0!Aux2(Tri%SizeUp,Tri%SizeUp+1,2)
		!If 3D
		if (relieve) XYZ(3,:) = values
		!**Introduce data**
		!Coordinates
        E_IO = VTK_GEO(NN = Nn,X=XYZ(1,:),Y=XYZ(2,:),Z=XYZ(3,:))
        !Connectivity
        E_IO = VTK_CON(NC = Nf,connect = connect, cell_type = type)
        !Scalar values
        E_IO = VTK_DAT(NC_NN = Nn,var_location = 'node')
        E_IO = VTK_VAR(NC_NN = Nn, varname = 'Vector_values', var = values)

		!Introduce CC values as values in the faces of those triangles
		deallocate(values)
		allocate(values(Tri%SizeUp**2))
		!**Store data**
		!Up data
		pos = 1
		do j = 1, Tri%SizeUp
		    do i = j+1, Tri%SizeUp + 1
	            values(pos) = Aux2(j,i,a)
				pos = pos + 1
		    end do
		end do
		!down data
		do j = 1, Tri%SizeDown
			do i = j+1, Tri%SizeDown+1
			    values(pos) = Aux2Down(j,i,a)
			    pos = pos + 1
			end do
		end do

        E_IO = VTK_DAT(NC_NN = size(values),var_location = 'cell')
        E_IO = VTK_VAR(NC_NN = size(values), varname = 'Scalar_values', var = values)

        !Close file
        E_IO = VTK_END()

	end subroutine Tri2VTK_Mix


	!This subroutine creates an VTK archive to be read by ParaView
	!Only works for side centered nodes
	!Tri(in) :: type(triangle)
	!(Optional) var(in) :: Integer. Select the variable to be plotted
	!								1 -> U(default)
	!								2 -> F
	!								3 -> Residual
	!								4 -> Auxiliar
	!(Optional) relief(in) :: Logical. Choose whether plot in 3d(default) or 2d
	!(Optional) filename(in) :: Character(*). Name of the archive to save into
	!(Optional) T(in) :: Integer. Number to add to the filename. Useful to time dependant problems
	!					or to see diferent iterations. As paraview permits to move along time.
	!(Optional) title(in) :: Character(*). Title to show in paraview
	!(Optional) Coarse :: logical. If true then only prints the coarse mesh
	subroutine Tri2VTK_MeshMix(meshL, var, relief, filename, T, title, coarse)
		Implicit none
		!Global variables
		type(mesh), intent(inout), dimension(:) :: meshL
		logical, optional, intent(in) :: relief, coarse
		integer, optional, intent(in) :: var, T
		character (len = *), optional, intent(in) :: filename, title
		!Local variables
		integer :: i,j,k, E_IO, Nn, Nf, pos, var2, Mpos, a, NOUSAR
		logical :: relieve
		double precision, dimension(2) :: v1,v2, aux, point
		real, allocatable, dimension(:,:) :: XYZ
		integer, allocatable, dimension(:) :: connect, type
		double precision, allocatable, dimension(:) :: values
		double precision, allocatable, dimension(:,:,:) :: Aux2, Aux2Down
		character (len = 200) :: filename2, title2, cadena, path
		logical :: bol, coar
		!*****************Optional area*******************************
		filename2 = "Mesh_Mix.vtk"
		if (present(filename)) filename2 = filename
		if (index(trim(filename2),".vtk")==0) filename2 = filename2//".vtk"
		title2 = "Semi-structured triangle"
		if (present(title)) title2 = title
		relieve = .true.
		if (present(relief)) relieve = relief
		var2 = 1
		if (present(var)) var2 = var
		NOUSAR = var2
		coar = .false.
		if (present(Coarse)) coar = Coarse
		!*************************************************************
		a = size(meshL(1)%Tri(1)%UpUF,3)
		!Number of nodes
		Nn = TotalNodes(meshL(1)%Tri(1),.false.)
		!Number of figures = Up + down triangles
		Nf = meshL(1)%Tri(1)%SizeUp**2
		!Allocate arrays
		allocate(XYZ(3,Nn*size(meshL)))
		allocate(type(Nf*size(meshL)))
		allocate(connect(Nf*7*size(meshL)))
		allocate(values(Nn*size(meshL)))
		allocate(Aux2(meshL(1)%Tri(1)%SizeUp+1,meshL(1)%Tri(1)%SizeUp+1,4))
		allocate(Aux2Down(meshL(1)%Tri(1)%SizeUp+1,meshL(1)%Tri(1)%SizeUp+1,4:4))
		!**Store type**
		!As we are with side nodes only the type is 22, cuadratic triangle
		Type = 22
		!Create file		!output_format = 'ASCII'
		if (present(T)) then
			i = len_trim(filename2) - 3
			write(cadena,*) T
			call removebksl(cadena)
			call insertstr(filename2,cadena,i)
!****************************ONLY WORKS WITH LINUX AND INTEL FORTRAN************************!
!			inquire( DIRECTORY='VTK_Images', exist=bol)
!			if (.not. bol) call system ('mkdir VTK_Images')
!****************************DISABLE IF NECESSARY, AS IT IS NOT CRITICAL*********************
			filename2 = "VTK_Images/"//filename2

		end if
		!Get absolute path
		call getCWD (path)
		!To avoid using debug as working directory
		call delsubstr(path,"/Debug")
		path = trim(path)//"/"//trim(filename2)
        E_IO = VTK_INI(output_format = 'ASCII',filename = path,&
        title = title2, mesh_topology = 'UNSTRUCTURED_GRID')
		do Mpos = 1, size(MeshL)
		!**Select data to plot**
		select case (NOUSAR)
		    case (2)!F
		    	do k = 1, 3
				    do i = 2, meshL(Mpos)%Tri(1)%SizeUp+1
				        do j = 1, i-1
				            Aux2(j,i,k) = meshL(Mpos)%Tri(1)%UpUF(i,j,k)
				        end do
				    end do
			   	end do
		    case (3)!Residual
				do k = 1, 3
				    do i = 2, meshL(Mpos)%Tri(1)%SizeUp+1
				        do j = 1, i-1
				            Aux2(j,i,k) = meshL(Mpos)%Tri(1)%UpResAux(j,i,k)
				        end do
				    end do
			    end do
		    case (4)!Aux
		    	do k = 1, 3
				    do i = 2, meshL(Mpos)%Tri(1)%SizeUp+1
				        do j = 1, i-1
				            Aux2(j,i,k) = meshL(Mpos)%Tri(1)%UpResAux(i,j,k)
				        end do
				    end do
				end do
		    case default!U
!		        Aux2Up = Tri%UpUF
!		        Aux2Down = Tri%DownUF
				do k = 1, 3
				    do i = 2, meshL(Mpos)%Tri(1)%SizeUp+1
				        do j = 1, i-1
				            Aux2(j,i,k) = meshL(Mpos)%Tri(1)%UpUF(j,i,k)
				        end do
				    end do
			    end do
		end select

			!Start node for current triangle
			i = (Mpos - 1) * Nn + 1
			E_IO = (Mpos - 1) * Nf + 1
			!End node for current triangle
			j = Nn * Mpos
			aux = Nf * Mpos
			!Add one triangle to our file
			call DataCoordinates(meshL, Mpos, XYZ(:,i:j), values(i:j), relieve, Aux2,.false. ,var = var2)

			!**Connectivity**
			call Connectivity(meshL,Mpos,connect, .false.)

		end do

		!If 3D
		if (relieve) XYZ(3,:) = values

		!**Only Print coarse grid**
		if (coar) then
	        do Mpos = 1, size(meshL)
	        	!Start node for current triangle
				i = (Mpos - 1) * Nn + 1
				aux = (Mpos - 1) * 7 * meshL(Mpos)%Tri(1)%SizeUp**2
				!End node for current triangle
				j = Nn * Mpos
				XYZ(1,i:j)=meshL(Mpos)%Xp1(1)
				XYZ(2,i:j)=meshL(Mpos)%Xp1(2)
				XYZ(3,i:j)=0d0
				!Except one triangle
				XYZ(1,connect(aux+1)+1)=meshL(Mpos)%Xp1(1)
				XYZ(2,connect(aux+1)+1)=meshL(Mpos)%Xp1(2)

				XYZ(1,connect(aux+2)+1)=meshL(Mpos)%Xp2(1)
				XYZ(2,connect(aux+2)+1)=meshL(Mpos)%Xp2(2)

				XYZ(1,connect(aux+3)+1)=meshL(Mpos)%Xp3(1)
				XYZ(2,connect(aux+3)+1)=meshL(Mpos)%Xp3(2)

				XYZ(1,connect(aux+4)+1)=meshL(Mpos)%Xp3(1)
				XYZ(2,connect(aux+4)+1)=meshL(Mpos)%Xp3(2)

				XYZ(1,connect(aux+5)+1)=meshL(Mpos)%Xp3(1)
				XYZ(2,connect(aux+5)+1)=meshL(Mpos)%Xp3(2)

				XYZ(1,connect(aux+6)+1)=meshL(Mpos)%Xp1(1)
				XYZ(2,connect(aux+6)+1)=meshL(Mpos)%Xp1(2)
	        end do
		end if

		!**Introduce data**
		!Coordinates
        E_IO = VTK_GEO(NN = size(XYZ,2),X=XYZ(1,:),Y=XYZ(2,:),Z=XYZ(3,:))
        !Connectivity
        E_IO = VTK_CON(NC = size(connect)/7,connect = connect, cell_type = type)
        !Scalar values
        E_IO = VTK_DAT(NC_NN = size(XYZ,2),var_location = 'node')
        E_IO = VTK_VAR(NC_NN = size(XYZ,2), varname = 'Vectorial_values', var = values)

		!********************Add more data**********************
		!!!!You cannot store two values consecutively, of the same type: cell, node!!!!
		!***Introduce CC values as values in the faces of those triangles****
		deallocate(values)
		var2 = meshL(1)%Tri(1)%SizeUp**2
		allocate(values(size(meshL)*var2))

        do Mpos = 1, size(meshL)

			!**Select data to plot**
			select case (NOUSAR)
			    case (2)!F
				    do i = 2, meshL(Mpos)%Tri(1)%SizeUp+1
				        do j = 1, i-1
				            Aux2(j,i,a) = meshL(Mpos)%Tri(1)%UpUF(i,j,a)
				        end do
				    end do
				    do i =2, meshL(Mpos)%Tri(1)%SizeDown+1
				    	do j = 1, i-1
				    	    Aux2Down(j,i,a) = meshL(Mpos)%Tri(1)%DownUF(i,j,a)
				    	end do
				    end do
			    case (3)!Residual
				    do i = 2, meshL(Mpos)%Tri(1)%SizeUp+1
				        do j = 1, i-1
				            Aux2(j,i,a) = meshL(Mpos)%Tri(1)%UpResAux(j,i,a)
				        end do
				    end do
				    do i =2, meshL(Mpos)%Tri(1)%SizeDown+1
				    	do j = 1, i-1
				    	    Aux2Down(j,i,a) = meshL(Mpos)%Tri(1)%DownResAux(j,i,a)
				    	end do
				    end do
			    case (4)!Aux
				    do i = 2, meshL(Mpos)%Tri(1)%SizeUp+1
				        do j = 1, i-1
				            Aux2(j,i,a) = meshL(Mpos)%Tri(1)%UpResAux(i,j,a)
				        end do
				    end do
				    do i =2, meshL(Mpos)%Tri(1)%SizeDown+1
				    	do j = 1, i-1
				    	    Aux2Down(j,i,a) = meshL(Mpos)%Tri(1)%DownResAux(i,j,a)
				    	end do
				    end do
			    case default!U
				    do i = 2, meshL(Mpos)%Tri(1)%SizeUp+1
				        do j = 1, i-1
				            Aux2(j,i,a) = meshL(Mpos)%Tri(1)%UpUF(j,i,a)
				        end do
				    end do
				    do i =2, meshL(Mpos)%Tri(1)%SizeDown+1
				    	do j = 1, i-1
				    	    Aux2Down(j,i,a) = meshL(Mpos)%Tri(1)%DownUF(j,i,a)
				    	end do
				    end do
			end select

        	!Start node for current triangle
			aux = (Mpos - 1) * var2 + 1
			!**Store data**
			!Up data
			do j = 1, meshL(Mpos)%Tri(1)%SizeUp
			    do i = j+1, meshL(Mpos)%Tri(1)%SizeUp + 1
		            values(aux) = Aux2(j,i,a)
					aux = aux + 1
			    end do
			end do
			!down data
			do j = 1, meshL(Mpos)%Tri(1)%SizeDown
				do i = j+1, meshL(Mpos)%Tri(1)%SizeDown+1
				    values(aux) = Aux2Down(j,i,a)
				    aux = aux + 1
				end do
			end do

		end do
        E_IO = VTK_DAT(NC_NN = size(values),var_location = 'cell')
        E_IO = VTK_VAR(NC_NN = size(values), varname = 'Scalar_values', var = values)

		if (size(meshL(1)%Tri)<=5) then
			!K coefficient used in each triangle
	        do Mpos = 1, size(meshL)
	        	!Start node for current triangle
				i = (Mpos - 1) * Nn + 1
				!aux = (Mpos - 1) * 7 * meshL(Mpos)%Tri(1)%SizeUp**2
				!End node for current triangle
				j = Nn * Mpos
				XYZ(2,i:j)=meshL(Mpos)%K_coef
	        end do
		        E_IO = VTK_DAT(NC_NN = size(XYZ,2),var_location = 'node')!The other option is cell
		        E_IO = VTK_VAR(NC_NN = size(XYZ,2), varname = 'K_Coefficient', var = int(XYZ(2,:)))
			!Smoothers used in each triangle
			deallocate(values)
			var2 = meshL(1)%Tri(1)%SizeUp**2
			allocate(values(size(meshL)*var2))

	        do pos = 1, size(meshL)
	        	!Start node for current triangle
				i = (pos - 1) * var2 + 1
				!End node for current triangle
				j = var2 * pos
		        values(i:j) = int(meshL(pos)%method/10)*10
	        end do
	        E_IO = VTK_DAT(NC_NN = size(values),var_location = 'cell')
	        E_IO = VTK_VAR(NC_NN = size(values), varname = 'Smoother', var = int(values))
		end if

        !Close file
        E_IO = VTK_END()
	end subroutine Tri2VTK_MeshMix

	!This subroutine creates an VTK archive to be read by ParaView
	!Only works for side centered nodes
	!Tri(in) :: type(triangle)
	!(Optional) var(in) :: Integer. Select the variable to be plotted
	!								1 -> U(default)
	!								2 -> F
	!								3 -> Residual
	!								4 -> Auxiliar
	!(Optional) relief(in) :: Logical. Choose whether plot in 3d(default) or 2d
	!(Optional) filename(in) :: Character(*). Name of the archive to save into
	!(Optional) T(in) :: Integer. Number to add to the filename. Useful to time dependant problems
	!					or to see diferent iterations. As paraview permits to move along time.
	!(Optional) title(in) :: Character(*). Title to show in paraview
	subroutine Tri2VTK_SD(Tri, var, relief, filename, T, title)
		Implicit none
		!Global variables
		type(Triangle), intent(in) :: Tri
		logical, optional, intent(in) :: relief
		integer, optional, intent(in) :: var, T
		character (len = *), optional, intent(in) :: filename, title
		!Local variables
		integer :: i,j,k, E_IO, Nn, Nf, pos, var2, a
		logical :: relieve
		double precision, dimension(2) :: v1,v2, aux, point
		real, allocatable, dimension(:,:) :: XYZ
		integer, allocatable, dimension(:) :: connect, type
		double precision, allocatable, dimension(:) :: values
		double precision, allocatable, dimension(:,:,:) :: Aux2
		character (len = 200) :: filename2, title2, cadena, path
		logical :: bol
		!*****************Optional area*******************************
		filename2 = "Triangle_SD.vtk"
		if (present(filename)) filename2 = filename//".vtk"
		title2 = "Semi-structured triangle"
		if (present(title)) title2 = title
		relieve = .true.
		if (present(relief)) relieve = relief
		var2 = 1
		if (present(var)) var2 = var
		a = size(Tri%UpUF,3)
		!*************************************************************
		!Number of nodes
		Nn = TotalNodes(Tri,.false.)
		!Number of figures = Up + down triangles
		Nf = Tri%SizeUp**2
		!Allocate arrays
		allocate(XYZ(3,Nn))
		allocate(type(Nf))
		allocate(connect(Nf*7))
		allocate(values(Nn))
		allocate(Aux2(Tri%SizeUp+1,Tri%SizeUp+1,3))
		!**Select data to plot**
		select case (var2)
		    case (2)!F
				do k = 1, 3
				    do i = 2, Tri%SizeUp+1
				        do j = 1, i-1
				            Aux2(j,i,k) = Tri%UpUF(i,j,k)
				        end do
				    end do
				end do
		    case (3)!Residual
		    	Aux2 = Tri%UpResAux
		    case (4)!Aux
				do k = 1, 3
				    do i = 2, Tri%SizeUp+1
				        do j = 1, i-1
				            Aux2(j,i,k) = Tri%UpResAux(i,j,k)
				        end do
				    end do
				end do
		    case default!U
		        Aux2 = Tri%UpUF
		end select


		!**Store type**
		!As we are with side nodes only the type is 22, cuadratic triangle
		Type = 22
		!Create file		!output_format = 'ASCII'
		if (present(T)) then
			i = len_trim(filename2) - 3
			write(cadena,*) T
			call removebksl(cadena)
			call insertstr(filename2,cadena,i)
!****************************ONLY WORKS WITH LINUX AND INTEL FORTRAN************************!
!			inquire( DIRECTORY='VTK_Images', exist=bol)
!			if (.not. bol) call system ('mkdir VTK_Images')
!****************************DISABLE IF NECESSARY, AS IT IS NOT CRITICAL*********************
			filename2 = "VTK_Images/"//filename2

		end if
		!Get absolute path
		call getCWD (path)
		!To avoid using debug as working directory
		call delsubstr(path,"/Debug")
		path = trim(path)//"/"//trim(filename2)
        E_IO = VTK_INI(output_format = 'ASCII',filename = path,&
        title = title2, mesh_topology = 'UNSTRUCTURED_GRID')

		!**Position of the nodes**
		!Initial quantity of nodes in the L2 side
		i = 1 + Tri%SizeUP * 2
		!Store coordinates
		v1 = (Tri%Xp2-Tri%Xp1) / (2d0*Tri%SizeUp)
		v2 = (Tri%Xp3-Tri%Xp1) / (2d0*Tri%SizeUp)
		point = Tri%Xp1
		k = 1
		do i = i, 1, -1
			aux = point
			do j = 1, i
				XYZ(1,k) = point(1)
				XYZ(2,k) = point(2)
				!2D
				XYZ(3,k) = 0d0
				k = k + 1
				point = point + v1
			end do
			point = aux + v2
		end do

		!**Connectivity**
		!For up triangles
		pos = 1
		j = 3
		k = Tri%SizeUp * 2 + 1
		do i = 1, getUpTriangles(Tri)*7, 7 !7 points is a triangle. Therefore, one triangle per iteration
		    connect(i) = 6
			connect(i+1) = pos - 1
			connect(i+2) = pos + 1
			connect(i+3) = pos + 2 * ( k - 1)
			connect(i+4) = pos
			connect(i+5) = pos + k
			connect(i+6) = pos + k - 1

			if (j >= k) then
				k = k - 2
				j = 3
				pos = pos + k + 4
			else
				j = j + 2
				pos = pos + 2
			end if
		end do
		!For Down triangles
		pos = 3
		j = 3
		k = Tri%SizeUp * 2 + 1
		do i = getUpTriangles(Tri)*7 + 1, size(connect,1), 7
			connect(i) = 6
			connect(i+1) = pos - 1
			connect(i+2) = pos + 2 * (k - 1)
			connect(i+3) = pos + 2 * (k - 2)
			connect(i+4) = pos - 1 + k
			connect(i+5) = pos + 2 * (k - 1) - 1
			connect(i+6) = pos - 2 + k

			if (j >= k-2) then
				k = k - 2
				j = 3
				pos = pos + k + 6
			else
				j = j + 2
				pos = pos + 2
			end if
		end do
		!**********************The places with !Aux2(col,row,node)...*******************
		!**** are the boundaries, that now are stablished manually as zero!*************
		!**Store data**
		pos = 1
		k = Tri%SizeUp * 2 + 1
		do j = 1, Tri%SizeUp
		    do i = j+1, Tri%SizeUp + 1
		        if (i == j + 1) then
		            values(pos) = 0d0!Aux2(j,i,2)
		            values(pos+1) = Aux2(j,i,2)
		            values(pos+k) = Aux2(j,i,1)
		            values(pos+k+1) = Aux2(j,i,3)
		        else if (i == Tri%SizeUp + 1) then
		            values(pos) = (Aux2(j,i,2) + Aux2(j,i-1,2)) / 2d0
		            values(pos+1) = Aux2(j,i,2)
		            values(pos+k) = Aux2(j,i,1)
		            values(pos+k+1) = Aux2(j,i,3)
		            values(pos+2) = 0d0!Aux2(j,i,2)
		        else
		            values(pos) = (Aux2(j,i,2) + Aux2(j,i-1,2)) / 2d0
		            values(pos+1) = Aux2(j,i,2)
		            values(pos+k) = Aux2(j,i,1)
		            values(pos+k+1) = Aux2(j,i,3)

		        end if
			pos = pos + 2
		    end do
		pos = pos + k
		k = k - 2
		end do
		!manually 2 nodes from Xp3
		values(size(values,1)-3) = 0d0!Aux2(Tri%SizeUp,Tri%SizeUp+1,2)
		values(size(values,1)) = 0d0!Aux2(Tri%SizeUp,Tri%SizeUp+1,2)
		!If 3D
		if (relieve) XYZ(3,:) = values
		!**Introduce data**
		!Coordinates
        E_IO = VTK_GEO(NN = Nn,X=XYZ(1,:),Y=XYZ(2,:),Z=XYZ(3,:))
        !Connectivity
        E_IO = VTK_CON(NC = Nf,connect = connect, cell_type = type)
        !Scalar values
        E_IO = VTK_DAT(NC_NN = Nn,var_location = 'node')
        E_IO = VTK_VAR(NC_NN = Nn, varname = 'U_values', var = values)

        !Close file
        E_IO = VTK_END()
	end subroutine Tri2VTK_SD

	!This subroutine creates an VTK archive to be read by ParaView
	!Only works for side centered nodes
	!Tri(in) :: type(triangle)
	!(Optional) var(in) :: Integer. Select the variable to be plotted
	!								1 -> U(default)
	!								2 -> F
	!								3 -> Residual
	!								4 -> Auxiliar
	!(Optional) relief(in) :: Logical. Choose whether plot in 3d(default) or 2d
	!(Optional) filename(in) :: Character(*). Name of the archive to save into
	!(Optional) T(in) :: Integer. Number to add to the filename. Useful to time dependant problems
	!					or to see diferent iterations. As paraview permits to move along time.
	!(Optional) title(in) :: Character(*). Title to show in paraview
	!(Optional) Coarse :: logical. If true then only prints the coarse mesh
	subroutine Tri2VTK_MeshSD(meshL, var, relief, filename, T, title, coarse)
		Implicit none
		!Global variables
		type(mesh), intent(inout), dimension(:) :: meshL
		logical, optional, intent(in) :: relief, coarse
		integer, optional, intent(in) :: var, T
		character (len = *), optional, intent(in) :: filename, title
		!Local variables
		integer :: i,j,k, E_IO, Nn, Nf, pos, var2, Mpos
		logical :: relieve
		double precision, dimension(2) :: v1,v2, aux, point
		real, allocatable, dimension(:,:) :: XYZ
		integer, allocatable, dimension(:) :: connect, type
		double precision, allocatable, dimension(:) :: values
		double precision, allocatable, dimension(:,:,:) :: Aux2
		character (len = 200) :: filename2, title2, cadena, path
		logical :: bol, coar
		!*****************Optional area*******************************
		filename2 = "Mesh_SD.vtk"
		if (present(filename)) filename2 = filename//".vtk"
		title2 = "Semi-structured triangle"
		if (present(title)) title2 = title
		relieve = .true.
		if (present(relief)) relieve = relief
		var2 = 1
		if (present(var)) var2 = var
		coar = .false.
		if (present(Coarse)) coar = Coarse
		!*************************************************************
		!Number of nodes
		Nn = TotalNodes(meshL(1)%Tri(1),.false.)
		!Number of figures = Up + down triangles
		Nf = meshL(1)%Tri(1)%SizeUp**2
		!Allocate arrays
		allocate(XYZ(3,Nn*size(meshL)))
		allocate(type(Nf*size(meshL)))
		allocate(connect(Nf*7*size(meshL)))
		allocate(values(Nn*size(meshL)))
		allocate(Aux2(meshL(1)%Tri(1)%SizeUp+1,meshL(1)%Tri(1)%SizeUp+1,3))
		!**Store type**
		!As we are with side nodes only the type is 22, cuadratic triangle
		Type = 22
		!Create file		!output_format = 'ASCII'
		if (present(T)) then
			i = len_trim(filename2) - 3
			write(cadena,*) T
			call removebksl(cadena)
			call insertstr(filename2,cadena,i)
!****************************ONLY WORKS WITH LINUX AND INTEL FORTRAN************************!
!			inquire( DIRECTORY='VTK_Images', exist=bol)
!			if (.not. bol) call system ('mkdir VTK_Images')
!****************************DISABLE IF NECESSARY, AS IT IS NOT CRITICAL*********************
			filename2 = "VTK_Images/"//filename2

		end if
		!Get absolute path
		call getCWD (path)
		!To avoid using debug as working directory
		call delsubstr(path,"/Debug")
		path = trim(path)//"/"//trim(filename2)
        E_IO = VTK_INI(output_format = 'ASCII',filename = path,&
        title = title2, mesh_topology = 'UNSTRUCTURED_GRID')
		do Mpos = 1, size(MeshL)
			!**Select data to plot**
			select case (var2)
			    case (2)!F
					do k = 1, 3
					    do i = 2, meshL(Mpos)%Tri(1)%SizeUp+1
					        do j = 1, i-1
					            Aux2(j,i,k) = meshL(Mpos)%Tri(1)%UpUF(i,j,k)
					        end do
					    end do
					end do
			    case (3)!Residual
					do k = 1, 3
					    do i = 2, meshL(Mpos)%Tri(1)%SizeUp+1
					        do j = 1, i-1
					            Aux2(j,i,k) = meshL(Mpos)%Tri(1)%UpResAux(j,i,k)
					        end do
					    end do
					end do
			    case (4)!Aux
					do k = 1, 3
					    do i = 2, meshL(Mpos)%Tri(1)%SizeUp+1
					        do j = 1, i-1
					            Aux2(j,i,k) = meshL(Mpos)%Tri(1)%UpResAux(i,j,k)
					        end do
					    end do
					end do
			    case default!U
			    	do k = 1, 3
					    do i = 2, meshL(Mpos)%Tri(1)%SizeUp+1
					        do j = 1, i-1
					            Aux2(j,i,k) = meshL(Mpos)%Tri(1)%UpUF(j,i,k)
					        end do
					    end do
					end do
!			        Aux2 = meshL(Mpos)%Tri(1)%UpUF
			end select

			!Start node for current triangle
			i = (Mpos - 1) * Nn + 1
			E_IO = (Mpos - 1) * Nf + 1
			!End node for current triangle
			j = Nn * Mpos
			aux = Nf * Mpos
			!Add one triangle to our file
			call DataCoordinates(meshL, Mpos, XYZ(:,i:j), values(i:j), relieve, Aux2,.false. ,var = var2)

			!**Connectivity**
			call Connectivity(meshL,Mpos,connect, .false.)

		end do

		!If 3D
		if (relieve) XYZ(3,:) = values

		!**Only Print coarse grid**
		if (coar) then
	        do Mpos = 1, size(meshL)
	        	!Start node for current triangle
				i = (Mpos - 1) * Nn + 1
				aux = (Mpos - 1) * 7 * meshL(Mpos)%Tri(1)%SizeUp**2
				!End node for current triangle
				j = Nn * Mpos
				XYZ(1,i:j)=meshL(Mpos)%Xp1(1)
				XYZ(2,i:j)=meshL(Mpos)%Xp1(2)
				XYZ(3,i:j)=0d0
				!Except one triangle
				XYZ(1,connect(aux+1)+1)=meshL(Mpos)%Xp1(1)
				XYZ(2,connect(aux+1)+1)=meshL(Mpos)%Xp1(2)

				XYZ(1,connect(aux+2)+1)=meshL(Mpos)%Xp2(1)
				XYZ(2,connect(aux+2)+1)=meshL(Mpos)%Xp2(2)

				XYZ(1,connect(aux+3)+1)=meshL(Mpos)%Xp3(1)
				XYZ(2,connect(aux+3)+1)=meshL(Mpos)%Xp3(2)

				XYZ(1,connect(aux+4)+1)=meshL(Mpos)%Xp3(1)
				XYZ(2,connect(aux+4)+1)=meshL(Mpos)%Xp3(2)

				XYZ(1,connect(aux+5)+1)=meshL(Mpos)%Xp3(1)
				XYZ(2,connect(aux+5)+1)=meshL(Mpos)%Xp3(2)

				XYZ(1,connect(aux+6)+1)=meshL(Mpos)%Xp1(1)
				XYZ(2,connect(aux+6)+1)=meshL(Mpos)%Xp1(2)
	        end do
		end if

		!**Introduce data**
		!Coordinates
        E_IO = VTK_GEO(NN = size(XYZ,2),X=XYZ(1,:),Y=XYZ(2,:),Z=XYZ(3,:))
        !Connectivity
        E_IO = VTK_CON(NC = size(connect)/7,connect = connect, cell_type = type)
        !Scalar values
        E_IO = VTK_DAT(NC_NN = size(XYZ,2),var_location = 'node')
        E_IO = VTK_VAR(NC_NN = size(XYZ,2), varname = 'Vectorial_values', var = values)

		!********************Add more data**********************
		!!!!You cannot store two values consecutively, of the same type: cell, node!!!!
		deallocate(values)
		var2 = meshL(1)%Tri(1)%SizeUp**2
		allocate(values(size(meshL)*var2))
		!Smoothers used in each triangle
        do pos = 1, size(meshL)
        	!Start node for current triangle
			i = (pos - 1) * var2 + 1
			!End node for current triangle
			j = var2 * pos
	        values(i:j) = int(meshL(pos)%method/10)*10
        end do
        E_IO = VTK_DAT(NC_NN = size(values),var_location = 'cell')
        E_IO = VTK_VAR(NC_NN = size(values), varname = 'Smoother', var = int(values))


        !Close file
        E_IO = VTK_END()
	end subroutine Tri2VTK_MeshSD



	!This subroutine creates an VTK archive to be read by ParaView
	!Only works for cell centered nodes
	!Tri(in) :: type(triangle)
	!(Optional) var(in) :: Integer. Select the variable to be plotted
	!								1 -> U(default)
	!								2 -> F
	!								3 -> Residual
	!								4 -> Auxiliar
	!(Optional) relief(in) :: Logical. Choose whether plot in 3d(default) or 2d
	!(Optional) filename(in) :: Character(*). Name of the archive to save into
	!(Optional) T(in) :: Integer. Number to add to the filename. Useful to time dependant problems
	!					or to see diferent iterations. As paraview permits to move along time.
	!(Optional) title(in) :: Character(*). Title to show in paraview
	!(Optional) CellColor :: logical. If true then the triangles will have a value for all the surface
	!								  if false(default) the data will be assigned to the nodes.
	!(Optional) Sys (in) :: integer. If 1 then 2 x 2 system of scalar
	subroutine Tri2VTK_CC(Tri, var, relief, filename,T , title, CellColor, Sys)
		Implicit none
		!Global variables
		type(Triangle), intent(in) :: Tri
		logical, optional, intent(in) :: relief, CellColor
		integer, optional, intent(in) :: var, T, Sys
		character (len = *), optional, intent(in) :: filename, title
		!Local variables
		integer :: i,j,k, E_IO, Nn, Nf, pos, var2, a, Sis
		logical :: relieve, Celula
		double precision, dimension(2) :: v1,v2, aux, point, Xvor
		real, allocatable, dimension(:,:) :: XYZ, Z
		integer, allocatable, dimension(:) :: connect, type
		double precision, allocatable, dimension(:) :: values
		double precision, allocatable, dimension(:,:,:) :: Aux2Up, Aux2Down
		character (len = 250) :: filename2, title2, cadena
		character (len=500) :: path
		logical :: bol
		!*****************Optional area*******************************
		filename2 = "Triangle_CC.vtk"
		if (present(filename)) filename2 = filename//".vtk"
		title2 = "Structured triangle"
		if (present(title)) title2 = title
		relieve = .true.
		if (present(relief)) relieve = relief
		var2 = 1
		if (present(var)) var2 = var
		Celula = .false.
		if (present(CellColor)) Celula = CellColor
		Sis = 0
		if (present(Sys)) Sis = Sys
		!*************************************************************
		a = size(Tri%UpUF,3)
		if (Sis==1) a = 1
		!Number of nodes
		Nn = TotalNodes(Tri,.true.) + Tri%SizeUp**2
		!Number of figures = Up + down triangles
		Nf = Tri%SizeUp**2
		!Allocate arrays
		allocate(XYZ(3,Nn))
		allocate(Z(1,Nn))
		allocate(type(Nf))
		allocate(connect(Nf*5))
		allocate(values(Nf))
		allocate(Aux2Up(Tri%SizeUp+1,Tri%SizeUp+1,1))
		allocate(Aux2Down(0:Tri%SizeUp+1,0:Tri%SizeUp+1,1))

		!**Select data to plot**
		select case (var2)
		    case (2)!F
			    do i = 2, Tri%SizeUp+1
			        do j = 1, i-1
			            Aux2Up(j,i,1) = Tri%UpUF(i,j,a)
			        end do
			    end do
			    do i =2, Tri%SizeDown+1
			    	do j = 1, i-1
			    	    Aux2Down(j,i,1) = Tri%DownUF(i,j,a)
			    	end do
			    end do
		    case (3)!Residual
!		    	Aux2Up = Tri%UpResAux
!		    	Aux2Down = Tri%DownResAux
			    do i = 2, Tri%SizeUp+1
			        do j = 1, i-1
			            Aux2Up(j,i,1) = Tri%UpResAux(j,i,a)
			        end do
			    end do
			    do i =2, Tri%SizeDown+1
			    	do j = 1, i-1
			    	    Aux2Down(j,i,1) = Tri%DownResAux(j,i,a)
			    	end do
			    end do
		    case (4)!Aux
			    do i = 2, Tri%SizeUp+1
			        do j = 1, i-1
			            Aux2Up(j,i,1) = Tri%UpResAux(i,j,a)
			        end do
			    end do
			    do i =2, Tri%SizeDown+1
			    	do j = 1, i-1
			    	    Aux2Down(j,i,1) = Tri%DownResAux(i,j,a)
			    	end do
			    end do
		    case default!U
			    do i = 2, Tri%SizeUp+1
			        do j = 1, i-1
			            Aux2Up(j,i,1) = Tri%UpUF(j,i,a)
			        end do
			    end do
			    do i =2, Tri%SizeDown+1
			    	do j = 1, i-1
			    	    Aux2Down(j,i,1) = Tri%DownUF(j,i,a)
			    	end do
			    end do
		end select

		!**Store type**
		!As we are with center nodes only, the type is 10, tetrahedron
		Type = 10
		!Create file		!output_format = 'BINARY'
		if (present(T)) then
			i = len_trim(filename2) - 3
			write(cadena,*) T
			call removebksl(cadena)
			call insertstr(filename2,cadena,i)
!****************************ONLY WORKS WITH LINUX AND INTEL FORTRAN************************!
!			inquire( DIRECTORY='VTK_Images', exist=bol)
!			if (.not. bol) call system ('mkdir VTK_Images')
!****************************DISABLE IF NECESSARY, AS IT IS NOT CRITICAL*********************
			filename2 = "VTK_Images/"//filename2

		end if
		!Get absolute path
		call getCWD (path)
		!To avoid using debug as working directory
		call delsubstr(path,"/Debug")
		path = trim(path)//"/"//filename2
        E_IO = VTK_INI(output_format = 'ASCII',filename = filename2,&
        title = title2, mesh_topology = 'UNSTRUCTURED_GRID')
		!**Position of the nodes**
		i = 1 + Tri%SizeUP
		!Store coordinates
		v1 = (Tri%Xp2-Tri%Xp1) / (Tri%SizeUp)
		v2 = (Tri%Xp3-Tri%Xp1) / (Tri%SizeUp)
		point = Tri%Xp1
		k = 1
		do i = i, 1, -1
			aux = point
			do j = 1, i
				XYZ(1,k) = point(1)
				XYZ(2,k) = point(2)
				!2D
				XYZ(3,k) = 0d0
				k = k + 1
				point = point + v1
			end do
			point = aux + v2
		end do
		!Position of the voronoi center for up triangles
		point = Tri%Xp1
		do j = 1, Tri%SizeUP
			aux = point
			do i = j+1, Tri%SizeUP+1
				call getVoronoiCenter(point,point+v1,point+v2,Xvor)
				XYZ(1,k) = Xvor(1)
				XYZ(2,k) = Xvor(2)

				!2D
				XYZ(3,k) = 0d0
				k = k + 1
				point = point + v1
			end do
			point = aux + v2
		end do
		!Position of the voronoi center for down triangles
		point = Tri%Xp1 + v1
		do j = 1, Tri%SizeDown
			aux = point
			do i = j+1, Tri%SizeDown+1
				call getVoronoiCenter(point,point+v2,point+v2-v1,Xvor)
				XYZ(1,k) = Xvor(1)
				XYZ(2,k) = Xvor(2)
				!2D
				XYZ(3,k) = 0d0
				k = k + 1
				point = point + v1
			end do
			point = aux + v2
		end do
		!**Connectivity**
		!For up triangles
		Var2 = TotalNodes(Tri,.true.)
		pos = 1
		j = 2
		k = Tri%SizeUp + 1
		do i = 1, getUpTriangles(Tri)*5, 5 !5 points is a triangle. Therefore, one triangle per iteration
		    connect(i) = 4
			connect(i+1) = pos - 1
			connect(i+2) = pos
			connect(i+3) = pos - 1 + k
			connect(i+4) = Var2
			Var2 = Var2 + 1
			if (j >= k) then
				k = k - 1
				j = 2
				pos = pos + 2
			else
				j = j + 1
				pos = pos + 1
			end if
		end do
		!For Down triangles
		pos = 2
		j = 2
		k = Tri%SizeUp
		do i = getUpTriangles(Tri)*5 + 1, size(connect,1), 5
			connect(i) = 4
			connect(i+1) = pos - 1
			connect(i+2) = pos + k
			connect(i+3) = pos + (k - 1)
			connect(i+4) = Var2
			Var2 = Var2 + 1
			if (j >= k) then
				k = k - 1
				j = 2
				pos = pos + 3
			else
				j = j + 1
				pos = pos + 1
			end if
		end do
		!**Store data**
		!Up data
		pos = 1
		do j = 1, Tri%SizeUp
		    do i = j+1, Tri%SizeUp + 1
	            values(pos) = Aux2Up(j,i,1)
				pos = pos + 1
		    end do
		end do
		!down data
		do j = 1, Tri%SizeDown
			do i = j+1, Tri%SizeDown+1
			    values(pos) = Aux2Down(j,i,1)
			    pos = pos + 1
			end do
		end do


		k = 1
		var2 = TotalNodes(Tri,.true.) + 1
		do j = 1, Tri%SizeUp
		    do i = j+1, Tri%SizeUp + 1
		    	!Vertex node
	    		XYZ(3,k) = (Aux2Up(j,i,1)+Aux2Down(j,i-1,1))/2d0
				k = k + 1
				!Voronoi node
				XYZ(3,var2) = Aux2Up(j,i,1)
				var2 = var2 + 1
		    	if (i == Tri%SizeUp + 1) then
		        	XYZ(3,k) = Aux2Up(j,i,1)
		        	k = k + 1
		        	if (j == Tri%SizeUp) then
			        	XYZ(3,k) = Aux2Up(j,i,1)
			        	k = k + 1
		        	end if
		        end if
		    end do
		end do
		!Finally down triangles
		do j = 1, Tri%SizeDown
		    do i = j+1, Tri%SizeDown + 1
		        XYZ(3,var2) = Aux2Down(j,i,1)
		        var2 = var2 + 1
		    end do
		end do
		!If 3D
		if (relieve) then
			Z(1,:) = XYZ(3,:)
		else
			Z = 0d0
		end if
		!**Introduce data**
		!Coordinates
        E_IO = VTK_GEO(NN = Nn,X=XYZ(1,:),Y=XYZ(2,:),Z=Z(1,:))
        !Connectivity
        E_IO = VTK_CON(NC = Nf,connect = connect, cell_type = type)
        if (Celula) then
	        !Scalar values
	        E_IO = VTK_DAT(NC_NN = Nf,var_location = 'cell')!The other option is node
	        E_IO = VTK_VAR(NC_NN = Nf, varname = 'U_values', var = values)
        else
	        !Scalar values
	        E_IO = VTK_DAT(NC_NN = size(XYZ,2),var_location = 'node')!The other option is cell
	        E_IO = VTK_VAR(NC_NN = size(XYZ,2), varname = 'U_values', var = XYZ(3,:))
        end if

		!***************2x2 System Area***************
		!**Select data to plot**
		if(Sis==1) then
			a=2
			var2=1
			if(present(var)) var2 = var
			!**Select data to plot**
			select case (var2)
			    case (2)!F
				    do i = 2, Tri%SizeUp+1
				        do j = 1, i-1
				            Aux2Up(j,i,1) = Tri%UpUF(i,j,a)
				        end do
				    end do
				    do i =2, Tri%SizeDown+1
				    	do j = 1, i-1
				    	    Aux2Down(j,i,1) = Tri%DownUF(i,j,a)
				    	end do
				    end do
			    case (3)!Residual
			    	!Aux2Up = Tri%UpResAux
			    	!Aux2Down = Tri%DownResAux
				    do i = 2, Tri%SizeUp+1
				        do j = 1, i-1
				            Aux2Up(j,i,1) = Tri%UpResAux(j,i,a)
				        end do
				    end do
				    do i =2, Tri%SizeDown+1
				    	do j = 1, i-1
				    	    Aux2Down(j,i,1) = Tri%DownResAux(j,i,a)
				    	end do
				    end do
			    case (4)!Aux
				    do i = 2, Tri%SizeUp+1
				        do j = 1, i-1
				            Aux2Up(j,i,1) = Tri%UpResAux(i,j,a)
				        end do
				    end do
				    do i =2, Tri%SizeDown+1
				    	do j = 1, i-1
				    	    Aux2Down(j,i,1) = Tri%DownResAux(i,j,a)
				    	end do
				    end do
			    case default!U
				    do i = 2, Tri%SizeUp+1
				        do j = 1, i-1
				            Aux2Up(j,i,1) = Tri%UpUF(j,i,a)
				        end do
				    end do
				    do i =2, Tri%SizeDown+1
				    	do j = 1, i-1
				    	    Aux2Down(j,i,1) = Tri%DownUF(j,i,a)
				    	end do
				    end do
			end select
			!Introduce CC values as values in the faces of those triangles
			deallocate(values)
			allocate(values(Tri%SizeUp**2))
			!**Store data**
			!Up data
			pos = 1
			do j = 1, Tri%SizeUp
			    do i = j+1, Tri%SizeUp + 1
		            values(pos) = Aux2Up(j,i,1)
					pos = pos + 1
			    end do
			end do
			!down data
			do j = 1, Tri%SizeDown
				do i = j+1, Tri%SizeDown+1
				    values(pos) = Aux2Down(j,i,1)
				    pos = pos + 1
				end do
			end do

	        E_IO = VTK_DAT(NC_NN = size(values),var_location = 'cell')
	        E_IO = VTK_VAR(NC_NN = size(values), varname = 'U2_Values', var = values)
		end if

        !Close file
        E_IO = VTK_END()
	end subroutine Tri2VTK_CC
	!This subroutine creates an VTK archive to be read by ParaView
	!Only works for cell centered nodes
	!meshL(inout) :: type(mesh)(:)
	!(Optional) var(in) :: Integer. Select the variable to be plotted
	!								1 -> U(default)
	!								2 -> F
	!								3 -> Residual
	!								4 -> Auxiliar
	!(Optional) relief(in) :: Logical. Choose whether plot in 3d(default) or 2d
	!(Optional) filename(in) :: Character(*). Name of the archive to save into
	!(Optional) T(in) :: Integer. Number to add to the filename. Useful to time dependant problems
	!					or to see diferent iterations. As paraview permits to move along time.
	!(Optional) title(in) :: Character(*). Title to show in paraview
	!(Optional) CellColor :: logical. If true then the triangles will have a value for all the surface
	!								  if false(default) the data will be assigned to the nodes.
	!(Optional) Coarse :: logical. If true then only prints the coarse mesh
	!(Optional) Sys :: integer. If == 1 then we suppose that meshL contains two scalar values
	subroutine Tri2VTK_MeshCC(meshL, var, relief, filename,T , title, CellColor, coarse, Sys)
		Implicit none
		!Global variables
		type(mesh), intent(inout), dimension(:) :: meshL
		logical, optional, intent(in) :: relief, CellColor, coarse
		integer, optional, intent(in) :: var, T, Sys
		character (len = *), optional, intent(in) :: filename, title
		!Local variables
		integer :: i,j, E_IO, Nn, Nf, var2, malla, aux,a, sis, Mpos, NOUSAR
		logical :: relieve, Celula
		real, allocatable, dimension(:,:) :: XYZ
		integer, allocatable, dimension(:) :: connect, type
		double precision, allocatable, dimension(:) :: values
		double precision, allocatable, dimension(:,:,:) :: Aux2Up, Aux2Down
		character (len = 250) :: filename2, title2, cadena, path
		logical :: bol, coar
		!*****************Optional area*******************************
		filename2 = "Mesh_CC.vtk"
		if (present(filename)) filename2 = filename//".vtk"
		title2 = "Semi-structured triangle"
		if (present(title)) title2 = title
		relieve = .true.
		if (present(relief)) relieve = relief
		var2 = 1
		if (present(var)) var2 = var
		NOUSAR = var2
		Celula = .false.
		if (present(CellColor)) Celula = CellColor
		coar = .false.
		if (present(Coarse)) coar = Coarse
		sis = 0
		if (present(Sys)) sis = Sys
		!*************************************************************
		a = size(meshL(1)%Tri(1)%UpUF,3)
		!Number of nodes
		Nn = TotalNodes(meshL(1)%Tri(1),.true.) + meshL(1)%Tri(1)%SizeUp**2
		!Number of figures = Up + down triangles
		Nf = meshL(1)%Tri(1)%SizeUp**2
		!Allocate arrays
		allocate(XYZ(3,Nn*size(meshL)))
		allocate(type(Nf*size(meshL)))
		allocate(connect(Nf*5*size(meshL)))
		allocate(values(Nn*size(meshL)))
		allocate(Aux2Up(meshL(1)%Tri(1)%SizeUp+1,meshL(1)%Tri(1)%SizeUp+1,1))
		allocate(Aux2Down(0:meshL(1)%Tri(1)%SizeUp+1,0:meshL(1)%Tri(1)%SizeUp+1,1))
		!Prepare Data
		XYZ = 0.0
		values = 0d0
		!**Store type**
		!As we are with center nodes only the type is 10, tetrahedron
		Type = 10

		!Create file		!output_format = 'BINARY'
		if (present(T)) then
			i = len_trim(filename2) - 3
			write(cadena,*) T
			call removebksl(cadena)
			call insertstr(filename2,cadena,i)
!****************************ONLY WORKS WITH LINUX AND INTEL FORTRAN************************!
!			inquire( DIRECTORY='VTK_Images', exist=bol)
!			if (.not. bol) call system ('mkdir VTK_Images')
!****************************DISABLE IF NECESSARY, AS IT IS NOT CRITICAL*********************
			filename2 = "VTK_Images/"//filename2
		end if
		!Get absolute path
		call getCWD (path)
		!To avoid using debug as working directory
		call delsubstr(path,"/Debug")
		path = trim(path)//"/"//trim(filename2)
        E_IO = VTK_INI(output_format = 'ASCII',filename = path,&
        title = title2, mesh_topology = 'UNSTRUCTURED_GRID')

		do malla = 1, size(meshL)
			!**Select data to plot**
			select case (var2)
			    case (2)!F
				    do i = 2, meshL(malla)%Tri(1)%SizeUp+1
				        do j = 1, i-1
				            Aux2Up(j,i,1) = meshL(malla)%Tri(1)%UpUF(i,j,a)
				        end do
				    end do
				    do i =2, meshL(malla)%Tri(1)%SizeDown+1
				    	do j = 1, i-1
				    	    Aux2Down(j,i,1) = meshL(malla)%Tri(1)%DownUF(i,j,a)
				    	end do
				    end do
			    case (3)!Residual
			    	Aux2Up = meshL(malla)%Tri(1)%UpResAux

			    	do i = 2, meshL(malla)%Tri(1)%SizeDown+1
			    	    do j = 1, i-1
			    	        Aux2Down(j,i,1) = meshL(malla)%Tri(1)%DownResAux(j,i,a)
			    	    end do
			    	end do
			    	!Aux2Down = meshL(malla)%Tri(1)%DownResAux
			    case (4)!Aux
				    do i = 2, meshL(malla)%Tri(1)%SizeUp+1
				        do j = 1, i-1
				            Aux2Up(j,i,1) = meshL(malla)%Tri(1)%UpResAux(i,j,a)
				        end do
				    end do
				    do i =2, meshL(malla)%Tri(1)%SizeDown+1
				    	do j = 1, i-1
				    	    Aux2Down(j,i,1) = meshL(malla)%Tri(1)%DownResAux(i,j,a)
				    	end do
				    end do
			    case default!U
			        Aux2Up = meshL(malla)%Tri(1)%UpUF
			        Aux2Down = meshL(malla)%Tri(1)%DownUF
			end select
			!Start node for current triangle
			i = (malla - 1) * Nn + 1
			E_IO = (malla - 1) * Nf + 1
			!End node for current triangle
			j = Nn * malla
			aux = Nf * malla
			!Add one triangle to our file
			call DataCoordinates(meshL, malla, XYZ(:,i:j), values(i:j), relieve, Aux2Up,.true., Aux2Down, var2)
			!Create connectivity vector
			call Connectivity(meshL,malla,connect,.true.)
		end do

		!**Introduce data**
		!Coordinates
		values = XYZ(3,:)

		if (.not.relieve) XYZ(3,:) = 0d0
		!**Only Print coarse grid**
		if (coar) then
	        do malla = 1, size(meshL)
	        	!Start node for current triangle
				i = (malla - 1) * Nn + 1
				aux = (malla - 1) * 5 * meshL(malla)%Tri(1)%SizeUp**2 + 1
				!End node for current triangle
				j = Nn * malla
				XYZ(1,i:j)=meshL(malla)%Xp1(1)
				XYZ(2,i:j)=meshL(malla)%Xp1(2)
				XYZ(3,i:j)=0d0
				!Except one triangle
				XYZ(1,connect(aux+1)+1)=meshL(malla)%Xp1(1)
				XYZ(2,connect(aux+1)+1)=meshL(malla)%Xp1(2)

				XYZ(1,connect(aux+2)+1)=meshL(malla)%Xp2(1)
				XYZ(2,connect(aux+2)+1)=meshL(malla)%Xp2(2)

				XYZ(1,connect(aux+3)+1)=meshL(malla)%Xp3(1)
				XYZ(2,connect(aux+3)+1)=meshL(malla)%Xp3(2)

				XYZ(1,connect(aux+4)+1)=meshL(malla)%Xp3(1)
				XYZ(2,connect(aux+4)+1)=meshL(malla)%Xp3(2)
	        end do
		end if


        E_IO = VTK_GEO(NN = size(XYZ,2),X=XYZ(1,:),Y=XYZ(2,:),Z=XYZ(3,:))
        !Connectivity
        E_IO = VTK_CON(NC = size(type),connect = connect, cell_type = type)

        if (Celula) then
	        !Scalar values
	        E_IO = VTK_DAT(NC_NN = size(values),var_location = 'cell')!The other option is node
	        E_IO = VTK_VAR(NC_NN = size(values), varname = 'U_values', var = values)
        else
        	if (present(Sys)) then
		        !Scalar values
		        E_IO = VTK_DAT(NC_NN = size(XYZ,2),var_location = 'node')!The other option is cell
		        E_IO = VTK_VAR(NC_NN = size(XYZ,2), varname = 'Scalar_2', var = values)
        	else
		        !Scalar values
		        E_IO = VTK_DAT(NC_NN = size(XYZ,2),var_location = 'node')!The other option is cell
		        E_IO = VTK_VAR(NC_NN = size(XYZ,2), varname = 'U_values', var = values)
        	end if
        end if
		!********************Add more data**********************
		!!!!You cannot store two values consecutively, of the same type: cell, node!!!!
	    !If Sys then save the second data
		if (Sis==1) then
			a = 1
			deallocate(values)
			var2 = meshL(1)%Tri(1)%SizeUp**2
			allocate(values(size(meshL)*var2))

	        do Mpos = 1, size(meshL)

				!**Select data to plot**
				select case (NOUSAR)
				    case (2)!F
					    do i = 2, meshL(Mpos)%Tri(1)%SizeUp+1
					        do j = 1, i-1
					            Aux2Up(j,i,a) = meshL(Mpos)%Tri(1)%UpUF(i,j,a)
					        end do
					    end do
					    do i =2, meshL(Mpos)%Tri(1)%SizeDown+1
					    	do j = 1, i-1
					    	    Aux2Down(j,i,a) = meshL(Mpos)%Tri(1)%DownUF(i,j,a)
					    	end do
					    end do
				    case (3)!Residual
					    do i = 2, meshL(Mpos)%Tri(1)%SizeUp+1
					        do j = 1, i-1
					            Aux2Up(j,i,a) = meshL(Mpos)%Tri(1)%UpResAux(j,i,a)
					        end do
					    end do
					    do i =2, meshL(Mpos)%Tri(1)%SizeDown+1
					    	do j = 1, i-1
					    	    Aux2Down(j,i,a) = meshL(Mpos)%Tri(1)%DownResAux(j,i,a)
					    	end do
					    end do
				    case (4)!Aux
					    do i = 2, meshL(Mpos)%Tri(1)%SizeUp+1
					        do j = 1, i-1
					            Aux2Up(j,i,a) = meshL(Mpos)%Tri(1)%UpResAux(i,j,a)
					        end do
					    end do
					    do i =2, meshL(Mpos)%Tri(1)%SizeDown+1
					    	do j = 1, i-1
					    	    Aux2Down(j,i,a) = meshL(Mpos)%Tri(1)%DownResAux(i,j,a)
					    	end do
					    end do
				    case default!U
					    do i = 2, meshL(Mpos)%Tri(1)%SizeUp+1
					        do j = 1, i-1
					            Aux2Up(j,i,a) = meshL(Mpos)%Tri(1)%UpUF(j,i,a)
					        end do
					    end do
					    do i =2, meshL(Mpos)%Tri(1)%SizeDown+1
					    	do j = 1, i-1
					    	    Aux2Down(j,i,a) = meshL(Mpos)%Tri(1)%DownUF(j,i,a)
					    	end do
					    end do
				end select

	        	!Start node for current triangle
				aux = (Mpos - 1) * var2 + 1
				!**Store data**
				!Up data
				do j = 1, meshL(Mpos)%Tri(1)%SizeUp
				    do i = j+1, meshL(Mpos)%Tri(1)%SizeUp + 1
			            values(aux) = Aux2Up(j,i,a)
						aux = aux + 1
				    end do
				end do
				!down data
				do j = 1, meshL(Mpos)%Tri(1)%SizeDown
					do i = j+1, meshL(Mpos)%Tri(1)%SizeDown+1
					    values(aux) = Aux2Down(j,i,a)
					    aux = aux + 1
					end do
				end do

			end do
	        E_IO = VTK_DAT(NC_NN = size(values),var_location = 'cell')
	        E_IO = VTK_VAR(NC_NN = size(values), varname = 'Scalar_1', var = values)

			!K coefficient used in each triangle
	        do malla = 1, size(meshL)
	        	!Start node for current triangle
				i = (malla - 1) * Nn + 1
				!End node for current triangle
				j = Nn * malla
				XYZ(2,i:j)=meshL(malla)%K_coef
	        end do
	        E_IO = VTK_DAT(NC_NN = size(XYZ,2),var_location = 'node')!The other option is cell
	        E_IO = VTK_VAR(NC_NN = size(XYZ,2), varname = 'K_Coefficient', var = int(XYZ(2,:)))

		end if

		if (size(meshL(1)%Tri)<=5) then
		    !Store smoothers
			deallocate(values)
			aux = meshL(1)%Tri(1)%SizeUp**2
			allocate(values(size(meshL)*aux))
			!Smoothers used in each triangle
	        do malla = 1, size(meshL)
	        	!Start node for current triangle
				i = (malla - 1) * aux + 1
				!End node for current triangle
				j = aux * malla
		        values(i:j) = int(meshL(malla)%method/10)*10
	        end do
	        E_IO = VTK_DAT(NC_NN = size(values),var_location = 'cell')
	        E_IO = VTK_VAR(NC_NN = size(values), varname = 'Smoother', var = int(values))
	        if (Sis/=1) then
	        	!K coefficient used in each triangle
		        do malla = 1, size(meshL)
		        	!Start node for current triangle
					i = (malla - 1) * Nn + 1
					!End node for current triangle
					j = Nn * malla
					XYZ(2,i:j)=meshL(malla)%K_coef
		        end do
		        E_IO = VTK_DAT(NC_NN = size(XYZ,2),var_location = 'node')!The other option is cell
		        E_IO = VTK_VAR(NC_NN = size(XYZ,2), varname = 'K_Coefficient', var = int(XYZ(2,:)))
	        end if
		end if
        !Close file
        E_IO = VTK_END()
	end subroutine

	!Introduce the coordinates for a given triangle
	!meshL(inout) :: type(mesh)(:)
	!Mpos(in) :: Integer. Current position
	!XYZ(inout) :: Real(3,:). Position of the coordinates
	!Values(inout) :: Double(:). Value in the voronoi centers
	!Relieve(in) :: Logical. Select whether 3d plot or no
	!Aux2Up(in) :: Double(:,:,:). Contains the data that will be introduced in XYZ and values
	!CheckCC(in) :: logical
	!(Optional, but necessary for CC)Aux2Down(in) :: Double(:,:,:). Contains the data that will be introduced in XYZ and values
	!(Optional)var(in) :: Integer. Select the variable to be plotted
	!								1 -> U(default)
	!								2 -> F
	!								3 -> Residual
	!								4 -> Auxiliar
	subroutine DataCoordinates(meshL,Mpos,XYZ, values, relieve, Aux2Up, CheckCC, Aux2Down, var)
		Implicit none
		!Global variables
		type(mesh), intent(inout), dimension(:) :: meshL
		logical, optional, intent(in) :: relieve
		real, intent(inout), dimension(:,:) :: XYZ
		double precision,intent(inout), dimension(:) :: values
		double precision,intent(in), dimension(:,:,:) :: Aux2Up
		double precision, optional,  intent(in), dimension(0:,0:,:) :: Aux2Down
		integer, intent(in) :: Mpos, var
		logical, intent(in) :: CheckCC
		!Local variables
		integer :: i,j,k, pos, var2, m, Mpos2,cas
		integer, allocatable, dimension(:) :: Mprev
		double precision, dimension(2) :: v1,v2, aux, point, Xvor
		logical :: relief

		relief = .true.
		if (present(relieve)) relief = relieve

		!As maximum the number of triangles with a common vertex can be the total amount
		!of triangles plus three, in order to avoid that it will try to check three triangles
		!but there are only 2 triangles
		allocate(Mprev(size(meshL)+3))

		if (checkCC) then
			!**Position of the nodes**
			i = 1 + meshL(Mpos)%Tri(1)%SizeUP
			!Store coordinates
			v1 = (meshL(Mpos)%Tri(1)%Xp2-meshL(Mpos)%Tri(1)%Xp1) / (meshL(Mpos)%Tri(1)%SizeUp)
			v2 = (meshL(Mpos)%Tri(1)%Xp3-meshL(Mpos)%Tri(1)%Xp1) / (meshL(Mpos)%Tri(1)%SizeUp)
			point = meshL(Mpos)%Tri(1)%Xp1
			k = 1
			do i = i, 1, -1
				aux = point
				do j = 1, i
					XYZ(1,k) = point(1)
					XYZ(2,k) = point(2)
					!2D
					XYZ(3,k) = 0d0
					k = k + 1
					point = point + v1
				end do
				point = aux + v2
			end do
			!Position of the voronoi center for up triangles
			point = meshL(Mpos)%Tri(1)%Xp1
			do j = 1, meshL(Mpos)%Tri(1)%SizeUP
				aux = point
				do i = j+1, meshL(Mpos)%Tri(1)%SizeUP+1
					call getVoronoiCenter(point,point+v1,point+v2,Xvor)
					if (relief) then
					    XYZ(1,k) = Xvor(1)
						XYZ(2,k) = Xvor(2)
					else
					 	XYZ(1,k) = point(1)
						XYZ(2,k) = point(2)
					end if
					!2D
					XYZ(3,k) = 0d0
					k = k + 1
					point = point + v1
				end do
				point = aux + v2
			end do
			!Position of the voronoi center for down triangles
			point = meshL(Mpos)%Tri(1)%Xp1 + v1
			do j = 1, meshL(Mpos)%Tri(1)%SizeDown
				aux = point
				do i = j+1, meshL(Mpos)%Tri(1)%SizeDown+1
					call getVoronoiCenter(point,point+v2,point+v2-v1,Xvor)
					if (relief) then
					    XYZ(1,k) = Xvor(1)
						XYZ(2,k) = Xvor(2)
					else
					 	XYZ(1,k) = point(1)
						XYZ(2,k) = point(2)
					end if
					!2D
					XYZ(3,k) = 0d0
					k = k + 1
					point = point + v1
				end do
				point = aux + v2
			end do
			!**Store data**
			!Up data
			pos = 1
			do j = 1, meshL(Mpos)%Tri(1)%SizeUp
			    do i = j+1, meshL(Mpos)%Tri(1)%SizeUp + 1
		            values(pos) = Aux2Up(j,i,1)
					pos = pos + 1
			    end do
			end do
			!down data
			do j = 1, meshL(Mpos)%Tri(1)%SizeDown
				do i = j+1, meshL(Mpos)%Tri(1)%SizeDown+1
				    values(pos) = Aux2Down(j,i,1)
				    pos = pos + 1
				end do
			end do

			!Store 3D data
			k = 1
			var2 = TotalNodes(meshL(Mpos)%Tri(1),CheckCC) + 1

			do j = 1, meshL(Mpos)%Tri(1)%SizeUp
			    do i = j+1, meshL(Mpos)%Tri(1)%SizeUp + 1
			        cas = 0
			        !For boundary values we will put the same data for both sides, the average
			        if (j+1==i.or.j==1.or.i==meshL(Mpos)%Tri(1)%SizeUp+1) then!Boundary side
			            !**If vertex**
			            if ((j==1.and.i==2).or.(j==1.and.i==meshL(Mpos)%Tri(1)%SizeUp+1)&
			            .or.(j==meshL(Mpos)%Tri(1)%SizeUp.and. i == meshL(Mpos)%Tri(1)%SizeUp+1)) then

			                Mpos2 = Mpos
			                Mprev = 0d0
			                XYZ(3,var2) = vertexAverage(meshL, j, i, Mpos2,1, var)
!			                XYZ(3,var2) = vertexAverage(meshL, j, i, Mpos2, Mprev,  var, .true.)
			                XYZ(3,k) = XYZ(3,var2)
			                k = k + 1
			                var2 = var2 + 1
			                !Just to get the case number
			                m = getNeigValue(meshL, Mpos, j,i, var,cas)
			            else!If side
			                !print *,Mpos,j,i,var
			                !print *, XYZ(3,k),XYZ(3,k)*2d0-Aux2Up(j,i,1), Aux2Up(j,i,1)
			                XYZ(3,var2) = (getNeigValue(meshL, Mpos, j,i, var,cas) + Aux2Up(j,i,1)) / 2d0
			                XYZ(3,k) = XYZ(3,var2)
			                k = k + 1
			                var2 = var2 + 1

			            end if
			        else
					    !Vertex node
			            XYZ(3,k) = (Aux2Up(j,i,1)+Aux2Down(j,i-1,1))/2d0
			            k = k + 1
			            !Voronoi node
			            XYZ(3,var2) = Aux2Up(j,i,1)
			            var2 = var2 + 1
			        end if
				    			        !****************************************
			        !***Extra points or special situations***
			        !For l1 side
			        select case (cas)
			            case (10)!l1->l1
			                if (j/=1.and.(meshL(Mpos)%Neig(1)>Mpos)) then
			                    XYZ(3,k-1) = XYZ(3,var2-1-(meshL(Mpos)%Tri(1)%SizeUp-j+2))
			                end if
			            case (12)
			                if (meshL(Mpos)%Neig(1)>Mpos.and.j/=1) then
			                    XYZ(3,k-1) = XYZ(3,var2-1-(meshL(Mpos)%Tri(1)%SizeUp-j+2))
			                end if
			            case (13)
			                if (i/=2) then
			                    XYZ(3,k-1) = XYZ(3,var2-2)
			                end if
			            case (14)!l2->l2
			                if (meshL(Mpos)%Neig(2)>Mpos.and.i/=2) then
			                    XYZ(3,k-1) = XYZ(3,var2-2)
			                end if
			            case (15)
			                if (meshL(Mpos)%Neig(2)>Mpos.and.i/=2) then
			                    XYZ(3,k-1) = XYZ(3,var2-2)
			                end if
			        end select

!						!Modify the vertex
!						!Xp1
!						if (k==meshL(Mpos)%Tri(1)%SizeUp + 2) then
!							 XYZ(3,k-1) = XYZ(3,TotalNodes(meshL(Mpos)%Tri(1),CheckCC) + 1)
!						end if
!						!Xp2
!						if (k==2*meshL(Mpos)%Tri(1)%SizeUp + 1) then
!							 XYZ(3,k-1) = XYZ(3,TotalNodes(meshL(Mpos)%Tri(1),CheckCC) + meshL(Mpos)%Tri(1)%SizeUp)
!						end if

				        			        !For L3 and Xp3 we also have to add more points
			        if (i == meshL(Mpos)%Tri(1)%SizeUp + 1) then
			            if (meshL(Mpos)%Neig(3)==0 .or. meshL(Mpos)%Dir(3).or.j==1.or.&
			            meshL(Mpos)%Neig(3)<Mpos) then
			                XYZ(3,k) = XYZ(3,var2-1)
			            else
			                XYZ(3,k) = XYZ(3,var2-1-(meshL(Mpos)%Tri(1)%SizeUp-j+1))
			            end if
			            k = k + 1
			            if (j == meshL(Mpos)%Tri(1)%SizeUp) then
			                !				        		XYZ(3,k-1) = XYZ(3,var2-1)
			                XYZ(3,k) = XYZ(3,var2-1)
			                k = k + 1
			            end if
			        end if
			    !************************************************
			    end do
			end do
			!Finally down triangles
			do j = 1, meshL(Mpos)%Tri(1)%SizeDown
			    do i = j+1, meshL(Mpos)%Tri(1)%SizeDown + 1
			        XYZ(3,var2) = Aux2Down(j,i,1)
			        var2 = var2 + 1
			    end do
			end do
		else
		!*****************SD area*******************
			!**Position of the nodes**
			!Initial quantity of nodes in the L2 side
			i = 1 + meshL(Mpos)%Tri(1)%SizeUP * 2
			!Store coordinates
			v1 = (meshL(Mpos)%Tri(1)%Xp2-meshL(Mpos)%Tri(1)%Xp1) / (2d0*meshL(Mpos)%Tri(1)%SizeUp)
			v2 = (meshL(Mpos)%Tri(1)%Xp3-meshL(Mpos)%Tri(1)%Xp1) / (2d0*meshL(Mpos)%Tri(1)%SizeUp)
			point = meshL(Mpos)%Tri(1)%Xp1
			k = 1
			do i = i, 1, -1
				aux = point
				do j = 1, i
					XYZ(1,k) = point(1)
					XYZ(2,k) = point(2)
					!2D
					XYZ(3,k) = 0d0
					k = k + 1
					point = point + v1
				end do
				point = aux + v2
			end do

			!**Store data**
			pos = 1
			k = meshL(Mpos)%Tri(1)%SizeUp * 2 + 1
			do j = 1, meshL(Mpos)%Tri(1)%SizeUp
			    do i = j+1, meshL(Mpos)%Tri(1)%SizeUp + 1
			        if (i == j + 1) then
			            values(pos) = Aux2Up(j,i,2)
			            values(pos+1) = Aux2Up(j,i,2)
			            values(pos+k) = Aux2Up(j,i,1)
			            values(pos+k+1) = Aux2Up(j,i,3)
			        else if (i == meshL(Mpos)%Tri(1)%SizeUp + 1) then
			            values(pos) = (Aux2Up(j,i,2) + Aux2Up(j,i-1,2)) / 2d0
			            values(pos+1) = Aux2Up(j,i,2)
			            values(pos+k) = Aux2Up(j,i,1)
			            values(pos+k+1) = Aux2Up(j,i,3)
			            values(pos+2) = Aux2Up(j,i,2)
			        else
			            values(pos) = (Aux2Up(j,i,2) + Aux2Up(j,i-1,2)) / 2d0
			            values(pos+1) = Aux2Up(j,i,2)
			            values(pos+k) = Aux2Up(j,i,1)
			            values(pos+k+1) = Aux2Up(j,i,3)

			        end if
				pos = pos + 2
			    end do
			pos = pos + k
			k = k - 2
			end do
			!manually 2 nodes from Xp3
			values(size(values,1)-3) = Aux2Up(meshL(Mpos)%Tri(1)%SizeUp,meshL(Mpos)%Tri(1)%SizeUp+1,2)
			values(size(values,1)) = Aux2Up(meshL(Mpos)%Tri(1)%SizeUp,meshL(Mpos)%Tri(1)%SizeUp+1,2)
		end if
	end subroutine DataCoordinates

	!Creates connectivity matrix for a group of triangles
	!meshL(inout) :: type(mesh)(:)
	!Mpos(in) :: Integer. Current position inside the mesh
	!Connect(inout) :: integer(:). Connectivity vector
	!CheckCC(in) :: logical
	subroutine Connectivity(meshL,Mpos,connect, checkCC)
		Implicit none
		!Global variables
		type(mesh), intent(inout), dimension(:) :: meshL
		integer, intent(in) :: Mpos
		integer, intent(inout), dimension(:) :: connect
		logical, intent(in) :: CheckCC
		!Local variables
		integer :: i,j,k, pos, var2, initial, fin, col, row, aux


		if (checkCC) then
			!Prepare data
			!Start node for current triangle
			Initial = (Mpos - 1) * 5 * meshL(Mpos)%Tri(1)%SizeUp**2
			!End node for current triangle
			fin = Initial + getUpTriangles(meshL(1)%Tri(1))*5
			!**Connectivity**
			!For up triangles
			aux = (Mpos - 1) * (meshL(Mpos)%Tri(1)%SizeUp**2 + getUpTriangles(meshL(Mpos)%Tri(1)) + meshL(Mpos)%Tri(1)%SizeUp + 1)
			pos =  aux + 1
			Var2 = aux + TotalNodes(meshL(Mpos)%Tri(1),CheckCC)
			j = 2
			col = 1
			row = 1
			k = meshL(Mpos)%Tri(1)%SizeUp + 1
			do i = Initial + 1, fin, 5 !5 points is a triangle. Therefore, one triangle per iteration
			    connect(i) = 4
				connect(i+1) = pos - 1
				connect(i+2) = pos
				connect(i+3) = pos - 1 + k
				connect(i+4) = Var2
				Var2 = Var2 + 1
				if (j >= k) then
					k = k - 1
					j = 2
					pos = pos + 2
					row = row + 1
				else
					j = j + 1
					pos = pos + 1
					col = col + 1
				end if
			end do
			!For Down triangles
			j = Initial
			k = meshL(Mpos)%Tri(1)%SizeDown + 1
			Initial = fin + 1
			fin = j + 5 * meshL(Mpos)%Tri(1)%SizeUp**2
			pos = aux + 2
			j = 2
			do i = Initial, fin, 5
				connect(i) = 4
				connect(i+1) = pos - 1
				connect(i+2) = pos + k
				connect(i+3) = pos + (k - 1)
				connect(i+4) = Var2
				Var2 = Var2 + 1
				if (j >= k) then
					k = k - 1
					j = 2
					pos = pos + 3
				else
					j = j + 1
					pos = pos + 1
				end if
			end do

		else!SD part
			!For up triangles
			!Prepare data
			!Start node for current triangle
			Initial = (Mpos - 1) * 7 * meshL(Mpos)%Tri(1)%SizeUp**2
			!End node for current triangle
			fin = Initial + getUpTriangles(meshL(1)%Tri(1))*7
			!**Connectivity**
			!For up triangles
			aux = (Mpos - 1) * (TotalNodes(meshL(Mpos)%Tri(1),CheckCC))
			pos =  aux + 1
			j = 3
			k = meshL(Mpos)%Tri(1)%SizeUp * 2 + 1
			do i = Initial+1 , fin, 7 !7 points is a triangle. Therefore, one triangle per iteration
			    connect(i) = 6
				connect(i+1) = pos - 1
				connect(i+2) = pos + 1
				connect(i+3) = pos + 2 * ( k - 1)
				connect(i+4) = pos
				connect(i+5) = pos + k
				connect(i+6) = pos + k - 1

				if (j >= k) then
					k = k - 2
					j = 3
					pos = pos + k + 4
				else
					j = j + 2
					pos = pos + 2
				end if
			end do
			!For Down triangles
			j = Initial
			k = meshL(Mpos)%Tri(1)%SizeDown + 1
			Initial = fin + 1
			fin = j + 7 * meshL(Mpos)%Tri(1)%SizeUp**2
			pos = aux + 3
			j = 3
			k = meshL(Mpos)%Tri(1)%SizeUp * 2 + 1
			do i = Initial, fin , 7
				connect(i) = 6
				connect(i+1) = pos - 1
				connect(i+2) = pos + 2 * (k - 1)
				connect(i+3) = pos + 2 * (k - 2)
				connect(i+4) = pos - 1 + k
				connect(i+5) = pos + 2 * (k - 1) - 1
				connect(i+6) = pos - 2 + k

				if (j >= k-2) then
					k = k - 2
					j = 3
					pos = pos + k + 6
				else
					j = j + 2
					pos = pos + 2
				end if
			end do
		end if

	end subroutine Connectivity

	!Depending on the input the result is:
	!-> Number of nodes for a quadratic type triangle(22 for VTK)
	!	This means nodes in the vextex and in the middle of the sides
	!-> Number of nodes for simple triangle(5 for VTK)
	!	Nodes in the vertex only. For cell centered nodes
	!
	!Tri(in) :: type(Triangle)
	!CheckCC(in) :: logical
	!Return  :: Integer
	integer function TotalNodes(Tri, checkCC)
	    Implicit none
	    !Global variables
		type(Triangle), intent(in) :: Tri
		logical, intent(in) :: CheckCC
		!Local variables
		integer :: i,j

		j = 0
		do i = 1, Tri%SizeUp+1
		    j = j + i
		end do

		if (checkCC) then
			TotalNodes = getUpTriangles(Tri) + Tri%SizeUp + 1
		else
			TotalNodes = j + 3 * getUpTriangles(Tri)
		end if
	end function TotalNodes

	!Returns the voronoi data of a neighbour triangle, if there is no neighbour then returns the value of that side
	!This is to allow to produce a code in which it doesn't matter if that side has a neighbour or not. So you can do the average easly.
	!getNeigValue(meshL, Mpos, j,i, var)
	function getNeigValue(meshL, Mpos, j,i, var,cas)
		Implicit none
		!Global variables
		type(mesh), intent(inout), dimension(:) :: meshL
		integer, intent(in) :: Mpos, var, j, i
		integer, optional, intent(inout) :: cas
		!Local variables
		double precision :: getNeigValue
		integer :: aux, jj, ii, m, malla, aux2
!		double precision, allocatable, dimension(:,:,:) :: Aux2Up
!		allocate(Aux2Up(meshL(Mpos)%Tri(1)%SizeUp+1,meshL(Mpos)%Tri(1)%SizeUp+1,1))
		aux = MeshL(Mpos)%Tri(1)%SizeUp
		aux2 = 0

		m = 0
		if (j+1==i) m = 1

		if ((j==1 .or. (MeshL(Mpos)%Neig(1) == 0 .and.j==1)).and.MeshL(Mpos)%Neig(2)/= 0) then
			m = 2
		end if

		if ((i==aux+1 .or. (MeshL(Mpos)%Neig(1) == 0 .and. i == aux+1) &
			 .or. (MeshL(Mpos)%Neig(2) == 0 .and. i == aux+1)).and.MeshL(Mpos)%Neig(3) /= 0 ) then
			m = 3
		end if

		if (m /= 0) then
        	if (MeshL(Mpos)%Neig(m) /= 0) then
				malla = MeshL(Mpos)%Neig(m)
!				select case (var)
!				    case (2)!F
!					    do ii = 2, meshL(malla)%Tri(1)%SizeUp+1
!					        do jj = 1, ii-1
!					            Aux2Up(jj,ii,1) = meshL(malla)%Tri(1)%UpUF(ii,jj,1)
!					        end do
!					    end do
!				    case (3)!Residual
!				    	Aux2Up = meshL(malla)%Tri(1)%UpResAux
!				    case (4)!Aux
!					    do ii = 2, meshL(malla)%Tri(1)%SizeUp+1
!					        do jj = 1, ii-1
!					            Aux2Up(jj,ii,1) = meshL(malla)%Tri(1)%UpResAux(ii,jj,1)
!					        end do
!					    end do
!				    case default!U
!				        Aux2Up = meshL(malla)%Tri(1)%UpUF
!				end select
				call getEqvCoord(j,i,jj,ii,m,NumLoc(MeshL( malla )%Neig, Mpos ),&
						aux, MeshL(Mpos)%Dir(m),aux2)
					getNeigValue = getData(meshL, malla, var,jj,ii,1)!Aux2Up(jj,ii,1)
			else
				getNeigValue =getData(meshL, Mpos, var,j,i,1)!meshL(Mpos)%Tri(1)%UpUF(j,i,1)
	        end if
	    else
	    	getNeigValue = getData(meshL, Mpos, var,j,i,1)!meshL(Mpos)%Tri(1)%UpUF(j,i,1)
		end if
	!Return case situation of getEqvCoord
	if (present(cas)) cas = aux2
	end function getNeigValue

	!Returns the data for the given meshl.
	!Only work for up triangles
	!And it is thought to be used only in this module
	!
	!double -> getData(meshL, Mpos, var,j,i,node)
	!
	!var(in) :: Integer. Select the variable
	!								1 -> U(default)
	!								2 -> F
	!								3 -> Residual
	!								4 -> Auxiliar
	function getData(meshL, Mpos, var,j,i,node)
		Implicit none
		!Global variables
		type(mesh), intent(inout), dimension(:) :: meshL
		integer, intent(in) :: Mpos, var,j,i,node
		!Local variables
		double precision :: getData
		select case (var)
		    case (2)
				getData = meshL(Mpos)%Tri(1)%UpUF(i,j,node)
		    case (3)
				getData = meshL(Mpos)%Tri(1)%UpResAux(j,i,node)
		    case (4)
				getData = meshL(Mpos)%Tri(1)%UpResAux(i,j,node)
		    case default
				getData = meshL(Mpos)%Tri(1)%UpUF(j,i,node)
		end select

	end function getData

	!The result is an arithmetic average between all the values that have a common vertex
	!meshL(inout) :: type(meshL)(:)
	!j(in) :: integer.
	!i(in) :: integer.
	!Mpos(inout) :: integer. Position of the present triangle
	!level(in) :: integer(:).
	!var(in) :: integer. Type of variable.
	function vertexAverage(meshL, j, i, Mpos,level, var)
		Implicit none
		!Global variables
		type(mesh), intent(inout), dimension(:) :: meshL
		integer, intent(in) :: j,i, level, var
		integer, intent(inout) :: Mpos
!		double precision,intent(in), dimension(:,:,:) :: Aux2Up
		!Local variables
		double precision :: vertexAverage
		integer :: ii, tam, aux, jj, pos, a
		logical :: done
		integer, dimension(size(meshL)+1,4) :: res
		integer, dimension(size(meshL)+1) :: Mprev
		!Prepare data
		a = size(meshL(1)%Tri(1)%UpUF,3)
		call getCollideVertexes(meshL, j, i,level, Mpos, Mprev, .true., pos, res,.false.)
			aux = 1
			vertexAverage = 0d0

			do while (res(aux,1)/=0)
			    vertexAverage = vertexAverage + getData(meshL, res(aux,1), var,res(aux,2),res(aux,3),a)
			    aux=aux+1
			end do
				vertexAverage = vertexAverage /dble(aux-1)
			!If there are no neighbours
			if (aux==1) then
				vertexAverage =  getData(meshL, Mpos, var,j,i,a)
			end if
	end function vertexAverage
!
!	!OLD VERTEX AVERAGE SUBROUTINE, DOES NOT NEED SMOOTHERS.F90
!	!This is a recursive subroutine that iterates through all the collide triangles with a vertex that have the same coordinates
!	!The result is an arithmetic average between all the values
!	!meshL(inout) :: type(meshL)(:)
!	!j(in) :: integer.
!	!i(in) :: integer.
!	!Mpos(inout) :: integer. Position of the present triangle
!	!Mprev(inout) :: integer(:). List that stores all the triangles that have been already visited
!	!Original(in) :: Logical. This must be allways true. Is to identify the first call.
!	recursive function vertexAverage(meshL, j, i, Mpos, Mprev,  var, original) result (Average)
!		Implicit none
!		!Global variables
!		type(mesh), intent(inout), dimension(:) :: meshL
!		integer, intent(in) :: j,i, var
!		integer, intent(inout) :: Mpos
!		integer, intent(inout), dimension(:) :: Mprev
!		logical, intent(in) :: original
!		!Local variables
!		double precision :: Average
!		integer :: ii, tam, aux, jj
!		logical :: done, added
!
!		!prepare data
!		done = .false.
!		added = .false.
!		tam = MeshL(1)%Tri(1)%SizeUp
!
!		average = 0d0
!
!		if ((j==1.and.i==2).and..not.CheckVector(Mprev,Mpos) ) then
!		    do ii = 1, 3
!		    	aux = meshL(Mpos)%Neig(ii)
!				if (AreEqual(meshL(Mpos)%Xp1, meshL(aux)%Xp1).and..not.CheckVector(Mprev,aux)) then
!					if (.not.added) then
!						average = average + getData(meshL, Mpos, var,j,i,1)
!						do jj = 1, size(Mprev)
!						    if (Mprev(jj)==0) then
!						    	Mprev(jj) = Mpos
!						    	added = .true.
!						    	exit
!						    end if
!						end do
!					end if
!
!					average = average + vertexAverage(meshL,1,2,aux,Mprev, var, .false.)
!					done = .true.
!				else if (AreEqual(meshL(Mpos)%Xp1, meshL(aux)%Xp2).and..not.CheckVector(Mprev,aux)) then
!					if (.not.added) then
!						average = average + getData(meshL, Mpos, var,j,i,1)
!						do jj = 1, size(Mprev)
!						    if (Mprev(jj)==0) then
!						    	Mprev(jj) = Mpos
!						    	added = .true.
!						    	exit
!						    end if
!						end do
!					end if
!					average = average + vertexAverage(meshL,1,tam + 1,aux,Mprev, var, .false.)
!					done = .true.
!				else if (AreEqual(meshL(Mpos)%Xp1, meshL(aux)%Xp3).and..not.CheckVector(Mprev,aux)) then
!					if (.not.added) then
!						average = average + getData(meshL, Mpos, var,j,i,1)
!						do jj = 1, size(Mprev)
!						    if (Mprev(jj)==0) then
!						    	Mprev(jj) = Mpos
!						    	added = .true.
!						    	exit
!						    end if
!						end do
!					end if
!					average = average + vertexAverage(meshL,tam,tam + 1,aux,Mprev, var, .false.)
!					done = .true.
!				end if
!		    end do
!		end if
!
!		if ((j==1.and.i==tam+1).and..not.CheckVector(Mprev,Mpos) ) then
!		    do ii = 1, 3
!		    	aux = meshL(Mpos)%Neig(ii)
!				if (AreEqual(meshL(Mpos)%Xp2, meshL(aux)%Xp1).and..not.CheckVector(Mprev,aux)) then
!					if (.not.added) then
!						average = average + getData(meshL, Mpos, var,j,i,1)
!						do jj = 1, size(Mprev)
!						    if (Mprev(jj)==0) then
!						    	Mprev(jj) = Mpos
!						    	added = .true.
!						    	exit
!						    end if
!						end do
!					end if
!					average = average + vertexAverage(meshL,1,2,aux,Mprev, var, .false.)
!					done = .true.
!				else if (AreEqual(meshL(Mpos)%Xp2, meshL(aux)%Xp2).and..not.CheckVector(Mprev,aux)) then
!					if (.not.added) then
!						average = average + getData(meshL, Mpos, var,j,i,1)
!						do jj = 1, size(Mprev)
!						    if (Mprev(jj)==0) then
!						    	Mprev(jj) = Mpos
!						    	added = .true.
!						    	exit
!						    end if
!						end do
!					end if
!					average = average + vertexAverage(meshL,1,tam + 1,aux,Mprev, var, .false.)
!					done = .true.
!				else if (AreEqual(meshL(Mpos)%Xp2, meshL(aux)%Xp3).and..not.CheckVector(Mprev,aux)) then
!					if (.not.added) then
!						average = average + getData(meshL, Mpos, var,j,i,1)
!						do jj = 1, size(Mprev)
!						    if (Mprev(jj)==0) then
!						    	Mprev(jj) = Mpos
!						    	added = .true.
!						    	exit
!						    end if
!						end do
!					end if
!					average = average + vertexAverage(meshL,tam,tam + 1,aux,Mprev, var, .false.)
!					done = .true.
!				end if
!		    end do
!		end if
!
!		if ((j==tam.and.i==tam+1).and..not.CheckVector(Mprev,Mpos) ) then
!		    do ii = 1, 3
!		    	aux = meshL(Mpos)%Neig(ii)
!				if (AreEqual(meshL(Mpos)%Xp3, meshL(aux)%Xp1) .and..not.CheckVector(Mprev,aux)) then
!					if (.not.added) then
!						average = average + getData(meshL, Mpos, var,j,i,1)
!						do jj = 1, size(Mprev)
!						    if (Mprev(jj)==0) then
!						    	Mprev(jj) = Mpos
!						    	added = .true.
!						    	exit
!						    end if
!						end do
!					end if
!					average = average + vertexAverage(meshL,1,2,aux,Mprev, var, .false.)
!					done = .true.
!				else if (AreEqual(meshL(Mpos)%Xp3, meshL(aux)%Xp2).and..not.CheckVector(Mprev,aux)) then
!					if (.not.added) then
!						average = average + getData(meshL, Mpos, var,j,i,1)
!						do jj = 1, size(Mprev)
!						    if (Mprev(jj)==0) then
!						    	Mprev(jj) = Mpos
!						    	added = .true.
!						    	exit
!						    end if
!						end do
!					end if
!					average = average + vertexAverage(meshL,1,tam + 1,aux,Mprev, var, .false.)
!					done = .true.
!				else if (AreEqual(meshL(Mpos)%Xp3, meshL(aux)%Xp3).and..not.CheckVector(Mprev,aux)) then
!					if (.not.added) then
!						average = average + getData(meshL, Mpos, var,j,i,1)
!						do jj = 1, size(Mprev)
!						    if (Mprev(jj)==0) then
!						    	Mprev(jj) = Mpos
!						    	added = .true.
!						    	exit
!						    end if
!						end do
!					end if
!					average = average + vertexAverage(meshL,tam,tam + 1,aux,Mprev, var, .false.)
!					done = .true.
!				end if
!		    end do
!		end if
!
!		if (original) then
!			jj = 1
!			ii = 0
!			do while(jj /= 0)
!				ii = ii + 1
!			    jj = Mprev(ii)
!			end do
!			ii = ii - 1
!			if (ii<=0) then
!				Average = getData(meshL, Mpos, var,j,i,1)
!			else
!				Average = average / ii
!			end if
!		else if (.not. done) then
!				Average = getData(meshL, Mpos, var,j,i,1)
!				if (.not.added) then
!					do jj = 1, size(Mprev)
!					    if (Mprev(jj)==0) then
!					    	Mprev(jj) = Mpos
!					    	added = .true.
!					    	exit
!					    end if
!					end do
!				end if
!				return
!		end if
!
!	end function vertexAverage

end module Tri2VTK
