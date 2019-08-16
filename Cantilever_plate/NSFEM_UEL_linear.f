C User subroutine for NSFEM method with triangular type of elements
C Mesh is of traingular elements only
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
	SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------    
	INCLUDE 'ABA_PARAM.INC'
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------       
       parameter(zero=0.d0, half=0.5, one=1.d0, two=2.d0, three=3.d0, 
     1 four=4.d0, six=6.d0, eight=8.d0, twelve=12.d0)
C     
     	INTEGER :: numnode,numelem
	REAL(8), DIMENSION(:,:), ALLOCATABLE :: node
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: element
C	
     	open(10,file='D:\PhD_work\Work_June_Dec_2017\ABAQUS_UEL\NSFEM\
     1TRIA\mesh.dat')
      read(10,*) numnode
      if(.not.allocated(node)) allocate(node(numnode,2))
      do i=1,numnode
      read(10,*) (node(i,j),j=1,2)
      enddo
      read(10,*) numelem
      if(.not.allocated(element)) allocate(element(numelem,3))
      do i=1,numelem
      read(10,*) (element(i,j),j=1,3)
      enddo
      close(10)     
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C------- 
C Call the paricular element to perform the analysis
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C------- 
      if(jtype.eq.1) then
C This is for a single node(may be at the corner of global mesh) having 3 surounding nodes 
C (including the node itself) or one adjacent element
         call U1(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
      elseif(jtype.eq.2) then
C This is for a node on the boundries or corners or inside the global mesh surrounded by four
C nodes or having more than one adjacent elements
         call U2(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
	elseif(jtype.eq.3) then
C This is for a node on the boundries or corners or inside the global mesh surrounded by five
C nodes or having more than one adjacent elements
         call U3(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
	elseif(jtype.eq.4) then
C This is for a node on the boundries or corners or inside the global mesh surrounded by six
C nodes or having more than one adjacent elements
         call U4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
	elseif(jtype.eq.5) then
C This is for a node on the boundries or corners or inside the global mesh surrounded by seven
C nodes or having more than one adjacent elements
         call U5(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
	elseif(jtype.eq.6) then
C This is for a node on the boundries or corners or inside the global mesh surrounded by eight
C nodes or having more than one adjacent elements
         call U6(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
	elseif(jtype.eq.7) then
C This is for a node on the boundries or corners or inside the global mesh surrounded by eight
C nodes or having more than one adjacent elements
         call U7(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
	elseif(jtype.eq.8) then
C This is for a node on the boundries or corners or inside the global mesh surrounded by eight
C nodes or having more than one adjacent elements
         call U8(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
	elseif(jtype.eq.9) then
C This is for a node on the boundries or corners or inside the global mesh surrounded by eight
C nodes or having more than one adjacent elements
         call U9(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
	elseif(jtype.eq.10) then
C This is for a node on the boundries or corners or inside the global mesh surrounded by eight
C nodes or having more than one adjacent elements
         call U10(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
	elseif(jtype.eq.11) then
C This is for a node on the boundries or corners or inside the global mesh surrounded by eight
C nodes or having more than one adjacent elements
         call U11(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
      else
C
C         We have a problem...
C
         write(*,*) 'Element type not supported, jtype=',jtype
         write(80,*) 'Element type not supported, jtype=',jtype
         call xit
C
      endif
      return
      end subroutine UEL
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
	SUBROUTINE U1(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------    
	INCLUDE 'ABA_PARAM.INC'
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------       
	DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
	1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
	2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
	3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
	4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
	5     JPROPS(*)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C------- 
C This is for a single node having one adjacent element
C Surrounding nodes are arranged in the same order
C as they are in the original global mesh.
!UEL actually calculates nodal stiffness using the surroudning nodes
!UEL  U1 will have one adjacent element (E1) to a node(N1) under consideration 
!N1-N3: Element nodes
!-------------------------------------------------------------------------------
!                     N3
!                    /  \
!					/    \
!                  /      \
!                 /   E1   \
!                /          \
!               /            \
!              N1-------------N2
!%-------------------------------------------------------------------------------
! A1 One-third area of element E1
! B: Final B-Matrix for node 1
! BT: Matrix contains the transpose matrix of B matrix
! NSTRAIN: Calculated strain of node 1.
! NSTRESS: Calculated stress of node 1.
! ENSNA: Matrix contains strain/stress and area results for each node.
! This matrix is constructed for output purpose.
! THICK: Thickness of the planner.
! E: Young's modulus of the material.
! EMU: Poisson's ratio of the material.	
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Define parameters, dimensions and element type
	parameter(NEG=-1.D0,ZERO=0.D0, HALF=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
	INTEGER :: numnode,numelem
	REAL(8) :: node(numnode,2)
	INTEGER :: element(numelem,3)
	REAL(8) B(3,6), BT(6,3), D(3,3),
	1NSTRAIN(3,1), NSTRESS(3,1), N1(2), N2(2), 
	2N3(2), ENSNA(7,1),NArea,A1,B1(3,6)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
	call get_matzero(AMATRX,6,6)
	call get_matzero(ENSNA,7,1)
	call get_matzero(B,3,6)
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
	N1(1)=COORDS(1,1)
	N1(2)=COORDS(2,1)
	N2(1)=COORDS(1,2)
	N2(2)=COORDS(2,2)
	N3(1)=COORDS(1,3)
	N3(2)=COORDS(2,3)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
	call get_NArea_NBmat(N1,N2,N3,A1,B1)
	NArea=A1/3
	do i = 1,3
		do j = 1,6
			B(i,j) = B(i,j)+(1/NArea)*B1(i,j)
		enddo
	enddo
	call get_mattran(B,BT,3,6)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Assemble Constitutive matrix---D matrix
C---------- Pass material properties into variables
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition
	call get_dmat(E,EMU,D)
C---------- Calculate Kij and assemble K matrix
	call get_kmat(BT,D,B,6,3,3,6,NArea,AMATRX)
C---------- Calculate nodal strain and stress
	do i = 1,3
		do j = 1,6
			NSTRAIN(i,1) = B(i,j)*U(j)
		enddo
	enddo
	do i = 1,3
		do j = 1,3
			NSTRESS(i,1) = D(i,j)*NSTRAIN(j,1)
		enddo
	enddo
C---------- Calculate results for ENSNA
	do i = 1,3
		ENSNA(i,1)=ENSNA(i,1)+NSTRAIN(i,1)
	enddo
	do i = 1,3
	  do j = 1,3
		    ENSNA(i+3,1)=ENSNA(i+3,1)+D(i,j)*ENSNA(j,1)
	  enddo
	enddo
	ENSNA(7,1)=ENSNA(7,1)+NArea
C---------- Calculate the residual force vector
	do k1 = 1, NDOFEL
		do KRHS = 1, NRHS
			RHS(K1,KRHS) = ZERO
		enddo
	enddo
	do I=1,6
		do J=1,6
			RHS(I,1)=RHS(I,1)-AMATRX(I,J)*U(J)*THICK
		enddo
	enddo
C---------- Output the saved results           
	write (*,*) 'Element number'
	write (*,*) JELEM
	write (*,*) 'B matrix'
	write (*,*) B
!	write (*,*) U
!	write (*,*) AMATRX
!	do I=1,3
!		do J=1,7
!			write (*,*) ENSNA(J,I)
!		enddo
!	enddo
	return
	end	subroutine U1      
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
	SUBROUTINE U2(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------    
	INCLUDE 'ABA_PARAM.INC'
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------       
	DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
	1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
	2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
	3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
	4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
	5     JPROPS(*)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C------- 
!This is for a single node having two adjacent elements
!Surrounding nodes are arranged in counter clockwise order such that first node and 
!pair of next two nodess form original element in the global mesh. Next element will
!be formed using first node and pair of nodes formed using last node of earlier formed node
!and the next node in the sequence. The same process will be applicable to a node with more than
!two adjacent elements.
!UEL actually calculates nodal stiffness using the surroudning nodes
!UEL  U2 will have first adjacent element (E1) formed using N1-N2-N3 and second element (E2)
!will be formed using N1-N3-N4
!-------------------------------------------------------------------------------
!                     N3-------------N2
!                    /  \			 /
!					/    \	E1      /
!                  /      \		   /
!                 /   E2   \	  /
!                /          \	 /
!               /            \	/
!              N4-------------N1
!%-------------------------------------------------------------------------------
! C(1) and C(2) are ond-third of areas of elements E1 and E2
! B1,B2: Matrix contains strain-displacement matrix of element E1,E2
! B: Final B-Matrix for node N1
! BT: Matrix contains the transpose matrix of B matrix
! NSTRAIN: Calculated strain of node 1.
! NSTRESS: Calculated stress of node 1.
! ENSNA: Matrix contains strain/stress and area results for each node.
! This matrix is constructed for output purpose.
! THICK: Thickness of the planner.
! E: Young's modulus of the material.
! EMU: Poisson's ratio of the material.
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Define parameters, dimensions and element type
	parameter(NEG=-1.D0,ZERO=0.D0, HALF=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
	INTEGER :: numnode,numelem,Ncen,nn1,nn2,nn3
	REAL(8) :: node(numnode,2)
      INTEGER :: element(numelem,3),Nodnum(4)
	REAL(8) B(3,8), BT(8,3), D(3,3),
	1NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2), 
	2N3(2),ENSNA(7,1),NArea,A1,B1(3,6)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
	call get_matzero(AMATRX,8,8)
	call get_matzero(ENSNA,7,1)
	call get_matzero(B,3,8)
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) node
	do i = 1,4
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=0.d0
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data	
	do ie = 1,numelem
	  A1=0.d0
	  N1(1)=0.d0
	  N1(2)=0.d0
	  N2(1)=0.d0
	  N2(2)=0.d0
	  N3(1)=0.d0
	  N3(2)=0.d0
	  call get_matzero(B1,3,6)
		if (Ncen .EQ. element(ie,1) .OR. Ncen .EQ. element(ie,2)  
	1.OR. Ncen .EQ.  element(ie,3)) then
!	  write (*,*) element(ie,1)
!	  write (*,*) element(ie,2)
!	  write (*,*) element(ie,3)
	  nn1=element(ie,1)
	  nn2=element(ie,2)
	  nn3=element(ie,3)
	  N1(1)=node(nn1,1)
	  N1(2)=node(nn1,2)
	  N2(1)=node(nn2,1)
	  N2(2)=node(nn2,2)
	  N3(1)=node(nn3,1)
	  N3(2)=node(nn3,2)
!	  write (*,*) N1
!	  write (*,*) N3
	  call get_NArea_NBmat(N1,N2,N3,A1,B1)
!	  write (*,*) B1
	  do i=1,3
		    do j = 1,4
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			    B(ii,(2*(j-1)+k)) = B(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
		enddo
		NArea=NArea+A1/3
!		write (*,*) NArea
		else
	  endif
	enddo	
	do i = 1,3
		do j = 1,8
			B(i,j) = (1/NArea)*B(i,j)
		enddo
	enddo
	call get_mattran(B,BT,3,8)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Assemble Constitutive matrix---D matrix
C---------- Pass material properties into variables
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition
	call get_dmat(E,EMU,D)
C---------- Calculate Kij and assemble K matrix
	call get_kmat(BT,D,B,8,3,3,8,NArea,AMATRX)
C---------- Calculate nodal strain and stress
	do i = 1,3
		do j = 1,8
			NSTRAIN(i,1) = B(i,j)*U(j)
		enddo
	enddo
	do i = 1,3
		do j = 1,3
			NSTRESS(i,1) = D(i,j)*NSTRAIN(j,1)
		enddo
	enddo
C---------- Calculate results for ENSNA
	do i = 1,3
		ENSNA(i,1)=ENSNA(i,1)+NSTRAIN(i,1)
	enddo
	do i = 1,3
	  do j = 1,3
		    ENSNA(i+3,1)=ENSNA(i+3,1)+D(i,j)*ENSNA(j,1)
	  enddo
	enddo
	ENSNA(7,1)=ENSNA(7,1)+NArea
C---------- Calculate the residual force vector
	do k1 = 1, NDOFEL
		do KRHS = 1, NRHS
			RHS(K1,KRHS) = ZERO
		enddo
	enddo
	do I=1,8
		do J=1,8
			RHS(I,1)=RHS(I,1)-AMATRX(I,J)*U(J)*THICK
		enddo
	enddo
C---------- Output the saved results 
	write (*,*) 'Element number'
	write (*,*) JELEM
	write (*,*) 'B matrix'
	write (*,*) B
C
!	write (*,*) AMATRX
!	do I=1,3
!		do J=1,7
!			write (*,*) ENSNA(J,I)
!		enddo
!	enddo
	return
	end subroutine U2 
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
	SUBROUTINE U3(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------    
	INCLUDE 'ABA_PARAM.INC'
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------       
	DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
	1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
	2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
	3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
	4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
	5     JPROPS(*)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C------- 
!This is for a single node having two adjacent elements
!Surrounding nodes are arranged in counter clockwise order such that first node and 
!pair of next two nodess form original element in the global mesh. Next element will
!be formed using first node and pair of nodes formed using last node of earlier formed node
!and the next node in the sequence. The same process will be applicable to a node with more than
!two adjacent elements.
!UEL actually calculates nodal stiffness using the surroudning nodes
!UEL  U2 will have first adjacent element (E1) formed using N1-N2-N3 and second element (E2)
!will be formed using N1-N3-N4
!-------------------------------------------------------------------------------
!                     N3-------------N2
!                    /  \			 /
!					/    \	E1      /
!                  /      \		   /
!                 /   E2   \	  /
!                /          \	 /
!               /            \	/
!              N4-------------N1
!%-------------------------------------------------------------------------------
! B1,B2: Matrix contains strain-displacement matrix of element E1
! B: Final B-Matrix for node 1
! BT: Matrix contains the transpose matrix of B matrix
! NSTRAIN: Calculated strain of node 1.
! NSTRESS: Calculated stress of node 1.
! ENSNA: Matrix contains strain/stress and area results for each node.
! This matrix is constructed for output purpose.
! THICK: Thickness of the planner.
! E: Young's modulus of the material.
! EMU: Poisson's ratio of the material.
! Un (n=<10) denotes number of adjacent elements to a node under consideration
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Define parameters, dimensions and element type
	parameter(NEG=-1.D0,ZERO=0.D0, HALF=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
	INTEGER :: numnode,numelem,Ncen,nn1,nn2,nn3
	REAL(8) :: node(numnode,2)
      INTEGER :: element(numelem,3),Nodnum(5)
	REAL(8) B(3,10), BT(10,3), D(3,3),
	1NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2), 
	2N3(2),ENSNA(7,1),NArea,A1,B1(3,6)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
	call get_matzero(AMATRX,10,10)
	call get_matzero(ENSNA,7,1)
	call get_matzero(B,3,10)
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,5
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=0.d0
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data	
	do ie = 1,numelem
	  A1=0.d0
	  N1(1)=0.d0
	  N1(2)=0.d0
	  N2(1)=0.d0
	  N2(2)=0.d0
	  N3(1)=0.d0
	  N3(2)=0.d0
	  call get_matzero(B1,3,6)
		if (Ncen .EQ. element(ie,1) .OR. Ncen .EQ. element(ie,2)  
	1.OR. Ncen .EQ.  element(ie,3)) then
!	  write (*,*) element(ie,1)
!	  write (*,*) element(ie,2)
!	  write (*,*) element(ie,3)
	  nn1=element(ie,1)
	  nn2=element(ie,2)
	  nn3=element(ie,3)
	  N1(1)=node(nn1,1)
	  N1(2)=node(nn1,2)
	  N2(1)=node(nn2,1)
	  N2(2)=node(nn2,2)
	  N3(1)=node(nn3,1)
	  N3(2)=node(nn3,2)
!	  write (*,*) N1
!	  write (*,*) N3
	  call get_NArea_NBmat(N1,N2,N3,A1,B1)
!	  write (*,*) B1
	  do i=1,3
		    do j = 1,5
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			    B(ii,(2*(j-1)+k)) = B(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
		enddo
		NArea=NArea+A1/3
!		write (*,*) NArea
		else
	  endif
	enddo	
	do i = 1,3
		do j = 1,10
			B(i,j) = (1/NArea)*B(i,j)
		enddo
	enddo
	call get_mattran(B,BT,3,10)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Assemble Constitutive matrix---D matrix
C---------- Pass material properties into variables
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition
	call get_dmat(E,EMU,D)
C---------- Calculate Kij and assemble K matrix
	call get_kmat(BT,D,B,10,3,3,10,NArea,AMATRX)
C---------- Calculate nodal strain and stress
	do i = 1,3
		do j = 1,10
			NSTRAIN(i,1) = B(i,j)*U(j)
		enddo
	enddo
	do i = 1,3
		do j = 1,3
			NSTRESS(i,1) = D(i,j)*NSTRAIN(j,1)
		enddo
	enddo
C---------- Calculate results for ENSNA
	do i = 1,3
		ENSNA(i,1)=ENSNA(i,1)+NSTRAIN(i,1)
	enddo
	do i = 1,3
	  do j = 1,3
		    ENSNA(i+3,1)=ENSNA(i+3,1)+D(i,j)*ENSNA(j,1)
	  enddo
	enddo
	ENSNA(7,1)=ENSNA(7,1)+NArea
	
C---------- Calculate the residual force vector
	do k1 = 1, NDOFEL
		do KRHS = 1, NRHS
			RHS(K1,KRHS) = ZERO
		enddo
	enddo
	do I=1,10
		do J=1,10
			RHS(I,1)=RHS(I,1)-AMATRX(I,J)*U(J)*THICK
		enddo
	enddo
C---------- Output the saved results
	write (*,*) 'Element number'
	write (*,*) JELEM
	write (*,*) 'B matrix'
	write (*,*) B
!	write (*,*) NArea
	write (*,*) U
!	write (*,*) AMATRX
!	do I=1,3
!		do J=1,7
!			write (*,*) ENSNA(J,I)
!		enddo
!	enddo

	return
	end	subroutine U3 
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
	SUBROUTINE U4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------    
	INCLUDE 'ABA_PARAM.INC'
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------       
	DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
	1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
	2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
	3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
	4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
	5     JPROPS(*)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C------- 
!This is for a single node having two adjacent elements
!Surrounding nodes are arranged in counter clockwise order such that first node and 
!pair of next two nodess form original element in the global mesh. Next element will
!be formed using first node and pair of nodes formed using last node of earlier formed node
!and the next node in the sequence. The same process will be applicable to a node with more than
!two adjacent elements.
!UEL actually calculates nodal stiffness using the surroudning nodes
!UEL  U2 will have first adjacent element (E1) formed using N1-N2-N3 and second element (E2)
!will be formed using N1-N3-N4
!-------------------------------------------------------------------------------
!                     N3-------------N2
!                    /  \			 /
!					/    \	E1      /
!                  /      \		   /
!                 /   E2   \	  /
!                /          \	 /
!               /            \	/
!              N4-------------N1
!%-------------------------------------------------------------------------------
! B1,B2: Matrix contains strain-displacement matrix of element E1
! B: Final B-Matrix for node 1
! BT: Matrix contains the transpose matrix of B matrix
! NSTRAIN: Calculated strain of node 1.
! NSTRESS: Calculated stress of node 1.
! ENSNA: Matrix contains strain/stress and area results for each node.
! This matrix is constructed for output purpose.
! THICK: Thickness of the planner.
! E: Young's modulus of the material.
! EMU: Poisson's ratio of the material.
! Un (n=<10) denotes number of adjacent elements to a node under consideration
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Define parameters, dimensions and element type
	parameter(NEG=-1.D0,ZERO=0.D0, HALF=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
	INTEGER :: numnode,numelem,Ncen,nn1,nn2,nn3
	REAL(8) :: node(numnode,2)
      INTEGER :: element(numelem,3),Nodnum(6)
	REAL(8) B(3,12), BT(12,3), D(3,3),
	1NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2), 
	2N3(2),ENSNA(7,1),NArea,A1,B1(3,6)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
	call get_matzero(AMATRX,12,12)
	call get_matzero(ENSNA,7,1)
	call get_matzero(B,3,12)
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,6
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=0.d0
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data	
	do ie = 1,numelem
	  A1=0.d0
	  N1(1)=0.d0
	  N1(2)=0.d0
	  N2(1)=0.d0
	  N2(2)=0.d0
	  N3(1)=0.d0
	  N3(2)=0.d0
	  call get_matzero(B1,3,6)
		if (Ncen .EQ. element(ie,1) .OR. Ncen .EQ. element(ie,2)  
	1.OR. Ncen .EQ.  element(ie,3)) then
!	  write (*,*) element(ie,1)
!	  write (*,*) element(ie,2)
!	  write (*,*) element(ie,3)
	  nn1=element(ie,1)
	  nn2=element(ie,2)
	  nn3=element(ie,3)
	  N1(1)=node(nn1,1)
	  N1(2)=node(nn1,2)
	  N2(1)=node(nn2,1)
	  N2(2)=node(nn2,2)
	  N3(1)=node(nn3,1)
	  N3(2)=node(nn3,2)
!	  write (*,*) N1
!	  write (*,*) N3
	  call get_NArea_NBmat(N1,N2,N3,A1,B1)
!	  write (*,*) B1
	  do i=1,3
		    do j = 1,6
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			    B(ii,(2*(j-1)+k)) = B(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
		enddo
		NArea=NArea+A1/3
!		write (*,*) NArea
		else
	  endif
	enddo	
	do i = 1,3
		do j = 1,12
			B(i,j) = (1/NArea)*B(i,j)
		enddo
	enddo
	call get_mattran(B,BT,3,12)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Assemble Constitutive matrix---D matrix
C---------- Pass material properties into variables
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition
	call get_dmat(E,EMU,D)
C---------- Calculate Kij and assemble K matrix
	call get_kmat(BT,D,B,12,3,3,12,NArea,AMATRX)
C---------- Calculate nodal strain and stress
	do i = 1,3
		do j = 1,12
			NSTRAIN(i,1) = B(i,j)*U(j)
		enddo
	enddo
	do i = 1,3
		do j = 1,3
			NSTRESS(i,1) = D(i,j)*NSTRAIN(j,1)
		enddo
	enddo
C---------- Calculate results for ENSNA
	do i = 1,3
		ENSNA(i,1)=ENSNA(i,1)+NSTRAIN(i,1)
	enddo
	do i = 1,3
	  do j = 1,3
		    ENSNA(i+3,1)=ENSNA(i+3,1)+D(i,j)*ENSNA(j,1)
	  enddo
	enddo
	ENSNA(7,1)=ENSNA(7,1)+NArea
	
C---------- Calculate the residual force vector
	do k1 = 1, NDOFEL
		do KRHS = 1, NRHS
			RHS(K1,KRHS) = ZERO
		enddo
	enddo
	do I=1,12
		do J=1,12
			RHS(I,1)=RHS(I,1)-AMATRX(I,J)*U(J)*THICK
		enddo
	enddo
C---------- Output the saved results
	write (*,*) 'Element number'
	write (*,*) JELEM
	write (*,*) 'B matrix'
	write (*,*) B
	write (*,*) U
!	write (*,*) AMATRX
!	do I=1,3
!		do J=1,7
!			write (*,*) ENSNA(J,I)
!		enddo
!	enddo
	return
	end subroutine U4 
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
	SUBROUTINE U5(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------    
	INCLUDE 'ABA_PARAM.INC'
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------       
	DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
	1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
	2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
	3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
	4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
	5     JPROPS(*)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C------- 
!This is for a single node having two adjacent elements
!Surrounding nodes are arranged in counter clockwise order such that first node and 
!pair of next two nodess form original element in the global mesh. Next element will
!be formed using first node and pair of nodes formed using last node of earlier formed node
!and the next node in the sequence. The same process will be applicable to a node with more than
!two adjacent elements.
!UEL actually calculates nodal stiffness using the surroudning nodes
!UEL  U2 will have first adjacent element (E1) formed using N1-N2-N3 and second element (E2)
!will be formed using N1-N3-N4
!-------------------------------------------------------------------------------
!                     N3-------------N2
!                    /  \			 /
!					/    \	E1      /
!                  /      \		   /
!                 /   E2   \	  /
!                /          \	 /
!               /            \	/
!              N4-------------N1
!%-------------------------------------------------------------------------------
! B1,B2: Matrix contains strain-displacement matrix of element E1
! B: Final B-Matrix for node 1
! BT: Matrix contains the transpose matrix of B matrix
! NSTRAIN: Calculated strain of node 1.
! NSTRESS: Calculated stress of node 1.
! ENSNA: Matrix contains strain/stress and area results for each node.
! This matrix is constructed for output purpose.
! THICK: Thickness of the planner.
! E: Young's modulus of the material.
! EMU: Poisson's ratio of the material.
! Un (n=<10) denotes number of adjacent elements to a node under consideration
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Define parameters, dimensions and element type
	parameter(NEG=-1.D0,ZERO=0.D0, HALF=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
	INTEGER :: numnode,numelem,Ncen,nn1,nn2,nn3
	REAL(8) :: node(numnode,2)
      INTEGER :: element(numelem,3),Nodnum(7)
	REAL(8) B(3,14), BT(14,3), D(3,3),
	1NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2), 
	2N3(2),ENSNA(7,1),NArea,A1,B1(3,6)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
	call get_matzero(AMATRX,14,14)
	call get_matzero(ENSNA,7,1)
	call get_matzero(B,3,14)
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,7
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=0.d0
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data	
	do ie = 1,numelem
	  A1=0.d0
	  N1(1)=0.d0
	  N1(2)=0.d0
	  N2(1)=0.d0
	  N2(2)=0.d0
	  N3(1)=0.d0
	  N3(2)=0.d0
	  call get_matzero(B1,3,6)
		if (Ncen .EQ. element(ie,1) .OR. Ncen .EQ. element(ie,2)  
	1.OR. Ncen .EQ.  element(ie,3)) then
!	  write (*,*) element(ie,1)
!	  write (*,*) element(ie,2)
!	  write (*,*) element(ie,3)
	  nn1=element(ie,1)
	  nn2=element(ie,2)
	  nn3=element(ie,3)
	  N1(1)=node(nn1,1)
	  N1(2)=node(nn1,2)
	  N2(1)=node(nn2,1)
	  N2(2)=node(nn2,2)
	  N3(1)=node(nn3,1)
	  N3(2)=node(nn3,2)
!	  write (*,*) N1
!	  write (*,*) N3
	  call get_NArea_NBmat(N1,N2,N3,A1,B1)
!	  write (*,*) B1
	  do i=1,3
		    do j = 1,7
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			    B(ii,(2*(j-1)+k)) = B(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
		enddo
		NArea=NArea+A1/3
!		write (*,*) NArea
		else
	  endif
	enddo	
	do i = 1,3
		do j = 1,14
			B(i,j) = (1/NArea)*B(i,j)
		enddo
	enddo
	call get_mattran(B,BT,3,14)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Assemble Constitutive matrix---D matrix
C---------- Pass material properties into variables
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition
	call get_dmat(E,EMU,D)
C---------- Calculate Kij and assemble K matrix
	call get_kmat(BT,D,B,14,3,3,14,NArea,AMATRX)
C---------- Calculate nodal strain and stress
	do i = 1,3
		do j = 1,14
			NSTRAIN(i,1) = B(i,j)*U(j)
		enddo
	enddo
	do i = 1,3
		do j = 1,3
			NSTRESS(i,1) = D(i,j)*NSTRAIN(j,1)
		enddo
	enddo
C---------- Calculate results for ENSNA
	do i = 1,3
		ENSNA(i,1)=ENSNA(i,1)+NSTRAIN(i,1)
	enddo
	do i = 1,3
	  do j = 1,3
		    ENSNA(i+3,1)=ENSNA(i+3,1)+D(i,j)*ENSNA(j,1)
	  enddo
	enddo
	ENSNA(7,1)=ENSNA(7,1)+NArea
	
C---------- Calculate the residual force vector
	do k1 = 1, NDOFEL
		do KRHS = 1, NRHS
			RHS(K1,KRHS) = ZERO
		enddo
	enddo
	do I=1,14
		do J=1,14
			RHS(I,1)=RHS(I,1)-AMATRX(I,J)*U(J)*THICK
		enddo
	enddo
C---------- Output the saved results
	write (*,*) 'Element number'
	write (*,*) JELEM
	write (*,*) 'B matrix'
	write (*,*) B
!	write (*,*) NArea
	write (*,*) U
!	write (*,*) AMATRX
!	do I=1,3
!		do J=1,7
!			write (*,*) ENSNA(J,I)
!		enddo
!	enddo
	return
	end subroutine U5 
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
	SUBROUTINE U6(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------    
	INCLUDE 'ABA_PARAM.INC'
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------       
	DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
	1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
	2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
	3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
	4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
	5     JPROPS(*)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C------- 
!This is for a single node having two adjacent elements
!Surrounding nodes are arranged in counter clockwise order such that first node and 
!pair of next two nodess form original element in the global mesh. Next element will
!be formed using first node and pair of nodes formed using last node of earlier formed node
!and the next node in the sequence. The same process will be applicable to a node with more than
!two adjacent elements.
!UEL actually calculates nodal stiffness using the surroudning nodes
!UEL  U2 will have first adjacent element (E1) formed using N1-N2-N3 and second element (E2)
!will be formed using N1-N3-N4
!-------------------------------------------------------------------------------
!                     N3-------------N2
!                    /  \			 /
!					/    \	E1      /
!                  /      \		   /
!                 /   E2   \	  /
!                /          \	 /
!               /            \	/
!              N4-------------N1
!%-------------------------------------------------------------------------------
! B1,B2: Matrix contains strain-displacement matrix of element E1
! B: Final B-Matrix for node 1
! BT: Matrix contains the transpose matrix of B matrix
! NSTRAIN: Calculated strain of node 1.
! NSTRESS: Calculated stress of node 1.
! ENSNA: Matrix contains strain/stress and area results for each node.
! This matrix is constructed for output purpose.
! THICK: Thickness of the planner.
! E: Young's modulus of the material.
! EMU: Poisson's ratio of the material.
! Un (n=<10) denotes number of adjacent elements to a node under consideration
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Define parameters, dimensions and element type
	parameter(NEG=-1.D0,ZERO=0.D0, HALF=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
	INTEGER :: numnode,numelem,Ncen,nn1,nn2,nn3
	REAL(8) :: node(numnode,2)
      INTEGER :: element(numelem,3),Nodnum(8)
	REAL(8) B(3,16), BT(16,3), D(3,3),
	1NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2), 
	2N3(2),ENSNA(7,1),NArea,A1,B1(3,6)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
	call get_matzero(AMATRX,16,16)
	call get_matzero(ENSNA,7,1)
	call get_matzero(B,3,16)
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,8
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=0.d0
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data	
	do ie = 1,numelem
	  A1=0.d0
	  N1(1)=0.d0
	  N1(2)=0.d0
	  N2(1)=0.d0
	  N2(2)=0.d0
	  N3(1)=0.d0
	  N3(2)=0.d0
	  call get_matzero(B1,3,6)
		if (Ncen .EQ. element(ie,1) .OR. Ncen .EQ. element(ie,2)  
	1.OR. Ncen .EQ.  element(ie,3)) then
	  nn1=element(ie,1)
	  nn2=element(ie,2)
	  nn3=element(ie,3)
	  N1(1)=node(nn1,1)
	  N1(2)=node(nn1,2)
	  N2(1)=node(nn2,1)
	  N2(2)=node(nn2,2)
	  N3(1)=node(nn3,1)
	  N3(2)=node(nn3,2)
!	  write (*,*) N1
!	  write (*,*) N3
	  call get_NArea_NBmat(N1,N2,N3,A1,B1)
!	  write (*,*) B1
	  do i=1,3
		    do j = 1,8
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			    B(ii,(2*(j-1)+k)) = B(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
		enddo
		NArea=NArea+A1/3
!		write (*,*) NArea
		else
	  endif
	enddo	
	do i = 1,3
		do j = 1,16
			B(i,j) = (1/NArea)*B(i,j)
		enddo
	enddo
	call get_mattran(B,BT,3,16)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Assemble Constitutive matrix---D matrix
C---------- Pass material properties into variables
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition
	call get_dmat(E,EMU,D)
C---------- Calculate Kij and assemble K matrix
	call get_kmat(BT,D,B,16,3,3,16,NArea,AMATRX)
C---------- Calculate nodal strain and stress
	do i = 1,3
		do j = 1,16
			NSTRAIN(i,1) = B(i,j)*U(j)
		enddo
	enddo
	do i = 1,3
		do j = 1,3
			NSTRESS(i,1) = D(i,j)*NSTRAIN(j,1)
		enddo
	enddo
C---------- Calculate results for ENSNA
	do i = 1,3
		ENSNA(i,1)=ENSNA(i,1)+NSTRAIN(i,1)
	enddo
	do i = 1,3
	  do j = 1,3
		    ENSNA(i+3,1)=ENSNA(i+3,1)+D(i,j)*ENSNA(j,1)
	  enddo
	enddo
	ENSNA(7,1)=ENSNA(7,1)+NArea
	
C---------- Calculate the residual force vector
	do k1 = 1, NDOFEL
		do KRHS = 1, NRHS
			RHS(K1,KRHS) = ZERO
		enddo
	enddo
	do I=1,16
		do J=1,16
			RHS(I,1)=RHS(I,1)-AMATRX(I,J)*U(J)*THICK
		enddo
	enddo
C---------- Output the saved results
!	write (*,*) element
!	write (*,*) node
	write (*,*) 'Element number'
	write (*,*) JELEM
	write (*,*) 'B matrix'
	write (*,*) B
!	write (*,*) NArea
	write (*,*) U
!	write (*,*) AMATRX
!	do I=1,3
!		do J=1,7
!			write (*,*) ENSNA(J,I)
!		enddo
!	enddo
	return
	end subroutine U6 
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------	
	SUBROUTINE U7(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------    
	INCLUDE 'ABA_PARAM.INC'
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------       
	DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
	1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
	2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
	3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
	4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
	5     JPROPS(*)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C------- 
!This is for a single node having two adjacent elements
!Surrounding nodes are arranged in counter clockwise order such that first node and 
!pair of next two nodess form original element in the global mesh. Next element will
!be formed using first node and pair of nodes formed using last node of earlier formed node
!and the next node in the sequence. The same process will be applicable to a node with more than
!two adjacent elements.
!UEL actually calculates nodal stiffness using the surroudning nodes
!UEL  U2 will have first adjacent element (E1) formed using N1-N2-N3 and second element (E2)
!will be formed using N1-N3-N4
!-------------------------------------------------------------------------------
!                     N3-------------N2
!                    /  \			 /
!					/    \	E1      /
!                  /      \		   /
!                 /   E2   \	  /
!                /          \	 /
!               /            \	/
!              N4-------------N1
!%-------------------------------------------------------------------------------
! B1,B2: Matrix contains strain-displacement matrix of element E1
! B: Final B-Matrix for node 1
! BT: Matrix contains the transpose matrix of B matrix
! NSTRAIN: Calculated strain of node 1.
! NSTRESS: Calculated stress of node 1.
! ENSNA: Matrix contains strain/stress and area results for each node.
! This matrix is constructed for output purpose.
! THICK: Thickness of the planner.
! E: Young's modulus of the material.
! EMU: Poisson's ratio of the material.
! Un (n=<10) denotes number of adjacent elements to a node under consideration
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Define parameters, dimensions and element type
	parameter(NEG=-1.D0,ZERO=0.D0, HALF=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
	INTEGER :: numnode,numelem,Ncen,nn1,nn2,nn3
	REAL(8) :: node(numnode,2)
      INTEGER :: element(numelem,3),Nodnum(9)
	REAL(8) B(3,18), BT(18,3), D(3,3),
	1NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2), 
	2N3(2),ENSNA(7,1),NArea,A1,B1(3,6)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
	call get_matzero(AMATRX,18,18)
	call get_matzero(ENSNA,7,1)
	call get_matzero(B,3,18)
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,9
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=0.d0
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data	
	do ie = 1,numelem
	  A1=0.d0
	  N1(1)=0.d0
	  N1(2)=0.d0
	  N2(1)=0.d0
	  N2(2)=0.d0
	  N3(1)=0.d0
	  N3(2)=0.d0
	  call get_matzero(B1,3,6)
		if (Ncen .EQ. element(ie,1) .OR. Ncen .EQ. element(ie,2)  
	1.OR. Ncen .EQ.  element(ie,3)) then
	  nn1=element(ie,1)
	  nn2=element(ie,2)
	  nn3=element(ie,3)
	  N1(1)=node(nn1,1)
	  N1(2)=node(nn1,2)
	  N2(1)=node(nn2,1)
	  N2(2)=node(nn2,2)
	  N3(1)=node(nn3,1)
	  N3(2)=node(nn3,2)
!	  write (*,*) N1
!	  write (*,*) N3
	  call get_NArea_NBmat(N1,N2,N3,A1,B1)
!	  write (*,*) B1
	  do i=1,3
		    do j = 1,9
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			    B(ii,(2*(j-1)+k)) = B(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
		enddo
		NArea=NArea+A1/3
!		write (*,*) NArea
		else
	  endif
	enddo	
	do i = 1,3
		do j = 1,18
			B(i,j) = (1/NArea)*B(i,j)
		enddo
	enddo
	call get_mattran(B,BT,3,18)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Assemble Constitutive matrix---D matrix
C---------- Pass material properties into variables
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition
	call get_dmat(E,EMU,D)
C---------- Calculate Kij and assemble K matrix
	call get_kmat(BT,D,B,18,3,3,18,NArea,AMATRX)
C---------- Calculate nodal strain and stress
	do i = 1,3
		do j = 1,18
			NSTRAIN(i,1) = B(i,j)*U(j)
		enddo
	enddo
	do i = 1,3
		do j = 1,3
			NSTRESS(i,1) = D(i,j)*NSTRAIN(j,1)
		enddo
	enddo
C---------- Calculate results for ENSNA
	do i = 1,3
		ENSNA(i,1)=ENSNA(i,1)+NSTRAIN(i,1)
	enddo
	do i = 1,3
	  do j = 1,3
		    ENSNA(i+3,1)=ENSNA(i+3,1)+D(i,j)*ENSNA(j,1)
	  enddo
	enddo
	ENSNA(7,1)=ENSNA(7,1)+NArea
	
C---------- Calculate the residual force vector
	do k1 = 1, NDOFEL
		do KRHS = 1, NRHS
			RHS(K1,KRHS) = ZERO
		enddo
	enddo
	do I=1,18
		do J=1,18
			RHS(I,1)=RHS(I,1)-AMATRX(I,J)*U(J)*THICK
		enddo
	enddo
C---------- Output the saved results
	write (*,*) JELEM
	write (*,*) B
!	write (*,*) NArea
	write (*,*) U
!	write (*,*) AMATRX
!	do I=1,3
!		do J=1,7
!			write (*,*) ENSNA(J,I)
!		enddo
!	enddo
	return
	end subroutine U7 
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
	SUBROUTINE U8(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------    
	INCLUDE 'ABA_PARAM.INC'
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------       
	DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
	1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
	2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
	3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
	4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
	5     JPROPS(*)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C------- 
!This is for a single node having two adjacent elements
!Surrounding nodes are arranged in counter clockwise order such that first node and 
!pair of next two nodess form original element in the global mesh. Next element will
!be formed using first node and pair of nodes formed using last node of earlier formed node
!and the next node in the sequence. The same process will be applicable to a node with more than
!two adjacent elements.
!UEL actually calculates nodal stiffness using the surroudning nodes
!UEL  U2 will have first adjacent element (E1) formed using N1-N2-N3 and second element (E2)
!will be formed using N1-N3-N4
!-------------------------------------------------------------------------------
!                     N3-------------N2
!                    /  \			 /
!					/    \	E1      /
!                  /      \		   /
!                 /   E2   \	  /
!                /          \	 /
!               /            \	/
!              N4-------------N1
!%-------------------------------------------------------------------------------
! B1,B2: Matrix contains strain-displacement matrix of element E1
! B: Final B-Matrix for node 1
! BT: Matrix contains the transpose matrix of B matrix
! NSTRAIN: Calculated strain of node 1.
! NSTRESS: Calculated stress of node 1.
! ENSNA: Matrix contains strain/stress and area results for each node.
! This matrix is constructed for output purpose.
! THICK: Thickness of the planner.
! E: Young's modulus of the material.
! EMU: Poisson's ratio of the material.
! Un (n=<10) denotes number of adjacent elements to a node under consideration
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Define parameters, dimensions and element type
	parameter(NEG=-1.D0,ZERO=0.D0, HALF=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
	INTEGER :: numnode,numelem,Ncen,nn1,nn2,nn3
	REAL(8) :: node(numnode,2)
      INTEGER :: element(numelem,3),Nodnum(10)
	REAL(8) B(3,20), BT(20,3), D(3,3),
	1NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2), 
	2N3(2),ENSNA(7,1),NArea,A1,B1(3,6)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
	call get_matzero(AMATRX,20,20)
	call get_matzero(ENSNA,7,1)
	call get_matzero(B,3,20)
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,10
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=0.d0
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data	
	do ie = 1,numelem
	  A1=0.d0
	  N1(1)=0.d0
	  N1(2)=0.d0
	  N2(1)=0.d0
	  N2(2)=0.d0
	  N3(1)=0.d0
	  N3(2)=0.d0
	  call get_matzero(B1,3,6)
		if (Ncen .EQ. element(ie,1) .OR. Ncen .EQ. element(ie,2)  
	1.OR. Ncen .EQ.  element(ie,3)) then
	  nn1=element(ie,1)
	  nn2=element(ie,2)
	  nn3=element(ie,3)
	  N1(1)=node(nn1,1)
	  N1(2)=node(nn1,2)
	  N2(1)=node(nn2,1)
	  N2(2)=node(nn2,2)
	  N3(1)=node(nn3,1)
	  N3(2)=node(nn3,2)
!	  write (*,*) N1
!	  write (*,*) N3
	  call get_NArea_NBmat(N1,N2,N3,A1,B1)
!	  write (*,*) B1
	  do i=1,3
		    do j = 1,10
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			    B(ii,(2*(j-1)+k)) = B(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
		enddo
		NArea=NArea+A1/3
!		write (*,*) NArea
		else
	  endif
	enddo	
	do i = 1,3
		do j = 1,20
			B(i,j) = (1/NArea)*B(i,j)
		enddo
	enddo
	call get_mattran(B,BT,3,20)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Assemble Constitutive matrix---D matrix
C---------- Pass material properties into variables
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition
	call get_dmat(E,EMU,D)
C---------- Calculate Kij and assemble K matrix
	call get_kmat(BT,D,B,20,3,3,20,NArea,AMATRX)
C---------- Calculate nodal strain and stress
	do i = 1,3
		do j = 1,20
			NSTRAIN(i,1) = B(i,j)*U(j)
		enddo
	enddo
	do i = 1,3
		do j = 1,3
			NSTRESS(i,1) = D(i,j)*NSTRAIN(j,1)
		enddo
	enddo
C---------- Calculate results for ENSNA
	do i = 1,3
		ENSNA(i,1)=ENSNA(i,1)+NSTRAIN(i,1)
	enddo
	do i = 1,3
	  do j = 1,3
		    ENSNA(i+3,1)=ENSNA(i+3,1)+D(i,j)*ENSNA(j,1)
	  enddo
	enddo
	ENSNA(7,1)=ENSNA(7,1)+NArea
C---------- Calculate the residual force vector
	do k1 = 1, NDOFEL
		do KRHS = 1, NRHS
			RHS(K1,KRHS) = ZERO
		enddo
	enddo
	do I=1,20
		do J=1,20
			RHS(I,1)=RHS(I,1)-AMATRX(I,J)*U(J)*THICK
		enddo
	enddo
C---------- Output the saved results
	write (*,*) JELEM
	write (*,*) B
!	write (*,*) NArea
	write (*,*) U
!	write (*,*) AMATRX
!	do I=1,3
!		do J=1,7
!			write (*,*) ENSNA(J,I)
!		enddo
!	enddo
	return
	end subroutine U8
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
	SUBROUTINE U9(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------    
	INCLUDE 'ABA_PARAM.INC'
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------       
	DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
	1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
	2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
	3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
	4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
	5     JPROPS(*)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C------- 
!This is for a single node having two adjacent elements
!Surrounding nodes are arranged in counter clockwise order such that first node and 
!pair of next two nodess form original element in the global mesh. Next element will
!be formed using first node and pair of nodes formed using last node of earlier formed node
!and the next node in the sequence. The same process will be applicable to a node with more than
!two adjacent elements.
!UEL actually calculates nodal stiffness using the surroudning nodes
!UEL  U2 will have first adjacent element (E1) formed using N1-N2-N3 and second element (E2)
!will be formed using N1-N3-N4
!-------------------------------------------------------------------------------
!                     N3-------------N2
!                    /  \			 /
!					/    \	E1      /
!                  /      \		   /
!                 /   E2   \	  /
!                /          \	 /
!               /            \	/
!              N4-------------N1
!%-------------------------------------------------------------------------------
! B1,B2: Matrix contains strain-displacement matrix of element E1
! B: Final B-Matrix for node 1
! BT: Matrix contains the transpose matrix of B matrix
! NSTRAIN: Calculated strain of node 1.
! NSTRESS: Calculated stress of node 1.
! ENSNA: Matrix contains strain/stress and area results for each node.
! This matrix is constructed for output purpose.
! THICK: Thickness of the planner.
! E: Young's modulus of the material.
! EMU: Poisson's ratio of the material.
! Un (n=<10) denotes number of adjacent elements to a node under consideration
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Define parameters, dimensions and element type
	parameter(NEG=-1.D0,ZERO=0.D0, HALF=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
	INTEGER :: numnode,numelem,Ncen,nn1,nn2,nn3
	REAL(8) :: node(numnode,2)
      INTEGER :: element(numelem,3),Nodnum(11)
	REAL(8) B(3,22), BT(22,3), D(3,3),
	1NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2), 
	2N3(2),ENSNA(7,1),NArea,A1,B1(3,6)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
	call get_matzero(AMATRX,22,22)
	call get_matzero(ENSNA,7,1)
	call get_matzero(B,3,22)
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,11
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=0.d0
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data	
	do ie = 1,numelem
	  A1=0.d0
	  N1(1)=0.d0
	  N1(2)=0.d0
	  N2(1)=0.d0
	  N2(2)=0.d0
	  N3(1)=0.d0
	  N3(2)=0.d0
	  call get_matzero(B1,3,6)
		if (Ncen .EQ. element(ie,1) .OR. Ncen .EQ. element(ie,2)  
	1.OR. Ncen .EQ.  element(ie,3)) then
	  nn1=element(ie,1)
	  nn2=element(ie,2)
	  nn3=element(ie,3)
	  N1(1)=node(nn1,1)
	  N1(2)=node(nn1,2)
	  N2(1)=node(nn2,1)
	  N2(2)=node(nn2,2)
	  N3(1)=node(nn3,1)
	  N3(2)=node(nn3,2)
!	  write (*,*) N1
!	  write (*,*) N3
	  call get_NArea_NBmat(N1,N2,N3,A1,B1)
!	  write (*,*) B1
	  do i=1,3
		    do j = 1,11
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			    B(ii,(2*(j-1)+k)) = B(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
		enddo
		NArea=NArea+A1/3
!		write (*,*) NArea
		else
	  endif
	enddo	
	do i = 1,3
		do j = 1,22
			B(i,j) = (1/NArea)*B(i,j)
		enddo
	enddo
	call get_mattran(B,BT,3,22)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Assemble Constitutive matrix---D matrix
C---------- Pass material properties into variables
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition
	call get_dmat(E,EMU,D)
C---------- Calculate Kij and assemble K matrix
	call get_kmat(BT,D,B,22,3,3,22,NArea,AMATRX)
C---------- Calculate nodal strain and stress
	do i = 1,3
		do j = 1,22
			NSTRAIN(i,1) = B(i,j)*U(j)
		enddo
	enddo
	do i = 1,3
		do j = 1,3
			NSTRESS(i,1) = D(i,j)*NSTRAIN(j,1)
		enddo
	enddo
C---------- Calculate results for ENSNA
	do i = 1,3
		ENSNA(i,1)=ENSNA(i,1)+NSTRAIN(i,1)
	enddo
	do i = 1,3
	  do j = 1,3
		    ENSNA(i+3,1)=ENSNA(i+3,1)+D(i,j)*ENSNA(j,1)
	  enddo
	enddo
	ENSNA(7,1)=ENSNA(7,1)+NArea
C---------- Calculate the residual force vector
	do k1 = 1, NDOFEL
		do KRHS = 1, NRHS
			RHS(K1,KRHS) = ZERO
		enddo
	enddo
	do I=1,22
		do J=1,22
			RHS(I,1)=RHS(I,1)-AMATRX(I,J)*U(J)*THICK
		enddo
	enddo
C---------- Output the saved results
	write (*,*) JELEM
	write (*,*) B
!	write (*,*) NArea
	write (*,*) U
!	write (*,*) AMATRX
!	do I=1,3
!		do J=1,7
!			write (*,*) ENSNA(J,I)
!		enddo
!	enddo
	return
	end subroutine U9 
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
	SUBROUTINE U10(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------    
	INCLUDE 'ABA_PARAM.INC'
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------       
	DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
	1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
	2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
	3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
	4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
	5     JPROPS(*)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C------- 
!This is for a single node having two adjacent elements
!Surrounding nodes are arranged in counter clockwise order such that first node and 
!pair of next two nodess form original element in the global mesh. Next element will
!be formed using first node and pair of nodes formed using last node of earlier formed node
!and the next node in the sequence. The same process will be applicable to a node with more than
!two adjacent elements.
!UEL actually calculates nodal stiffness using the surroudning nodes
!UEL  U2 will have first adjacent element (E1) formed using N1-N2-N3 and second element (E2)
!will be formed using N1-N3-N4
!-------------------------------------------------------------------------------
!                     N3-------------N2
!                    /  \			 /
!					/    \	E1      /
!                  /      \		   /
!                 /   E2   \	  /
!                /          \	 /
!               /            \	/
!              N4-------------N1
!%-------------------------------------------------------------------------------
! B1,B2: Matrix contains strain-displacement matrix of element E1
! B: Final B-Matrix for node 1
! BT: Matrix contains the transpose matrix of B matrix
! NSTRAIN: Calculated strain of node 1.
! NSTRESS: Calculated stress of node 1.
! ENSNA: Matrix contains strain/stress and area results for each node.
! This matrix is constructed for output purpose.
! THICK: Thickness of the planner.
! E: Young's modulus of the material.
! EMU: Poisson's ratio of the material.
! Un (n=<10) denotes number of adjacent elements to a node under consideration
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Define parameters, dimensions and element type
	parameter(NEG=-1.D0,ZERO=0.D0, HALF=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
	INTEGER :: numnode,numelem,Ncen,nn1,nn2,nn3
	REAL(8) :: node(numnode,2)
      INTEGER :: element(numelem,3),Nodnum(12)
	REAL(8) B(3,24), BT(24,3), D(3,3),
	1NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2), 
	2N3(2),ENSNA(7,1),NArea,A1,B1(3,6)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
	call get_matzero(AMATRX,24,24)
	call get_matzero(ENSNA,7,1)
	call get_matzero(B,3,24)
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,12
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=0.d0
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data	
	do ie = 1,numelem
	  A1=0.d0
	  N1(1)=0.d0
	  N1(2)=0.d0
	  N2(1)=0.d0
	  N2(2)=0.d0
	  N3(1)=0.d0
	  N3(2)=0.d0
	  call get_matzero(B1,3,6)
		if (Ncen .EQ. element(ie,1) .OR. Ncen .EQ. element(ie,2)  
	1.OR. Ncen .EQ.  element(ie,3)) then
	  nn1=element(ie,1)
	  nn2=element(ie,2)
	  nn3=element(ie,3)
	  N1(1)=node(nn1,1)
	  N1(2)=node(nn1,2)
	  N2(1)=node(nn2,1)
	  N2(2)=node(nn2,2)
	  N3(1)=node(nn3,1)
	  N3(2)=node(nn3,2)
!	  write (*,*) N1
!	  write (*,*) N3
	  call get_NArea_NBmat(N1,N2,N3,A1,B1)
!	  write (*,*) B1
	  do i=1,3
		    do j = 1,12
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			    B(ii,(2*(j-1)+k)) = B(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
		enddo
		NArea=NArea+A1/3
!		write (*,*) NArea
		else
	  endif
	enddo	
	do i = 1,3
		do j = 1,24
			B(i,j) = (1/NArea)*B(i,j)
		enddo
	enddo
	call get_mattran(B,BT,3,24)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Assemble Constitutive matrix---D matrix
C---------- Pass material properties into variables
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition
	call get_dmat(E,EMU,D)
C---------- Calculate Kij and assemble K matrix
	call get_kmat(BT,D,B,24,3,3,24,NArea,AMATRX)
C---------- Calculate nodal strain and stress
	do i = 1,3
		do j = 1,24
			NSTRAIN(i,1) = B(i,j)*U(j)
		enddo
	enddo
	do i = 1,3
		do j = 1,3
			NSTRESS(i,1) = D(i,j)*NSTRAIN(j,1)
		enddo
	enddo
C---------- Calculate results for ENSNA
	do i = 1,3
		ENSNA(i,1)=ENSNA(i,1)+NSTRAIN(i,1)
	enddo
	do i = 1,3
	  do j = 1,3
		    ENSNA(i+3,1)=ENSNA(i+3,1)+D(i,j)*ENSNA(j,1)
	  enddo
	enddo
	ENSNA(7,1)=ENSNA(7,1)+NArea
C---------- Calculate the residual force vector
	do k1 = 1, NDOFEL
		do KRHS = 1, NRHS
			RHS(K1,KRHS) = ZERO
		enddo
	enddo
	do I=1,24
		do J=1,24
			RHS(I,1)=RHS(I,1)-AMATRX(I,J)*U(J)*THICK
		enddo
	enddo
C---------- Output the saved results
	write (*,*) JELEM
	write (*,*) B
!	write (*,*) NArea
	write (*,*) U
!	write (*,*) AMATRX
!	do I=1,3
!		do J=1,7
!			write (*,*) ENSNA(J,I)
!		enddo
!	enddo
	return
	end subroutine U10 
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
	SUBROUTINE U11(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------    
	INCLUDE 'ABA_PARAM.INC'
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------       
	DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
	1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
	2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
	3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
	4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
	5     JPROPS(*)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C------- 
!This is for a single node having two adjacent elements
!Surrounding nodes are arranged in counter clockwise order such that first node and 
!pair of next two nodess form original element in the global mesh. Next element will
!be formed using first node and pair of nodes formed using last node of earlier formed node
!and the next node in the sequence. The same process will be applicable to a node with more than
!two adjacent elements.
!UEL actually calculates nodal stiffness using the surroudning nodes
!UEL  U2 will have first adjacent element (E1) formed using N1-N2-N3 and second element (E2)
!will be formed using N1-N3-N4
!-------------------------------------------------------------------------------
!                     N3-------------N2
!                    /  \			 /
!					/    \	E1      /
!                  /      \		   /
!                 /   E2   \	  /
!                /          \	 /
!               /            \	/
!              N4-------------N1
!%-------------------------------------------------------------------------------
! B1,B2: Matrix contains strain-displacement matrix of element E1
! B: Final B-Matrix for node 1
! BT: Matrix contains the transpose matrix of B matrix
! NSTRAIN: Calculated strain of node 1.
! NSTRESS: Calculated stress of node 1.
! ENSNA: Matrix contains strain/stress and area results for each node.
! This matrix is constructed for output purpose.
! THICK: Thickness of the planner.
! E: Young's modulus of the material.
! EMU: Poisson's ratio of the material.
! Un (n=<10) denotes number of adjacent elements to a node under consideration
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Define parameters, dimensions and element type
	parameter(NEG=-1.D0,ZERO=0.D0, HALF=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
	INTEGER :: numnode,numelem,Ncen,nn1,nn2,nn3
	REAL(8) :: node(numnode,2)
      INTEGER :: element(numelem,3),Nodnum(13)
	REAL(8) B(3,26), BT(26,3), D(3,3),
	1NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2), 
	2N3(2),ENSNA(7,1),NArea,A1,B1(3,6)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
	call get_matzero(AMATRX,26,26)
	call get_matzero(ENSNA,7,1)
	call get_matzero(B,3,26)
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,13
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=0.d0
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data	
	do ie = 1,numelem
	  A1=0.d0
	  N1(1)=0.d0
	  N1(2)=0.d0
	  N2(1)=0.d0
	  N2(2)=0.d0
	  N3(1)=0.d0
	  N3(2)=0.d0
	  call get_matzero(B1,3,6)
		if (Ncen .EQ. element(ie,1) .OR. Ncen .EQ. element(ie,2)  
	1.OR. Ncen .EQ.  element(ie,3)) then
	  nn1=element(ie,1)
	  nn2=element(ie,2)
	  nn3=element(ie,3)
	  N1(1)=node(nn1,1)
	  N1(2)=node(nn1,2)
	  N2(1)=node(nn2,1)
	  N2(2)=node(nn2,2)
	  N3(1)=node(nn3,1)
	  N3(2)=node(nn3,2)
!	  write (*,*) N1
!	  write (*,*) N3
	  call get_NArea_NBmat(N1,N2,N3,A1,B1)
!	  write (*,*) B1
	  do i=1,3
		    do j = 1,13
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			    B(ii,(2*(j-1)+k)) = B(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
		enddo
		NArea=NArea+A1/3
!		write (*,*) NArea
		else
	  endif
	enddo	
	do i = 1,3
		do j = 1,26
			B(i,j) = (1/NArea)*B(i,j)
		enddo
	enddo
	call get_mattran(B,BT,3,26)
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Assemble Constitutive matrix---D matrix
C---------- Pass material properties into variables
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition
	call get_dmat(E,EMU,D)
C---------- Calculate Kij and assemble K matrix
	call get_kmat(BT,D,B,26,3,3,26,NArea,AMATRX)
C---------- Calculate nodal strain and stress
	do i = 1,3
		do j = 1,26
			NSTRAIN(i,1) = B(i,j)*U(j)
		enddo
	enddo
	do i = 1,3
		do j = 1,3
			NSTRESS(i,1) = D(i,j)*NSTRAIN(j,1)
		enddo
	enddo
C---------- Calculate results for ENSNA
	do i = 1,3
		ENSNA(i,1)=ENSNA(i,1)+NSTRAIN(i,1)
	enddo
	do i = 1,3
	  do j = 1,3
		    ENSNA(i+3,1)=ENSNA(i+3,1)+D(i,j)*ENSNA(j,1)
	  enddo
	enddo
	ENSNA(7,1)=ENSNA(7,1)+NArea
C---------- Calculate the residual force vector
	do k1 = 1, NDOFEL
		do KRHS = 1, NRHS
			RHS(K1,KRHS) = ZERO
		enddo
	enddo
	do I=1,26
		do J=1,26
			RHS(I,1)=RHS(I,1)-AMATRX(I,J)*U(J)*THICK
		enddo
	enddo
C---------- Output the saved results
	write (*,*) JELEM
	write (*,*) B
!	write (*,*) NArea
	write (*,*) U
!	write (*,*) AMATRX
!	do I=1,3
!		do J=1,7
!			write (*,*) ENSNA(J,I)
!		enddo
!	enddo
	return
	end subroutine U11 
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------  
C SUBROUTINES C
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C-Subroutine to calculate Narea and B matrix for node based smoothening domain
C
	subroutine get_NArea_NBmat(p1,p2,p3,area,b)
	INCLUDE 'ABA_PARAM.INC'
	real(8) x(3),y(3),
	1p1(2),p2(2),p3(2),area,tt(6),
	2b(3,6),A(3,3),DetA,A1
C---------------------------------------------------------------------------C
C Variables
C p1-p3: Coordinates of the three nodes forming the original traingular element.
C x: X coordinates of the three nodes forming the original traingular element.
C y: Y coordinates of the three nodes forming the original traingular element.
C area: Area of the node based smoothening domain
C b: B matrix the node based smoothening domain
C---------------------------------------------------------------------------C
!	write (*,*) p1
!	write (*,*) p2
!	write (*,*) p3
	x(1)=p1(1)
	y(1)=p1(2)
	x(2)=p2(1)
	y(2)=p2(2)
	x(3)=p3(1)
	y(3)=p3(2)
C	
	do i = 1, 3
		A(1,i)=1
		A(2,i)=x(i)
		A(3,i)=y(i)
	enddo
	DetA = A(2,2)*A(3,3)-A(3,2)*A(2,3)-A(2,1)*A(3,3)+A(2,3)*A(3,1)+
	1A(2,1)*A(3,2)-A(2,2)*A(3,1)
	A1 = abs(DetA)
	area = A1/2
C
	tt(1)=x(3)-x(2)
	tt(2)=y(2)-y(3)
	tt(3)=x(1)-x(3)
	tt(4)=y(3)-y(1)
	tt(5)=x(2)-x(1)
	tt(6)=y(1)-y(2)
C
	call get_matzero(b,3,6)
C
	do i = 1, 3
	  b(1,2*i-1)=b(1,2*i-1)+tt(2*i)/6
	  b(2,2*i)=b(2,2*i)+tt(2*i-1)/6
	enddo
	do i = 1, 6
		b(3,i)=b(3,i)+tt(i)/6
	enddo
!	write (*,*) tt(1)
!	write (*,*) tt(3)
!	write (*,*) tt(5)
!	write (*,*) b
C	
	return
	end
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C############### Subroutine to calculate material matrix
	subroutine get_dmat(e,mu,d)
	INCLUDE 'ABA_PARAM.INC'
	real(8) e,mu,d(3,3)
C---------- Plane Stress condition
	d(1,1)=e/(1.0-mu*mu)
	d(1,2)=(e*mu)/(1.0-mu*mu)
	d(1,3)=0.0
	d(2,1)=(e*mu)/(1.0-mu*mu)
	d(2,2)=e/(1.0-mu*mu)
	d(2,3)=0.0
	d(3,1)=0.0
	d(3,2)=0.0
	d(3,3)=e/(2*(1.0+mu))
	return
	end
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C############### Subroutine to make zero matrix
	subroutine get_matzero(A,n,m)
	INCLUDE 'ABA_PARAM.INC'
	real(8) A(n,m)
	do i = 1, n
		do j = 1, m
			A(i,j)=0.d0
		enddo
	enddo
	return
	end
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C############### Subroutine to multiply matrix
	subroutine get_matmul(A,B,C,l,n,m)
	INCLUDE 'ABA_PARAM.INC'
	real(8) A(l,n),B(n,m),C(l,m)
	call get_matzero(C,l,m)
	do i = 1, l
		do j = 1, m
			do k = 1, n
				C(i,j)=c(i,j)+A(i,k)*B(k,j)
			enddo
		enddo
	enddo
	return
	end
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C############### Subroutine to Calculate KIJ
	subroutine get_kmat(A,B,C,k,l,m,n,area,D)
	INCLUDE 'ABA_PARAM.INC'
	real(8) A(k,l),B(l,m),C(m,n)
	real(8) D(k,n),E(k,m)
	real(8) area
	call get_matzero(D,k,n)
	call get_matzero(E,k,m)
	call get_matmul(A,B,E,k,l,m)
	call get_matmul(E,C,D,k,m,n)
	do i=1,k
		do j=1,n
			D(i,j)=D(i,j)*area
		enddo
	enddo
	return
	end
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C############### Subroutine to do matrix transpose
	subroutine get_mattran(A,B,m,n)
	INCLUDE 'ABA_PARAM.INC'
	real(8) A(m,n),B(n,m)
	call get_matzero(B,n,m)
	do i=1,n
		do j=1,m
			B(i,j)=A(j,i)
		enddo
	enddo
	return
	end
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------