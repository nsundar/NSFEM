C User subroutine for NSFEM method with base mesh of triangular elements
C RHS is calculated using AMATRIX and U multiplication
C NDOFEL is used to make each subroutine consistent
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
     	open(10,file='D:\PhD_work\Work_Jan_June_2018\Oct_2018\ABQ_UEL\
     1NSFEM_UEL\Bevel_plt\mesh.dat')
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
C This is for a single node(may be at the corner of global mesh) having three surounding nodes 
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
C This is for a node on the boundries or corners or inside the global mesh surrounded by nine
C nodes or having more than one adjacent elements
         call U7(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
	elseif(jtype.eq.8) then
C This is for a node on the boundries or corners or inside the global mesh surrounded by ten
C nodes or having more than one adjacent elements
         call U8(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
	elseif(jtype.eq.9) then
C This is for a node on the boundries or corners or inside the global mesh surrounded by eleven
C nodes or having more than one adjacent elements
         call U9(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
	elseif(jtype.eq.10) then
C This is for a node on the boundries or corners or inside the global mesh surrounded by twelve
C nodes or having more than one adjacent elements
         call U10(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
	1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
	2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
	3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
	4     PERIOD,numnode,numelem,node,element)
	elseif(jtype.eq.11) then
C This is for a node on the boundries or corners or inside the global mesh surrounded by thirteen
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
         call exit
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
! NArea is the area of each node cell.
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Define parameters, dimensions and element type
		parameter(NEG=-1.D0,zero=0.D0, half=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
		INTEGER :: numnode,numelem,Ncen
		REAL(8) :: node(numnode,2)
		INTEGER :: element(numelem,3),Nodnum(NDOFEL*half) 
C
		REAL(8) NN(6),NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2)
        REAL(8) N3(2),NArea,A11,B1(3,6),UNL(NDOFEL),detCC
		REAL(8) Bl(3,NDOFEL),pb(3,3),S(4,4),GT(3,NDOFEL),mu1
        REAL(8) lam,mu,F(3,3),CC(3,3),INVCC(3,3),dd1(3,3,3,3)
        REAL(8) CMAT(3,3),Fbar(3,4),BB(NDOFEL,4),pbs(3,1)
        REAL(8) GC(NDOFEL,3),GCGT(NDOFEL,NDOFEL),BBS(NDOFEL,4)
        REAL(8) BBSBT(NDOFEL,NDOFEL),GB(NDOFEL,1)
        REAL(8) KSMATRX(NDOFEL,1),KTMATRX(NDOFEL,NDOFEL) 
C
        INTEGER i,j,kk1,kk2,i2,nn1,nn2,nn3,x1,j1,k1,l1
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
		AMATRX=zero
		Bl=zero
        KSMATRX=zero
        KTMATRX=zero
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
    	do I=1,NDOFEL
         UNL(I) = U(I)
         RHS(I,1)=ZERO
        enddo
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
C     write (*,*) Ncen
C     write (*,*) COORDS
C     write (*,*) node
	do i = 1,(NDOFEL*half)
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo
C     write (*,*) Nodnum
	NArea=ZERO
C---------- Assemble Constitutive matrix---D matrix-----C-------C-------C-------C-------
C---------- Pass material properties into variables-----C-------C-------C-------C-------
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition------C-------C-------C-------C-------C-------C-------
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data
! Loop for all elements		
	do ie = 1,numelem
	  A11=ZERO
	  N1=ZERO
	  N2=ZERO
	  N3=ZERO
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
	  call get_NArea_NBmat(N1,N2,N3,A11,B1)
	  do i=1,3
		    do j = 1,(NDOFEL*half)
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			  Bl(ii,(2*(j-1)+k)) = Bl(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
	  enddo
		NArea=NArea+A11/3
	  else
	  endif
	enddo
C	
	do i = 1,3
		do j = 1,NDOFEL
			Bl(i,j) = (1/NArea)*Bl(i,j)
		enddo
	enddo
c
c Calculation of initial stiffness and tangent stifness matrices
c
        pb=zero
		S=zero
        Fbar=zero
        BB=zero
		GT=zero
		GC=zero
        GCGT=zero
        BBS=zero
        BBSBT=zero
        GB=zero
        pbs=zero
        lam=zero
        mu=zero
        mu1=zero
        dd1=zero
        CMAT=zero
        F=zero
        detF=zero
        CC=zero
        INVCC=zero
        detCC=zero
        dudx = zero
        dudy = zero
        dvdx = zero
        dvdy = zero        
c
		do j = 1,(NDOFEL*half)
			dudx = dudx+Bl(1,(2*j-1))*UNL(2*j-1)
            dudy = dudy+Bl(2,2*j)*UNL(2*j-1)
            dvdx = dvdx+Bl(1,(2*j-1))*UNL(2*j)
            dvdy = dvdy+Bl(2,2*j)*UNL(2*j)
		enddo
c  
		lam=(EMU*E)/((1+EMU)*(1-2*EMU))	
        mu=E/(2*(1+EMU))	

		F(1,1)=dudx+1
        F(1,2)=dudy
        F(2,1)=dvdx
        F(2,2)=dvdy+1
        F(3,3)=1

        CC(1,1)=F(1,1)*F(1,1)+F(2,1)*F(2,1)
        CC(1,2)=F(1,2)*F(1,1)+F(2,1)*F(2,2)
        CC(2,1)=CC(1,2)
        CC(2,2)=F(1,2)*F(1,2)+F(2,2)*F(2,2)
        CC(3,3)=1

        detCC=(CC(1,1)*CC(2,2))-(CC(1,2)*CC(2,1))

        INVCC(1,1)=CC(2,2)/detCC
        INVCC(1,2)=-CC(1,2)/detCC
        INVCC(2,1)=-CC(2,1)/detCC
        INVCC(2,2)=CC(1,1)/detCC
        INVCC(3,3)=1
        
		mu1=mu-lam*half*log(detCC)

	   	DO ii=1,3
	  	DO jj=1,3
        	if (ii .EQ. jj)  then
	   		pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*(INVCC(ii,jj)-1)
        	else
            pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*INVCC(ii,jj)
          	endif
	  	ENDDO
        ENDDO
        
C Define CMAT matrix for Neohookian material model  
		
		DO x1=1,3
	  	DO j1=1,3
		DO k1=1,3
	  	DO l1=1,3
        dd1(x1,j1,k1,l1)=lam*INVCC(x1,j1)*INVCC(k1,l1)+
     1mu1*(INVCC(x1,k1)*INVCC(j1,l1)+INVCC(x1,l1)*INVCC(j1,k1))
	  	ENDDO
        ENDDO
        ENDDO
        ENDDO
		
	  	CMAT(1,1)=dd1(1,1,1,1)
	  	CMAT(1,2)=dd1(1,1,2,2)
	  	CMAT(1,3)=dd1(1,1,1,2)
	  	CMAT(2,1)=dd1(2,2,1,1)
      	CMAT(2,2)=dd1(2,2,2,2)
	  	CMAT(2,3)=dd1(2,2,1,2)  
      	CMAT(3,1)=dd1(1,2,1,1)
      	CMAT(3,2)=dd1(1,2,2,2)
	  	CMAT(3,3)=dd1(1,2,1,2)
c        
		S(1,1)=pb(1,1)
 		S(1,2)=pb(1,2)
        S(2,1)=pb(2,1)
        S(2,2)=pb(2,2)
        S(3,3)=pb(1,1)
        S(3,4)=pb(1,2)
        S(4,3)=pb(2,1)
        S(4,4)=pb(2,2)

        Fbar(1,1)=F(1,1)
 		Fbar(1,3)=F(2,1)
        Fbar(2,2)=F(1,2)
        Fbar(2,4)=F(2,2)
        Fbar(3,1)=F(1,2)
        Fbar(3,2)=F(1,1)
        Fbar(3,3)=F(2,2)
        Fbar(3,4)=F(2,1)
c
        DO i=1,(NDOFEL*half)
		BB((2*i-1),1)=Bl(1,(2*i-1))
        BB((2*i-1),2)=Bl(2,(2*i))
        BB((2*i),3)=Bl(1,(2*i-1))
        BB((2*i),4)=Bl(2,(2*i))
        ENDDO
c		
		DO i=1,3
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	GT(i,j)=GT(i,j)+Fbar(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
        pbs(1,1)=pb(1,1)
        pbs(2,1)=pb(2,2)
        pbs(3,1)=pb(1,2)
c
		DO i=1,NDOFEL
	  	DO j=1,3
        DO kk1=1,3
	   	GC(i,j)=GC(i,j)+GT(kk1,i)*CMAT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,3
	   	GCGT(i,j)=GCGT(i,j)+GC(i,kk1)*GT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c
	  	DO i=1,NDOFEL
        DO j=1,4
        DO kk1=1,4
	   	BBS(i,j)=BBS(i,j)+BB(i,kk1)*S(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	BBSBT(i,j)=BBSBT(i,j)+BBS(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
		DO i=1,NDOFEL
        DO kk2=1,3
	   	GB(i,1)=GB(i,1)+GT(kk2,i)*pbs(kk2,1)
	  	ENDDO
        ENDDO
c
      	DO i=1,NDOFEL
        KSMATRX(i,1)=KSMATRX(i,1)+GB(i,1)
        DO j=1,NDOFEL
		KTMATRX(i,j)=KTMATRX(i,j)+GCGT(i,j)+BBSBT(i,j)
        ENDDO
		ENDDO
C===============================C
C        calculate  RHS         C
C===============================C
  	  DO K1=1, NDOFEL              !***!   
	    RHS(K1, 1)=RHS(K1, 1)-NArea*KSMATRX(K1,1)
	  ENDDO
c
  	  DO K1=1,NDOFEL              !***!  DEFINING K MATRIX 
	     DO K2=1,NDOFEL              
	     AMATRX(K1, K2)=AMATRX(K1, K2)+NArea*KTMATRX(K1,K2)
	     ENDDO
	  ENDDO       
C---------- Output the saved results           
C         write (*,*) 'Node number',JELEM
c         write (*,*) 'Bmat matrix', GT
c         write (*,*) 'Bgeo matrix transpose', BB 
c         write (*,*) 'F matrix', F
c		  write (*,*) 'dd1 matrix', dd1
c         write (*,*) 'pb matrix', pb
c         write (*,*) 'Cmat matrix', CMAT
c         write (*,*) 'KSMATRX matrix', KSMATRX
c         write (*,*) 'KTMATRX matrix', KTMATRX
c         write (*,*) 'Displacement values', U
C           pause
	return
	end subroutine U1      
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
!This is for a single node having four surrounding nodes including the node itself
!UEL actually calculates nodal stiffness; not the element stiffness
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
		parameter(NEG=-1.D0,zero=0.D0, half=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
		INTEGER :: numnode,numelem,Ncen
		REAL(8) :: node(numnode,2)
		INTEGER :: element(numelem,3),Nodnum(NDOFEL*half) 
C
		REAL(8) NN(6),NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2)
        REAL(8) N3(2),NArea,A11,B1(3,6),UNL(NDOFEL),detCC
		REAL(8) Bl(3,NDOFEL),pb(3,3),S(4,4),GT(3,NDOFEL),mu1
        REAL(8) lam,mu,F(3,3),CC(3,3),INVCC(3,3),dd1(3,3,3,3)
        REAL(8) CMAT(3,3),Fbar(3,4),BB(NDOFEL,4),pbs(3,1)
        REAL(8) GC(NDOFEL,3),GCGT(NDOFEL,NDOFEL),BBS(NDOFEL,4)
        REAL(8) BBSBT(NDOFEL,NDOFEL),GB(NDOFEL,1)
        REAL(8) KSMATRX(NDOFEL,1),KTMATRX(NDOFEL,NDOFEL) 
C
        INTEGER i,j,kk1,kk2,i2,nn1,nn2,nn3,x1,j1,k1,l1
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
		AMATRX=zero
		Bl=zero
        KSMATRX=zero
        KTMATRX=zero
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
    	do I=1,NDOFEL
         UNL(I) = U(I)
         RHS(I,1)=ZERO
        enddo
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,(NDOFEL*half)
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=ZERO
C---------- Assemble Constitutive matrix---D matrix-----C-------C-------C-------C-------
C---------- Pass material properties into variables-----C-------C-------C-------C-------
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition------C-------C-------C-------C-------C-------C-------
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data
! Loop for all elements		
	do ie = 1,numelem
	  A11=ZERO
	  N1=ZERO
	  N2=ZERO
	  N3=ZERO
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
	  call get_NArea_NBmat(N1,N2,N3,A11,B1)
	  do i=1,3
		    do j = 1,(NDOFEL*half)
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			  Bl(ii,(2*(j-1)+k)) = Bl(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
	  enddo
		NArea=NArea+A11/3
	  else
	  endif
	enddo
C	
	do i = 1,3
		do j = 1,NDOFEL
			Bl(i,j) = (1/NArea)*Bl(i,j)
		enddo
	enddo
c
c Calculation of initial stiffness and tangent stifness matrices
c
        pb=zero
		S=zero
        Fbar=zero
        BB=zero
		GT=zero
		GC=zero
        GCGT=zero
        BBS=zero
        BBSBT=zero
        GB=zero
        pbs=zero
        lam=zero
        mu=zero
        mu1=zero
        dd1=zero
        CMAT=zero
        F=zero
        detF=zero
        CC=zero
        INVCC=zero
        detCC=zero
        dudx = zero
        dudy = zero
        dvdx = zero
        dvdy = zero        
c
		do j = 1,(NDOFEL*half)
			dudx = dudx+Bl(1,(2*j-1))*UNL(2*j-1)
            dudy = dudy+Bl(2,2*j)*UNL(2*j-1)
            dvdx = dvdx+Bl(1,(2*j-1))*UNL(2*j)
            dvdy = dvdy+Bl(2,2*j)*UNL(2*j)
		enddo
c  
		lam=(EMU*E)/((1+EMU)*(1-2*EMU))	
        mu=E/(2*(1+EMU))	

		F(1,1)=dudx+1
        F(1,2)=dudy
        F(2,1)=dvdx
        F(2,2)=dvdy+1
        F(3,3)=1

        CC(1,1)=F(1,1)*F(1,1)+F(2,1)*F(2,1)
        CC(1,2)=F(1,2)*F(1,1)+F(2,1)*F(2,2)
        CC(2,1)=CC(1,2)
        CC(2,2)=F(1,2)*F(1,2)+F(2,2)*F(2,2)
        CC(3,3)=1

        detCC=(CC(1,1)*CC(2,2))-(CC(1,2)*CC(2,1))

        INVCC(1,1)=CC(2,2)/detCC
        INVCC(1,2)=-CC(1,2)/detCC
        INVCC(2,1)=-CC(2,1)/detCC
        INVCC(2,2)=CC(1,1)/detCC
        INVCC(3,3)=1
        
		mu1=mu-lam*half*log(detCC)

	   	DO ii=1,3
	  	DO jj=1,3
        	if (ii .EQ. jj)  then
	   		pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*(INVCC(ii,jj)-1)
        	else
            pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*INVCC(ii,jj)
          	endif
	  	ENDDO
        ENDDO
        
C Define CMAT matrix for Neohookian material model  
		
		DO x1=1,3
	  	DO j1=1,3
		DO k1=1,3
	  	DO l1=1,3
        dd1(x1,j1,k1,l1)=lam*INVCC(x1,j1)*INVCC(k1,l1)+
     1mu1*(INVCC(x1,k1)*INVCC(j1,l1)+INVCC(x1,l1)*INVCC(j1,k1))
	  	ENDDO
        ENDDO
        ENDDO
        ENDDO
		
	  	CMAT(1,1)=dd1(1,1,1,1)
	  	CMAT(1,2)=dd1(1,1,2,2)
	  	CMAT(1,3)=dd1(1,1,1,2)
	  	CMAT(2,1)=dd1(2,2,1,1)
      	CMAT(2,2)=dd1(2,2,2,2)
	  	CMAT(2,3)=dd1(2,2,1,2)  
      	CMAT(3,1)=dd1(1,2,1,1)
      	CMAT(3,2)=dd1(1,2,2,2)
	  	CMAT(3,3)=dd1(1,2,1,2)
c        
		S(1,1)=pb(1,1)
 		S(1,2)=pb(1,2)
        S(2,1)=pb(2,1)
        S(2,2)=pb(2,2)
        S(3,3)=pb(1,1)
        S(3,4)=pb(1,2)
        S(4,3)=pb(2,1)
        S(4,4)=pb(2,2)

        Fbar(1,1)=F(1,1)
 		Fbar(1,3)=F(2,1)
        Fbar(2,2)=F(1,2)
        Fbar(2,4)=F(2,2)
        Fbar(3,1)=F(1,2)
        Fbar(3,2)=F(1,1)
        Fbar(3,3)=F(2,2)
        Fbar(3,4)=F(2,1)
c
        DO i=1,(NDOFEL*half)
		BB((2*i-1),1)=Bl(1,(2*i-1))
        BB((2*i-1),2)=Bl(2,(2*i))
        BB((2*i),3)=Bl(1,(2*i-1))
        BB((2*i),4)=Bl(2,(2*i))
        ENDDO
c		
		DO i=1,3
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	GT(i,j)=GT(i,j)+Fbar(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
        pbs(1,1)=pb(1,1)
        pbs(2,1)=pb(2,2)
        pbs(3,1)=pb(1,2)
c
		DO i=1,NDOFEL
	  	DO j=1,3
        DO kk1=1,3
	   	GC(i,j)=GC(i,j)+GT(kk1,i)*CMAT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,3
	   	GCGT(i,j)=GCGT(i,j)+GC(i,kk1)*GT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c
	  	DO i=1,NDOFEL
        DO j=1,4
        DO kk1=1,4
	   	BBS(i,j)=BBS(i,j)+BB(i,kk1)*S(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	BBSBT(i,j)=BBSBT(i,j)+BBS(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
		DO i=1,NDOFEL
        DO kk2=1,3
	   	GB(i,1)=GB(i,1)+GT(kk2,i)*pbs(kk2,1)
	  	ENDDO
        ENDDO
c
      	DO i=1,NDOFEL
        KSMATRX(i,1)=KSMATRX(i,1)+GB(i,1)
        DO j=1,NDOFEL
		KTMATRX(i,j)=KTMATRX(i,j)+GCGT(i,j)+BBSBT(i,j)
        ENDDO
		ENDDO
C===============================C
C        calculate  RHS         C
C===============================C
  	  DO K1=1, NDOFEL              !***!   
	    RHS(K1, 1)=RHS(K1, 1)-NArea*KSMATRX(K1,1)
	  ENDDO
c
  	  DO K1=1,NDOFEL              !***!  DEFINING K MATRIX 
	     DO K2=1,NDOFEL              
	     AMATRX(K1, K2)=AMATRX(K1, K2)+NArea*KTMATRX(K1,K2)
	     ENDDO
	  ENDDO       
C---------- Output the saved results           
c         write (*,*) 'Node number',JELEM
c         write (*,*) 'Bmat matrix', GT
c         write (*,*) 'Bgeo matrix transpose', BB 
c         write (*,*) 'F matrix', F
c		  write (*,*) 'dd1 matrix', dd1
c         write (*,*) 'pb matrix', pb
c         write (*,*) 'Cmat matrix', CMAT
c         write (*,*) 'KSMATRX matrix', KSMATRX
c         write (*,*) 'KTMATRX matrix', KTMATRX
c         write (*,*) 'Displacement values', U
c         pause
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
!This is for a single node having five surrounding nodes including the node itself
!UEL actually calculates nodal stiffness; not the element stiffness
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
		parameter(NEG=-1.D0,zero=0.D0, half=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
		INTEGER :: numnode,numelem,Ncen
		REAL(8) :: node(numnode,2)
		INTEGER :: element(numelem,3),Nodnum(NDOFEL*half) 
C
		REAL(8) NN(6),NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2)
        REAL(8) N3(2),NArea,A11,B1(3,6),UNL(NDOFEL),detCC
		REAL(8) Bl(3,NDOFEL),pb(3,3),S(4,4),GT(3,NDOFEL),mu1
        REAL(8) lam,mu,F(3,3),CC(3,3),INVCC(3,3),dd1(3,3,3,3)
        REAL(8) CMAT(3,3),Fbar(3,4),BB(NDOFEL,4),pbs(3,1)
        REAL(8) GC(NDOFEL,3),GCGT(NDOFEL,NDOFEL),BBS(NDOFEL,4)
        REAL(8) BBSBT(NDOFEL,NDOFEL),GB(NDOFEL,1)
        REAL(8) KSMATRX(NDOFEL,1),KTMATRX(NDOFEL,NDOFEL) 
C
        INTEGER i,j,kk1,kk2,i2,nn1,nn2,nn3,x1,j1,k1,l1
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
		AMATRX=zero
		Bl=zero
        KSMATRX=zero
        KTMATRX=zero
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
    	do I=1,NDOFEL
         UNL(I) = U(I)
         RHS(I,1)=ZERO
        enddo
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,(NDOFEL*half)
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=ZERO
C---------- Assemble Constitutive matrix---D matrix-----C-------C-------C-------C-------
C---------- Pass material properties into variables-----C-------C-------C-------C-------
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition------C-------C-------C-------C-------C-------C-------
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data
! Loop for all elements		
	do ie = 1,numelem
	  A11=ZERO
	  N1=ZERO
	  N2=ZERO
	  N3=ZERO
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
	  call get_NArea_NBmat(N1,N2,N3,A11,B1)
	  do i=1,3
		    do j = 1,(NDOFEL*half)
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			  Bl(ii,(2*(j-1)+k)) = Bl(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
	  enddo
		NArea=NArea+A11/3
	  else
	  endif
	enddo
C	
	do i = 1,3
		do j = 1,NDOFEL
			Bl(i,j) = (1/NArea)*Bl(i,j)
		enddo
	enddo
c
c Calculation of initial stiffness and tangent stifness matrices
c
        pb=zero
		S=zero
        Fbar=zero
        BB=zero
		GT=zero
		GC=zero
        GCGT=zero
        BBS=zero
        BBSBT=zero
        GB=zero
        pbs=zero
        lam=zero
        mu=zero
        mu1=zero
        dd1=zero
        CMAT=zero
        F=zero
        detF=zero
        CC=zero
        INVCC=zero
        detCC=zero
        dudx = zero
        dudy = zero
        dvdx = zero
        dvdy = zero        
c
		do j = 1,(NDOFEL*half)
			dudx = dudx+Bl(1,(2*j-1))*UNL(2*j-1)
            dudy = dudy+Bl(2,2*j)*UNL(2*j-1)
            dvdx = dvdx+Bl(1,(2*j-1))*UNL(2*j)
            dvdy = dvdy+Bl(2,2*j)*UNL(2*j)
		enddo
c  
		lam=(EMU*E)/((1+EMU)*(1-2*EMU))	
        mu=E/(2*(1+EMU))	

		F(1,1)=dudx+1
        F(1,2)=dudy
        F(2,1)=dvdx
        F(2,2)=dvdy+1
        F(3,3)=1

        CC(1,1)=F(1,1)*F(1,1)+F(2,1)*F(2,1)
        CC(1,2)=F(1,2)*F(1,1)+F(2,1)*F(2,2)
        CC(2,1)=CC(1,2)
        CC(2,2)=F(1,2)*F(1,2)+F(2,2)*F(2,2)
        CC(3,3)=1

        detCC=(CC(1,1)*CC(2,2))-(CC(1,2)*CC(2,1))

        INVCC(1,1)=CC(2,2)/detCC
        INVCC(1,2)=-CC(1,2)/detCC
        INVCC(2,1)=-CC(2,1)/detCC
        INVCC(2,2)=CC(1,1)/detCC
        INVCC(3,3)=1
        
		mu1=mu-lam*half*log(detCC)

	   	DO ii=1,3
	  	DO jj=1,3
        	if (ii .EQ. jj)  then
	   		pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*(INVCC(ii,jj)-1)
        	else
            pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*INVCC(ii,jj)
          	endif
	  	ENDDO
        ENDDO
        
C Define CMAT matrix for Neohookian material model  
		
		DO x1=1,3
	  	DO j1=1,3
		DO k1=1,3
	  	DO l1=1,3
        dd1(x1,j1,k1,l1)=lam*INVCC(x1,j1)*INVCC(k1,l1)+
     1mu1*(INVCC(x1,k1)*INVCC(j1,l1)+INVCC(x1,l1)*INVCC(j1,k1))
	  	ENDDO
        ENDDO
        ENDDO
        ENDDO
		
	  	CMAT(1,1)=dd1(1,1,1,1)
	  	CMAT(1,2)=dd1(1,1,2,2)
	  	CMAT(1,3)=dd1(1,1,1,2)
	  	CMAT(2,1)=dd1(2,2,1,1)
      	CMAT(2,2)=dd1(2,2,2,2)
	  	CMAT(2,3)=dd1(2,2,1,2)  
      	CMAT(3,1)=dd1(1,2,1,1)
      	CMAT(3,2)=dd1(1,2,2,2)
	  	CMAT(3,3)=dd1(1,2,1,2)
c        
		S(1,1)=pb(1,1)
 		S(1,2)=pb(1,2)
        S(2,1)=pb(2,1)
        S(2,2)=pb(2,2)
        S(3,3)=pb(1,1)
        S(3,4)=pb(1,2)
        S(4,3)=pb(2,1)
        S(4,4)=pb(2,2)

        Fbar(1,1)=F(1,1)
 		Fbar(1,3)=F(2,1)
        Fbar(2,2)=F(1,2)
        Fbar(2,4)=F(2,2)
        Fbar(3,1)=F(1,2)
        Fbar(3,2)=F(1,1)
        Fbar(3,3)=F(2,2)
        Fbar(3,4)=F(2,1)
c
        DO i=1,(NDOFEL*half)
		BB((2*i-1),1)=Bl(1,(2*i-1))
        BB((2*i-1),2)=Bl(2,(2*i))
        BB((2*i),3)=Bl(1,(2*i-1))
        BB((2*i),4)=Bl(2,(2*i))
        ENDDO
c		
		DO i=1,3
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	GT(i,j)=GT(i,j)+Fbar(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
        pbs(1,1)=pb(1,1)
        pbs(2,1)=pb(2,2)
        pbs(3,1)=pb(1,2)
c
		DO i=1,NDOFEL
	  	DO j=1,3
        DO kk1=1,3
	   	GC(i,j)=GC(i,j)+GT(kk1,i)*CMAT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,3
	   	GCGT(i,j)=GCGT(i,j)+GC(i,kk1)*GT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c
	  	DO i=1,NDOFEL
        DO j=1,4
        DO kk1=1,4
	   	BBS(i,j)=BBS(i,j)+BB(i,kk1)*S(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	BBSBT(i,j)=BBSBT(i,j)+BBS(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
		DO i=1,NDOFEL
        DO kk2=1,3
	   	GB(i,1)=GB(i,1)+GT(kk2,i)*pbs(kk2,1)
	  	ENDDO
        ENDDO
c
      	DO i=1,NDOFEL
        KSMATRX(i,1)=KSMATRX(i,1)+GB(i,1)
        DO j=1,NDOFEL
		KTMATRX(i,j)=KTMATRX(i,j)+GCGT(i,j)+BBSBT(i,j)
        ENDDO
		ENDDO
C===============================C
C        calculate  RHS         C
C===============================C
  	  DO K1=1, NDOFEL              !***!   
	    RHS(K1, 1)=RHS(K1, 1)-NArea*KSMATRX(K1,1)
	  ENDDO
c
  	  DO K1=1,NDOFEL              !***!  DEFINING K MATRIX 
	     DO K2=1,NDOFEL              
	     AMATRX(K1, K2)=AMATRX(K1, K2)+NArea*KTMATRX(K1,K2)
	     ENDDO
	  ENDDO       
C---------- Output the saved results           
c         write (*,*) 'Node number',JELEM
c         write (*,*) 'Bmat matrix', GT
c         write (*,*) 'Bgeo matrix transpose', BB 
c         write (*,*) 'F matrix', F
c		  write (*,*) 'dd1 matrix', dd1
c         write (*,*) 'pb matrix', pb
c         write (*,*) 'Cmat matrix', CMAT
c         write (*,*) 'KSMATRX matrix', KSMATRX
c         write (*,*) 'KTMATRX matrix', KTMATRX
c         write (*,*) 'Displacement values', U
c         pause
	return
	end subroutine U3 
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
!This is for a single node having six surrounding nodes including the node itself
!UEL actually calculates nodal stiffness; not the element stiffness
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
		parameter(NEG=-1.D0,zero=0.D0, half=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
		INTEGER :: numnode,numelem,Ncen
		REAL(8) :: node(numnode,2)
		INTEGER :: element(numelem,3),Nodnum(NDOFEL*half) 
C
		REAL(8) NN(6),NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2)
        REAL(8) N3(2),NArea,A11,B1(3,6),UNL(NDOFEL),detCC
		REAL(8) Bl(3,NDOFEL),pb(3,3),S(4,4),GT(3,NDOFEL),mu1
        REAL(8) lam,mu,F(3,3),CC(3,3),INVCC(3,3),dd1(3,3,3,3)
        REAL(8) CMAT(3,3),Fbar(3,4),BB(NDOFEL,4),pbs(3,1)
        REAL(8) GC(NDOFEL,3),GCGT(NDOFEL,NDOFEL),BBS(NDOFEL,4)
        REAL(8) BBSBT(NDOFEL,NDOFEL),GB(NDOFEL,1)
        REAL(8) KSMATRX(NDOFEL,1),KTMATRX(NDOFEL,NDOFEL) 
C
        INTEGER i,j,kk1,kk2,i2,nn1,nn2,nn3,x1,j1,k1,l1
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
		AMATRX=zero
		Bl=zero
        KSMATRX=zero
        KTMATRX=zero
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
    	do I=1,NDOFEL
         UNL(I) = U(I)
         RHS(I,1)=ZERO
        enddo
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,(NDOFEL*half)
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=ZERO
C---------- Assemble Constitutive matrix---D matrix-----C-------C-------C-------C-------
C---------- Pass material properties into variables-----C-------C-------C-------C-------
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition------C-------C-------C-------C-------C-------C-------
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data
! Loop for all elements		
	do ie = 1,numelem
	  A11=ZERO
	  N1=ZERO
	  N2=ZERO
	  N3=ZERO
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
	  call get_NArea_NBmat(N1,N2,N3,A11,B1)
	  do i=1,3
		    do j = 1,(NDOFEL*half)
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			  Bl(ii,(2*(j-1)+k)) = Bl(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
	  enddo
		NArea=NArea+A11/3
	  else
	  endif
	enddo
C	
	do i = 1,3
		do j = 1,NDOFEL
			Bl(i,j) = (1/NArea)*Bl(i,j)
		enddo
	enddo
c
c Calculation of initial stiffness and tangent stifness matrices
c
        pb=zero
		S=zero
        Fbar=zero
        BB=zero
		GT=zero
		GC=zero
        GCGT=zero
        BBS=zero
        BBSBT=zero
        GB=zero
        pbs=zero
        lam=zero
        mu=zero
        mu1=zero
        dd1=zero
        CMAT=zero
        F=zero
        detF=zero
        CC=zero
        INVCC=zero
        detCC=zero
        dudx = zero
        dudy = zero
        dvdx = zero
        dvdy = zero        
c
		do j = 1,(NDOFEL*half)
			dudx = dudx+Bl(1,(2*j-1))*UNL(2*j-1)
            dudy = dudy+Bl(2,2*j)*UNL(2*j-1)
            dvdx = dvdx+Bl(1,(2*j-1))*UNL(2*j)
            dvdy = dvdy+Bl(2,2*j)*UNL(2*j)
		enddo
c  
		lam=(EMU*E)/((1+EMU)*(1-2*EMU))	
        mu=E/(2*(1+EMU))	

		F(1,1)=dudx+1
        F(1,2)=dudy
        F(2,1)=dvdx
        F(2,2)=dvdy+1
        F(3,3)=1

        CC(1,1)=F(1,1)*F(1,1)+F(2,1)*F(2,1)
        CC(1,2)=F(1,2)*F(1,1)+F(2,1)*F(2,2)
        CC(2,1)=CC(1,2)
        CC(2,2)=F(1,2)*F(1,2)+F(2,2)*F(2,2)
        CC(3,3)=1

        detCC=(CC(1,1)*CC(2,2))-(CC(1,2)*CC(2,1))

        INVCC(1,1)=CC(2,2)/detCC
        INVCC(1,2)=-CC(1,2)/detCC
        INVCC(2,1)=-CC(2,1)/detCC
        INVCC(2,2)=CC(1,1)/detCC
        INVCC(3,3)=1
        
		mu1=mu-lam*half*log(detCC)

	   	DO ii=1,3
	  	DO jj=1,3
        	if (ii .EQ. jj)  then
	   		pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*(INVCC(ii,jj)-1)
        	else
            pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*INVCC(ii,jj)
          	endif
	  	ENDDO
        ENDDO
        
C Define CMAT matrix for Neohookian material model  
		
		DO x1=1,3
	  	DO j1=1,3
		DO k1=1,3
	  	DO l1=1,3
        dd1(x1,j1,k1,l1)=lam*INVCC(x1,j1)*INVCC(k1,l1)+
     1mu1*(INVCC(x1,k1)*INVCC(j1,l1)+INVCC(x1,l1)*INVCC(j1,k1))
	  	ENDDO
        ENDDO
        ENDDO
        ENDDO
		
	  	CMAT(1,1)=dd1(1,1,1,1)
	  	CMAT(1,2)=dd1(1,1,2,2)
	  	CMAT(1,3)=dd1(1,1,1,2)
	  	CMAT(2,1)=dd1(2,2,1,1)
      	CMAT(2,2)=dd1(2,2,2,2)
	  	CMAT(2,3)=dd1(2,2,1,2)  
      	CMAT(3,1)=dd1(1,2,1,1)
      	CMAT(3,2)=dd1(1,2,2,2)
	  	CMAT(3,3)=dd1(1,2,1,2)
c        
		S(1,1)=pb(1,1)
 		S(1,2)=pb(1,2)
        S(2,1)=pb(2,1)
        S(2,2)=pb(2,2)
        S(3,3)=pb(1,1)
        S(3,4)=pb(1,2)
        S(4,3)=pb(2,1)
        S(4,4)=pb(2,2)

        Fbar(1,1)=F(1,1)
 		Fbar(1,3)=F(2,1)
        Fbar(2,2)=F(1,2)
        Fbar(2,4)=F(2,2)
        Fbar(3,1)=F(1,2)
        Fbar(3,2)=F(1,1)
        Fbar(3,3)=F(2,2)
        Fbar(3,4)=F(2,1)
c
        DO i=1,(NDOFEL*half)
		BB((2*i-1),1)=Bl(1,(2*i-1))
        BB((2*i-1),2)=Bl(2,(2*i))
        BB((2*i),3)=Bl(1,(2*i-1))
        BB((2*i),4)=Bl(2,(2*i))
        ENDDO
c		
		DO i=1,3
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	GT(i,j)=GT(i,j)+Fbar(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
        pbs(1,1)=pb(1,1)
        pbs(2,1)=pb(2,2)
        pbs(3,1)=pb(1,2)
c
		DO i=1,NDOFEL
	  	DO j=1,3
        DO kk1=1,3
	   	GC(i,j)=GC(i,j)+GT(kk1,i)*CMAT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,3
	   	GCGT(i,j)=GCGT(i,j)+GC(i,kk1)*GT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c
	  	DO i=1,NDOFEL
        DO j=1,4
        DO kk1=1,4
	   	BBS(i,j)=BBS(i,j)+BB(i,kk1)*S(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	BBSBT(i,j)=BBSBT(i,j)+BBS(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
		DO i=1,NDOFEL
        DO kk2=1,3
	   	GB(i,1)=GB(i,1)+GT(kk2,i)*pbs(kk2,1)
	  	ENDDO
        ENDDO
c
      	DO i=1,NDOFEL
        KSMATRX(i,1)=KSMATRX(i,1)+GB(i,1)
        DO j=1,NDOFEL
		KTMATRX(i,j)=KTMATRX(i,j)+GCGT(i,j)+BBSBT(i,j)
        ENDDO
		ENDDO
C===============================C
C        calculate  RHS         C
C===============================C
  	  DO K1=1, NDOFEL              !***!   
	    RHS(K1, 1)=RHS(K1, 1)-NArea*KSMATRX(K1,1)
	  ENDDO
c
  	  DO K1=1,NDOFEL              !***!  DEFINING K MATRIX 
	     DO K2=1,NDOFEL              
	     AMATRX(K1, K2)=AMATRX(K1, K2)+NArea*KTMATRX(K1,K2)
	     ENDDO
	  ENDDO       
C---------- Output the saved results           
c         write (*,*) 'Node number',JELEM
c         write (*,*) 'Bmat matrix', GT
c         write (*,*) 'Bgeo matrix transpose', BB 
c         write (*,*) 'F matrix', F
c		  write (*,*) 'dd1 matrix', dd1
c         write (*,*) 'pb matrix', pb
c         write (*,*) 'Cmat matrix', CMAT
c         write (*,*) 'KSMATRX matrix', KSMATRX
c         write (*,*) 'KTMATRX matrix', KTMATRX
c         write (*,*) 'Displacement values', U
c         pause
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
!This is for a single node having seven surrounding nodes including the node itself
!UEL actually calculates nodal stiffness; not the element stiffness
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
		parameter(NEG=-1.D0,zero=0.D0, half=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
		INTEGER :: numnode,numelem,Ncen
		REAL(8) :: node(numnode,2)
		INTEGER :: element(numelem,3),Nodnum(NDOFEL*half) 
C
		REAL(8) NN(6),NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2)
        REAL(8) N3(2),NArea,A11,B1(3,6),UNL(NDOFEL),detCC
		REAL(8) Bl(3,NDOFEL),pb(3,3),S(4,4),GT(3,NDOFEL),mu1
        REAL(8) lam,mu,F(3,3),CC(3,3),INVCC(3,3),dd1(3,3,3,3)
        REAL(8) CMAT(3,3),Fbar(3,4),BB(NDOFEL,4),pbs(3,1)
        REAL(8) GC(NDOFEL,3),GCGT(NDOFEL,NDOFEL),BBS(NDOFEL,4)
        REAL(8) BBSBT(NDOFEL,NDOFEL),GB(NDOFEL,1)
        REAL(8) KSMATRX(NDOFEL,1),KTMATRX(NDOFEL,NDOFEL) 
C
        INTEGER i,j,kk1,kk2,i2,nn1,nn2,nn3,x1,j1,k1,l1
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
		AMATRX=zero
		Bl=zero
        KSMATRX=zero
        KTMATRX=zero
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
    	do I=1,NDOFEL
         UNL(I) = U(I)
         RHS(I,1)=ZERO
        enddo
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,(NDOFEL*half)
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=ZERO
C---------- Assemble Constitutive matrix---D matrix-----C-------C-------C-------C-------
C---------- Pass material properties into variables-----C-------C-------C-------C-------
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition------C-------C-------C-------C-------C-------C-------
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data
! Loop for all elements		
	do ie = 1,numelem
	  A11=ZERO
	  N1=ZERO
	  N2=ZERO
	  N3=ZERO
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
	  call get_NArea_NBmat(N1,N2,N3,A11,B1)
	  do i=1,3
		    do j = 1,(NDOFEL*half)
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			  Bl(ii,(2*(j-1)+k)) = Bl(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
	  enddo
		NArea=NArea+A11/3
	  else
	  endif
	enddo
C	
	do i = 1,3
		do j = 1,NDOFEL
			Bl(i,j) = (1/NArea)*Bl(i,j)
		enddo
	enddo
c
c Calculation of initial stiffness and tangent stifness matrices
c
        pb=zero
		S=zero
        Fbar=zero
        BB=zero
		GT=zero
		GC=zero
        GCGT=zero
        BBS=zero
        BBSBT=zero
        GB=zero
        pbs=zero
        lam=zero
        mu=zero
        mu1=zero
        dd1=zero
        CMAT=zero
        F=zero
        detF=zero
        CC=zero
        INVCC=zero
        detCC=zero
        dudx = zero
        dudy = zero
        dvdx = zero
        dvdy = zero        
c
		do j = 1,(NDOFEL*half)
			dudx = dudx+Bl(1,(2*j-1))*UNL(2*j-1)
            dudy = dudy+Bl(2,2*j)*UNL(2*j-1)
            dvdx = dvdx+Bl(1,(2*j-1))*UNL(2*j)
            dvdy = dvdy+Bl(2,2*j)*UNL(2*j)
		enddo
c  
		lam=(EMU*E)/((1+EMU)*(1-2*EMU))	
        mu=E/(2*(1+EMU))	

		F(1,1)=dudx+1
        F(1,2)=dudy
        F(2,1)=dvdx
        F(2,2)=dvdy+1
        F(3,3)=1

        CC(1,1)=F(1,1)*F(1,1)+F(2,1)*F(2,1)
        CC(1,2)=F(1,2)*F(1,1)+F(2,1)*F(2,2)
        CC(2,1)=CC(1,2)
        CC(2,2)=F(1,2)*F(1,2)+F(2,2)*F(2,2)
        CC(3,3)=1

        detCC=(CC(1,1)*CC(2,2))-(CC(1,2)*CC(2,1))

        INVCC(1,1)=CC(2,2)/detCC
        INVCC(1,2)=-CC(1,2)/detCC
        INVCC(2,1)=-CC(2,1)/detCC
        INVCC(2,2)=CC(1,1)/detCC
        INVCC(3,3)=1
        
		mu1=mu-lam*half*log(detCC)

	   	DO ii=1,3
	  	DO jj=1,3
        	if (ii .EQ. jj)  then
	   		pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*(INVCC(ii,jj)-1)
        	else
            pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*INVCC(ii,jj)
          	endif
	  	ENDDO
        ENDDO
        
C Define CMAT matrix for Neohookian material model  
		
		DO x1=1,3
	  	DO j1=1,3
		DO k1=1,3
	  	DO l1=1,3
        dd1(x1,j1,k1,l1)=lam*INVCC(x1,j1)*INVCC(k1,l1)+
     1mu1*(INVCC(x1,k1)*INVCC(j1,l1)+INVCC(x1,l1)*INVCC(j1,k1))
	  	ENDDO
        ENDDO
        ENDDO
        ENDDO
		
	  	CMAT(1,1)=dd1(1,1,1,1)
	  	CMAT(1,2)=dd1(1,1,2,2)
	  	CMAT(1,3)=dd1(1,1,1,2)
	  	CMAT(2,1)=dd1(2,2,1,1)
      	CMAT(2,2)=dd1(2,2,2,2)
	  	CMAT(2,3)=dd1(2,2,1,2)  
      	CMAT(3,1)=dd1(1,2,1,1)
      	CMAT(3,2)=dd1(1,2,2,2)
	  	CMAT(3,3)=dd1(1,2,1,2)
c        
		S(1,1)=pb(1,1)
 		S(1,2)=pb(1,2)
        S(2,1)=pb(2,1)
        S(2,2)=pb(2,2)
        S(3,3)=pb(1,1)
        S(3,4)=pb(1,2)
        S(4,3)=pb(2,1)
        S(4,4)=pb(2,2)

        Fbar(1,1)=F(1,1)
 		Fbar(1,3)=F(2,1)
        Fbar(2,2)=F(1,2)
        Fbar(2,4)=F(2,2)
        Fbar(3,1)=F(1,2)
        Fbar(3,2)=F(1,1)
        Fbar(3,3)=F(2,2)
        Fbar(3,4)=F(2,1)
c
        DO i=1,(NDOFEL*half)
		BB((2*i-1),1)=Bl(1,(2*i-1))
        BB((2*i-1),2)=Bl(2,(2*i))
        BB((2*i),3)=Bl(1,(2*i-1))
        BB((2*i),4)=Bl(2,(2*i))
        ENDDO
c		
		DO i=1,3
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	GT(i,j)=GT(i,j)+Fbar(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
        pbs(1,1)=pb(1,1)
        pbs(2,1)=pb(2,2)
        pbs(3,1)=pb(1,2)
c
		DO i=1,NDOFEL
	  	DO j=1,3
        DO kk1=1,3
	   	GC(i,j)=GC(i,j)+GT(kk1,i)*CMAT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,3
	   	GCGT(i,j)=GCGT(i,j)+GC(i,kk1)*GT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c
	  	DO i=1,NDOFEL
        DO j=1,4
        DO kk1=1,4
	   	BBS(i,j)=BBS(i,j)+BB(i,kk1)*S(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	BBSBT(i,j)=BBSBT(i,j)+BBS(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
		DO i=1,NDOFEL
        DO kk2=1,3
	   	GB(i,1)=GB(i,1)+GT(kk2,i)*pbs(kk2,1)
	  	ENDDO
        ENDDO
c
      	DO i=1,NDOFEL
        KSMATRX(i,1)=KSMATRX(i,1)+GB(i,1)
        DO j=1,NDOFEL
		KTMATRX(i,j)=KTMATRX(i,j)+GCGT(i,j)+BBSBT(i,j)
        ENDDO
		ENDDO
C===============================C
C        calculate  RHS         C
C===============================C
  	  DO K1=1, NDOFEL              !***!   
	    RHS(K1, 1)=RHS(K1, 1)-NArea*KSMATRX(K1,1)
	  ENDDO
c
  	  DO K1=1,NDOFEL              !***!  DEFINING K MATRIX 
	     DO K2=1,NDOFEL              
	     AMATRX(K1, K2)=AMATRX(K1, K2)+NArea*KTMATRX(K1,K2)
	     ENDDO
	  ENDDO       
C---------- Output the saved results           
c         write (*,*) 'Node number',JELEM
c         write (*,*) 'Bmat matrix', GT
c         write (*,*) 'Bgeo matrix transpose', BB 
c         write (*,*) 'F matrix', F
c		  write (*,*) 'dd1 matrix', dd1
c         write (*,*) 'pb matrix', pb
c         write (*,*) 'Cmat matrix', CMAT
c         write (*,*) 'KSMATRX matrix', KSMATRX
c         write (*,*) 'KTMATRX matrix', KTMATRX
c         write (*,*) 'Displacement values', U
c         pause
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
!This is for a single node having eight surrounding nodes including the node itself
!UEL actually calculates nodal stiffness; not the element stiffness
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
		parameter(NEG=-1.D0,zero=0.D0, half=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
		INTEGER :: numnode,numelem,Ncen
		REAL(8) :: node(numnode,2)
		INTEGER :: element(numelem,3),Nodnum(NDOFEL*half) 
C
		REAL(8) NN(6),NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2)
        REAL(8) N3(2),NArea,A11,B1(3,6),UNL(NDOFEL),detCC
		REAL(8) Bl(3,NDOFEL),pb(3,3),S(4,4),GT(3,NDOFEL),mu1
        REAL(8) lam,mu,F(3,3),CC(3,3),INVCC(3,3),dd1(3,3,3,3)
        REAL(8) CMAT(3,3),Fbar(3,4),BB(NDOFEL,4),pbs(3,1)
        REAL(8) GC(NDOFEL,3),GCGT(NDOFEL,NDOFEL),BBS(NDOFEL,4)
        REAL(8) BBSBT(NDOFEL,NDOFEL),GB(NDOFEL,1)
        REAL(8) KSMATRX(NDOFEL,1),KTMATRX(NDOFEL,NDOFEL) 
C
        INTEGER i,j,kk1,kk2,i2,nn1,nn2,nn3,x1,j1,k1,l1
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
		AMATRX=zero
		Bl=zero
        KSMATRX=zero
        KTMATRX=zero
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
    	do I=1,NDOFEL
         UNL(I) = U(I)
         RHS(I,1)=ZERO
        enddo
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,(NDOFEL*half)
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=ZERO
C---------- Assemble Constitutive matrix---D matrix-----C-------C-------C-------C-------
C---------- Pass material properties into variables-----C-------C-------C-------C-------
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition------C-------C-------C-------C-------C-------C-------
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data
! Loop for all elements		
	do ie = 1,numelem
	  A11=ZERO
	  N1=ZERO
	  N2=ZERO
	  N3=ZERO
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
	  call get_NArea_NBmat(N1,N2,N3,A11,B1)
	  do i=1,3
		    do j = 1,(NDOFEL*half)
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			  Bl(ii,(2*(j-1)+k)) = Bl(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
	  enddo
		NArea=NArea+A11/3
	  else
	  endif
	enddo
C	
	do i = 1,3
		do j = 1,NDOFEL
			Bl(i,j) = (1/NArea)*Bl(i,j)
		enddo
	enddo
c
c Calculation of initial stiffness and tangent stifness matrices
c
        pb=zero
		S=zero
        Fbar=zero
        BB=zero
		GT=zero
		GC=zero
        GCGT=zero
        BBS=zero
        BBSBT=zero
        GB=zero
        pbs=zero
        lam=zero
        mu=zero
        mu1=zero
        dd1=zero
        CMAT=zero
        F=zero
        detF=zero
        CC=zero
        INVCC=zero
        detCC=zero
        dudx = zero
        dudy = zero
        dvdx = zero
        dvdy = zero        
c
		do j = 1,(NDOFEL*half)
			dudx = dudx+Bl(1,(2*j-1))*UNL(2*j-1)
            dudy = dudy+Bl(2,2*j)*UNL(2*j-1)
            dvdx = dvdx+Bl(1,(2*j-1))*UNL(2*j)
            dvdy = dvdy+Bl(2,2*j)*UNL(2*j)
		enddo
c  
		lam=(EMU*E)/((1+EMU)*(1-2*EMU))	
        mu=E/(2*(1+EMU))	

		F(1,1)=dudx+1
        F(1,2)=dudy
        F(2,1)=dvdx
        F(2,2)=dvdy+1
        F(3,3)=1

        CC(1,1)=F(1,1)*F(1,1)+F(2,1)*F(2,1)
        CC(1,2)=F(1,2)*F(1,1)+F(2,1)*F(2,2)
        CC(2,1)=CC(1,2)
        CC(2,2)=F(1,2)*F(1,2)+F(2,2)*F(2,2)
        CC(3,3)=1

        detCC=(CC(1,1)*CC(2,2))-(CC(1,2)*CC(2,1))

        INVCC(1,1)=CC(2,2)/detCC
        INVCC(1,2)=-CC(1,2)/detCC
        INVCC(2,1)=-CC(2,1)/detCC
        INVCC(2,2)=CC(1,1)/detCC
        INVCC(3,3)=1
        
		mu1=mu-lam*half*log(detCC)

	   	DO ii=1,3
	  	DO jj=1,3
        	if (ii .EQ. jj)  then
	   		pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*(INVCC(ii,jj)-1)
        	else
            pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*INVCC(ii,jj)
          	endif
	  	ENDDO
        ENDDO
        
C Define CMAT matrix for Neohookian material model  
		
		DO x1=1,3
	  	DO j1=1,3
		DO k1=1,3
	  	DO l1=1,3
        dd1(x1,j1,k1,l1)=lam*INVCC(x1,j1)*INVCC(k1,l1)+
     1mu1*(INVCC(x1,k1)*INVCC(j1,l1)+INVCC(x1,l1)*INVCC(j1,k1))
	  	ENDDO
        ENDDO
        ENDDO
        ENDDO
		
	  	CMAT(1,1)=dd1(1,1,1,1)
	  	CMAT(1,2)=dd1(1,1,2,2)
	  	CMAT(1,3)=dd1(1,1,1,2)
	  	CMAT(2,1)=dd1(2,2,1,1)
      	CMAT(2,2)=dd1(2,2,2,2)
	  	CMAT(2,3)=dd1(2,2,1,2)  
      	CMAT(3,1)=dd1(1,2,1,1)
      	CMAT(3,2)=dd1(1,2,2,2)
	  	CMAT(3,3)=dd1(1,2,1,2)
c        
		S(1,1)=pb(1,1)
 		S(1,2)=pb(1,2)
        S(2,1)=pb(2,1)
        S(2,2)=pb(2,2)
        S(3,3)=pb(1,1)
        S(3,4)=pb(1,2)
        S(4,3)=pb(2,1)
        S(4,4)=pb(2,2)

        Fbar(1,1)=F(1,1)
 		Fbar(1,3)=F(2,1)
        Fbar(2,2)=F(1,2)
        Fbar(2,4)=F(2,2)
        Fbar(3,1)=F(1,2)
        Fbar(3,2)=F(1,1)
        Fbar(3,3)=F(2,2)
        Fbar(3,4)=F(2,1)
c
        DO i=1,(NDOFEL*half)
		BB((2*i-1),1)=Bl(1,(2*i-1))
        BB((2*i-1),2)=Bl(2,(2*i))
        BB((2*i),3)=Bl(1,(2*i-1))
        BB((2*i),4)=Bl(2,(2*i))
        ENDDO
c		
		DO i=1,3
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	GT(i,j)=GT(i,j)+Fbar(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
        pbs(1,1)=pb(1,1)
        pbs(2,1)=pb(2,2)
        pbs(3,1)=pb(1,2)
c
		DO i=1,NDOFEL
	  	DO j=1,3
        DO kk1=1,3
	   	GC(i,j)=GC(i,j)+GT(kk1,i)*CMAT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,3
	   	GCGT(i,j)=GCGT(i,j)+GC(i,kk1)*GT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c
	  	DO i=1,NDOFEL
        DO j=1,4
        DO kk1=1,4
	   	BBS(i,j)=BBS(i,j)+BB(i,kk1)*S(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	BBSBT(i,j)=BBSBT(i,j)+BBS(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
		DO i=1,NDOFEL
        DO kk2=1,3
	   	GB(i,1)=GB(i,1)+GT(kk2,i)*pbs(kk2,1)
	  	ENDDO
        ENDDO
c
      	DO i=1,NDOFEL
        KSMATRX(i,1)=KSMATRX(i,1)+GB(i,1)
        DO j=1,NDOFEL
		KTMATRX(i,j)=KTMATRX(i,j)+GCGT(i,j)+BBSBT(i,j)
        ENDDO
		ENDDO
C===============================C
C        calculate  RHS         C
C===============================C
  	  DO K1=1, NDOFEL              !***!   
	    RHS(K1, 1)=RHS(K1, 1)-NArea*KSMATRX(K1,1)
	  ENDDO
c
  	  DO K1=1,NDOFEL              !***!  DEFINING K MATRIX 
	     DO K2=1,NDOFEL              
	     AMATRX(K1, K2)=AMATRX(K1, K2)+NArea*KTMATRX(K1,K2)
	     ENDDO
	  ENDDO       
C---------- Output the saved results           
c         write (*,*) 'Node number',JELEM
c         write (*,*) 'Bmat matrix', GT
c         write (*,*) 'Bgeo matrix transpose', BB 
c         write (*,*) 'F matrix', F
c		  write (*,*) 'dd1 matrix', dd1
c         write (*,*) 'pb matrix', pb
c         write (*,*) 'Cmat matrix', CMAT
c         write (*,*) 'KSMATRX matrix', KSMATRX
c         write (*,*) 'KTMATRX matrix', KTMATRX
c         write (*,*) 'Displacement values', U
c         pause
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
!This is for a single node having nine surrounding nodes including the node itself
!UEL actually calculates nodal stiffness; not the element stiffness
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
		parameter(NEG=-1.D0,zero=0.D0, half=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
		INTEGER :: numnode,numelem,Ncen
		REAL(8) :: node(numnode,2)
		INTEGER :: element(numelem,3),Nodnum(NDOFEL*half) 
C
		REAL(8) NN(6),NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2)
        REAL(8) N3(2),NArea,A11,B1(3,6),UNL(NDOFEL),detCC
		REAL(8) Bl(3,NDOFEL),pb(3,3),S(4,4),GT(3,NDOFEL),mu1
        REAL(8) lam,mu,F(3,3),CC(3,3),INVCC(3,3),dd1(3,3,3,3)
        REAL(8) CMAT(3,3),Fbar(3,4),BB(NDOFEL,4),pbs(3,1)
        REAL(8) GC(NDOFEL,3),GCGT(NDOFEL,NDOFEL),BBS(NDOFEL,4)
        REAL(8) BBSBT(NDOFEL,NDOFEL),GB(NDOFEL,1)
        REAL(8) KSMATRX(NDOFEL,1),KTMATRX(NDOFEL,NDOFEL) 
C
        INTEGER i,j,kk1,kk2,i2,nn1,nn2,nn3,x1,j1,k1,l1
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
		AMATRX=zero
		Bl=zero
        KSMATRX=zero
        KTMATRX=zero
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
    	do I=1,NDOFEL
         UNL(I) = U(I)
         RHS(I,1)=ZERO
        enddo
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,(NDOFEL*half)
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=ZERO
C---------- Assemble Constitutive matrix---D matrix-----C-------C-------C-------C-------
C---------- Pass material properties into variables-----C-------C-------C-------C-------
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition------C-------C-------C-------C-------C-------C-------
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data
! Loop for all elements		
	do ie = 1,numelem
	  A11=ZERO
	  N1=ZERO
	  N2=ZERO
	  N3=ZERO
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
	  call get_NArea_NBmat(N1,N2,N3,A11,B1)
	  do i=1,3
		    do j = 1,(NDOFEL*half)
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			  Bl(ii,(2*(j-1)+k)) = Bl(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
	  enddo
		NArea=NArea+A11/3
	  else
	  endif
	enddo
C	
	do i = 1,3
		do j = 1,NDOFEL
			Bl(i,j) = (1/NArea)*Bl(i,j)
		enddo
	enddo
c
c Calculation of initial stiffness and tangent stifness matrices
c
        pb=zero
		S=zero
        Fbar=zero
        BB=zero
		GT=zero
		GC=zero
        GCGT=zero
        BBS=zero
        BBSBT=zero
        GB=zero
        pbs=zero
        lam=zero
        mu=zero
        mu1=zero
        dd1=zero
        CMAT=zero
        F=zero
        detF=zero
        CC=zero
        INVCC=zero
        detCC=zero
        dudx = zero
        dudy = zero
        dvdx = zero
        dvdy = zero        
c
		do j = 1,(NDOFEL*half)
			dudx = dudx+Bl(1,(2*j-1))*UNL(2*j-1)
            dudy = dudy+Bl(2,2*j)*UNL(2*j-1)
            dvdx = dvdx+Bl(1,(2*j-1))*UNL(2*j)
            dvdy = dvdy+Bl(2,2*j)*UNL(2*j)
		enddo
c  
		lam=(EMU*E)/((1+EMU)*(1-2*EMU))	
        mu=E/(2*(1+EMU))	

		F(1,1)=dudx+1
        F(1,2)=dudy
        F(2,1)=dvdx
        F(2,2)=dvdy+1
        F(3,3)=1

        CC(1,1)=F(1,1)*F(1,1)+F(2,1)*F(2,1)
        CC(1,2)=F(1,2)*F(1,1)+F(2,1)*F(2,2)
        CC(2,1)=CC(1,2)
        CC(2,2)=F(1,2)*F(1,2)+F(2,2)*F(2,2)
        CC(3,3)=1

        detCC=(CC(1,1)*CC(2,2))-(CC(1,2)*CC(2,1))

        INVCC(1,1)=CC(2,2)/detCC
        INVCC(1,2)=-CC(1,2)/detCC
        INVCC(2,1)=-CC(2,1)/detCC
        INVCC(2,2)=CC(1,1)/detCC
        INVCC(3,3)=1
        
		mu1=mu-lam*half*log(detCC)

	   	DO ii=1,3
	  	DO jj=1,3
        	if (ii .EQ. jj)  then
	   		pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*(INVCC(ii,jj)-1)
        	else
            pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*INVCC(ii,jj)
          	endif
	  	ENDDO
        ENDDO
        
C Define CMAT matrix for Neohookian material model  
		
		DO x1=1,3
	  	DO j1=1,3
		DO k1=1,3
	  	DO l1=1,3
        dd1(x1,j1,k1,l1)=lam*INVCC(x1,j1)*INVCC(k1,l1)+
     1mu1*(INVCC(x1,k1)*INVCC(j1,l1)+INVCC(x1,l1)*INVCC(j1,k1))
	  	ENDDO
        ENDDO
        ENDDO
        ENDDO
		
	  	CMAT(1,1)=dd1(1,1,1,1)
	  	CMAT(1,2)=dd1(1,1,2,2)
	  	CMAT(1,3)=dd1(1,1,1,2)
	  	CMAT(2,1)=dd1(2,2,1,1)
      	CMAT(2,2)=dd1(2,2,2,2)
	  	CMAT(2,3)=dd1(2,2,1,2)  
      	CMAT(3,1)=dd1(1,2,1,1)
      	CMAT(3,2)=dd1(1,2,2,2)
	  	CMAT(3,3)=dd1(1,2,1,2)
c        
		S(1,1)=pb(1,1)
 		S(1,2)=pb(1,2)
        S(2,1)=pb(2,1)
        S(2,2)=pb(2,2)
        S(3,3)=pb(1,1)
        S(3,4)=pb(1,2)
        S(4,3)=pb(2,1)
        S(4,4)=pb(2,2)

        Fbar(1,1)=F(1,1)
 		Fbar(1,3)=F(2,1)
        Fbar(2,2)=F(1,2)
        Fbar(2,4)=F(2,2)
        Fbar(3,1)=F(1,2)
        Fbar(3,2)=F(1,1)
        Fbar(3,3)=F(2,2)
        Fbar(3,4)=F(2,1)
c
        DO i=1,(NDOFEL*half)
		BB((2*i-1),1)=Bl(1,(2*i-1))
        BB((2*i-1),2)=Bl(2,(2*i))
        BB((2*i),3)=Bl(1,(2*i-1))
        BB((2*i),4)=Bl(2,(2*i))
        ENDDO
c		
		DO i=1,3
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	GT(i,j)=GT(i,j)+Fbar(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
        pbs(1,1)=pb(1,1)
        pbs(2,1)=pb(2,2)
        pbs(3,1)=pb(1,2)
c
		DO i=1,NDOFEL
	  	DO j=1,3
        DO kk1=1,3
	   	GC(i,j)=GC(i,j)+GT(kk1,i)*CMAT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,3
	   	GCGT(i,j)=GCGT(i,j)+GC(i,kk1)*GT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c
	  	DO i=1,NDOFEL
        DO j=1,4
        DO kk1=1,4
	   	BBS(i,j)=BBS(i,j)+BB(i,kk1)*S(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	BBSBT(i,j)=BBSBT(i,j)+BBS(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
		DO i=1,NDOFEL
        DO kk2=1,3
	   	GB(i,1)=GB(i,1)+GT(kk2,i)*pbs(kk2,1)
	  	ENDDO
        ENDDO
c
      	DO i=1,NDOFEL
        KSMATRX(i,1)=KSMATRX(i,1)+GB(i,1)
        DO j=1,NDOFEL
		KTMATRX(i,j)=KTMATRX(i,j)+GCGT(i,j)+BBSBT(i,j)
        ENDDO
		ENDDO
C===============================C
C        calculate  RHS         C
C===============================C
  	  DO K1=1, NDOFEL              !***!   
	    RHS(K1, 1)=RHS(K1, 1)-NArea*KSMATRX(K1,1)
	  ENDDO
c
  	  DO K1=1,NDOFEL              !***!  DEFINING K MATRIX 
	     DO K2=1,NDOFEL              
	     AMATRX(K1, K2)=AMATRX(K1, K2)+NArea*KTMATRX(K1,K2)
	     ENDDO
	  ENDDO       
C---------- Output the saved results           
c         write (*,*) 'Node number',JELEM
c         write (*,*) 'Bmat matrix', GT
c         write (*,*) 'Bgeo matrix transpose', BB 
c         write (*,*) 'F matrix', F
c		  write (*,*) 'dd1 matrix', dd1
c         write (*,*) 'pb matrix', pb
c         write (*,*) 'Cmat matrix', CMAT
c         write (*,*) 'KSMATRX matrix', KSMATRX
c         write (*,*) 'KTMATRX matrix', KTMATRX
c         write (*,*) 'Displacement values', U
c         pause
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
!This is for a single node having ten surrounding nodes including the node itself
!UEL actually calculates nodal stiffness; not the element stiffness
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
		parameter(NEG=-1.D0,zero=0.D0, half=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
		INTEGER :: numnode,numelem,Ncen
		REAL(8) :: node(numnode,2)
		INTEGER :: element(numelem,3),Nodnum(NDOFEL*half) 
C
		REAL(8) NN(6),NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2)
        REAL(8) N3(2),NArea,A11,B1(3,6),UNL(NDOFEL),detCC
		REAL(8) Bl(3,NDOFEL),pb(3,3),S(4,4),GT(3,NDOFEL),mu1
        REAL(8) lam,mu,F(3,3),CC(3,3),INVCC(3,3),dd1(3,3,3,3)
        REAL(8) CMAT(3,3),Fbar(3,4),BB(NDOFEL,4),pbs(3,1)
        REAL(8) GC(NDOFEL,3),GCGT(NDOFEL,NDOFEL),BBS(NDOFEL,4)
        REAL(8) BBSBT(NDOFEL,NDOFEL),GB(NDOFEL,1)
        REAL(8) KSMATRX(NDOFEL,1),KTMATRX(NDOFEL,NDOFEL) 
C
        INTEGER i,j,kk1,kk2,i2,nn1,nn2,nn3,x1,j1,k1,l1
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
		AMATRX=zero
		Bl=zero
        KSMATRX=zero
        KTMATRX=zero
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
    	do I=1,NDOFEL
         UNL(I) = U(I)
         RHS(I,1)=ZERO
        enddo
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,(NDOFEL*half)
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=ZERO
C---------- Assemble Constitutive matrix---D matrix-----C-------C-------C-------C-------
C---------- Pass material properties into variables-----C-------C-------C-------C-------
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition------C-------C-------C-------C-------C-------C-------
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data
! Loop for all elements		
	do ie = 1,numelem
	  A11=ZERO
	  N1=ZERO
	  N2=ZERO
	  N3=ZERO
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
	  call get_NArea_NBmat(N1,N2,N3,A11,B1)
	  do i=1,3
		    do j = 1,(NDOFEL*half)
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			  Bl(ii,(2*(j-1)+k)) = Bl(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
	  enddo
		NArea=NArea+A11/3
	  else
	  endif
	enddo
C	
	do i = 1,3
		do j = 1,NDOFEL
			Bl(i,j) = (1/NArea)*Bl(i,j)
		enddo
	enddo
c
c Calculation of initial stiffness and tangent stifness matrices
c
        pb=zero
		S=zero
        Fbar=zero
        BB=zero
		GT=zero
		GC=zero
        GCGT=zero
        BBS=zero
        BBSBT=zero
        GB=zero
        pbs=zero
        lam=zero
        mu=zero
        mu1=zero
        dd1=zero
        CMAT=zero
        F=zero
        detF=zero
        CC=zero
        INVCC=zero
        detCC=zero
        dudx = zero
        dudy = zero
        dvdx = zero
        dvdy = zero        
c
		do j = 1,(NDOFEL*half)
			dudx = dudx+Bl(1,(2*j-1))*UNL(2*j-1)
            dudy = dudy+Bl(2,2*j)*UNL(2*j-1)
            dvdx = dvdx+Bl(1,(2*j-1))*UNL(2*j)
            dvdy = dvdy+Bl(2,2*j)*UNL(2*j)
		enddo
c  
		lam=(EMU*E)/((1+EMU)*(1-2*EMU))	
        mu=E/(2*(1+EMU))	

		F(1,1)=dudx+1
        F(1,2)=dudy
        F(2,1)=dvdx
        F(2,2)=dvdy+1
        F(3,3)=1

        CC(1,1)=F(1,1)*F(1,1)+F(2,1)*F(2,1)
        CC(1,2)=F(1,2)*F(1,1)+F(2,1)*F(2,2)
        CC(2,1)=CC(1,2)
        CC(2,2)=F(1,2)*F(1,2)+F(2,2)*F(2,2)
        CC(3,3)=1

        detCC=(CC(1,1)*CC(2,2))-(CC(1,2)*CC(2,1))

        INVCC(1,1)=CC(2,2)/detCC
        INVCC(1,2)=-CC(1,2)/detCC
        INVCC(2,1)=-CC(2,1)/detCC
        INVCC(2,2)=CC(1,1)/detCC
        INVCC(3,3)=1
        
		mu1=mu-lam*half*log(detCC)

	   	DO ii=1,3
	  	DO jj=1,3
        	if (ii .EQ. jj)  then
	   		pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*(INVCC(ii,jj)-1)
        	else
            pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*INVCC(ii,jj)
          	endif
	  	ENDDO
        ENDDO
        
C Define CMAT matrix for Neohookian material model  
		
		DO x1=1,3
	  	DO j1=1,3
		DO k1=1,3
	  	DO l1=1,3
        dd1(x1,j1,k1,l1)=lam*INVCC(x1,j1)*INVCC(k1,l1)+
     1mu1*(INVCC(x1,k1)*INVCC(j1,l1)+INVCC(x1,l1)*INVCC(j1,k1))
	  	ENDDO
        ENDDO
        ENDDO
        ENDDO
		
	  	CMAT(1,1)=dd1(1,1,1,1)
	  	CMAT(1,2)=dd1(1,1,2,2)
	  	CMAT(1,3)=dd1(1,1,1,2)
	  	CMAT(2,1)=dd1(2,2,1,1)
      	CMAT(2,2)=dd1(2,2,2,2)
	  	CMAT(2,3)=dd1(2,2,1,2)  
      	CMAT(3,1)=dd1(1,2,1,1)
      	CMAT(3,2)=dd1(1,2,2,2)
	  	CMAT(3,3)=dd1(1,2,1,2)
c        
		S(1,1)=pb(1,1)
 		S(1,2)=pb(1,2)
        S(2,1)=pb(2,1)
        S(2,2)=pb(2,2)
        S(3,3)=pb(1,1)
        S(3,4)=pb(1,2)
        S(4,3)=pb(2,1)
        S(4,4)=pb(2,2)

        Fbar(1,1)=F(1,1)
 		Fbar(1,3)=F(2,1)
        Fbar(2,2)=F(1,2)
        Fbar(2,4)=F(2,2)
        Fbar(3,1)=F(1,2)
        Fbar(3,2)=F(1,1)
        Fbar(3,3)=F(2,2)
        Fbar(3,4)=F(2,1)
c
        DO i=1,(NDOFEL*half)
		BB((2*i-1),1)=Bl(1,(2*i-1))
        BB((2*i-1),2)=Bl(2,(2*i))
        BB((2*i),3)=Bl(1,(2*i-1))
        BB((2*i),4)=Bl(2,(2*i))
        ENDDO
c		
		DO i=1,3
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	GT(i,j)=GT(i,j)+Fbar(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
        pbs(1,1)=pb(1,1)
        pbs(2,1)=pb(2,2)
        pbs(3,1)=pb(1,2)
c
		DO i=1,NDOFEL
	  	DO j=1,3
        DO kk1=1,3
	   	GC(i,j)=GC(i,j)+GT(kk1,i)*CMAT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,3
	   	GCGT(i,j)=GCGT(i,j)+GC(i,kk1)*GT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c
	  	DO i=1,NDOFEL
        DO j=1,4
        DO kk1=1,4
	   	BBS(i,j)=BBS(i,j)+BB(i,kk1)*S(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	BBSBT(i,j)=BBSBT(i,j)+BBS(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
		DO i=1,NDOFEL
        DO kk2=1,3
	   	GB(i,1)=GB(i,1)+GT(kk2,i)*pbs(kk2,1)
	  	ENDDO
        ENDDO
c
      	DO i=1,NDOFEL
        KSMATRX(i,1)=KSMATRX(i,1)+GB(i,1)
        DO j=1,NDOFEL
		KTMATRX(i,j)=KTMATRX(i,j)+GCGT(i,j)+BBSBT(i,j)
        ENDDO
		ENDDO
C===============================C
C        calculate  RHS         C
C===============================C
  	  DO K1=1, NDOFEL              !***!   
	    RHS(K1, 1)=RHS(K1, 1)-NArea*KSMATRX(K1,1)
	  ENDDO
c
  	  DO K1=1,NDOFEL              !***!  DEFINING K MATRIX 
	     DO K2=1,NDOFEL              
	     AMATRX(K1, K2)=AMATRX(K1, K2)+NArea*KTMATRX(K1,K2)
	     ENDDO
	  ENDDO       
C---------- Output the saved results           
c         write (*,*) 'Node number',JELEM
c         write (*,*) 'Bmat matrix', GT
c         write (*,*) 'Bgeo matrix transpose', BB 
c         write (*,*) 'F matrix', F
c		  write (*,*) 'dd1 matrix', dd1
c         write (*,*) 'pb matrix', pb
c         write (*,*) 'Cmat matrix', CMAT
c         write (*,*) 'KSMATRX matrix', KSMATRX
c         write (*,*) 'KTMATRX matrix', KTMATRX
c         write (*,*) 'Displacement values', U
c         pause
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
!This is for a single node having eleven surrounding nodes including the node itself
!UEL actually calculates nodal stiffness; not the element stiffness
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
		parameter(NEG=-1.D0,zero=0.D0, half=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
		INTEGER :: numnode,numelem,Ncen
		REAL(8) :: node(numnode,2)
		INTEGER :: element(numelem,3),Nodnum(NDOFEL*half)  
C
		REAL(8) NN(6),NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2)
        REAL(8) N3(2),NArea,A11,B1(3,6),UNL(NDOFEL),detCC
		REAL(8) Bl(3,NDOFEL),pb(3,3),S(4,4),GT(3,NDOFEL),mu1
        REAL(8) lam,mu,F(3,3),CC(3,3),INVCC(3,3),dd1(3,3,3,3)
        REAL(8) CMAT(3,3),Fbar(3,4),BB(NDOFEL,4),pbs(3,1)
        REAL(8) GC(NDOFEL,3),GCGT(NDOFEL,NDOFEL),BBS(NDOFEL,4)
        REAL(8) BBSBT(NDOFEL,NDOFEL),GB(NDOFEL,1)
        REAL(8) KSMATRX(NDOFEL,1),KTMATRX(NDOFEL,NDOFEL) 
C
        INTEGER i,j,kk1,kk2,i2,nn1,nn2,nn3,x1,j1,k1,l1
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
		AMATRX=zero
		Bl=zero
        KSMATRX=zero
        KTMATRX=zero
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
    	do I=1,NDOFEL
         UNL(I) = U(I)
         RHS(I,1)=ZERO
        enddo
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,(NDOFEL*half)
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=ZERO
C---------- Assemble Constitutive matrix---D matrix-----C-------C-------C-------C-------
C---------- Pass material properties into variables-----C-------C-------C-------C-------
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition------C-------C-------C-------C-------C-------C-------
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data
! Loop for all elements		
	do ie = 1,numelem
	  A11=ZERO
	  N1=ZERO
	  N2=ZERO
	  N3=ZERO
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
	  call get_NArea_NBmat(N1,N2,N3,A11,B1)
	  do i=1,3
		    do j = 1,(NDOFEL*half)
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			  Bl(ii,(2*(j-1)+k)) = Bl(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
	  enddo
		NArea=NArea+A11/3
	  else
	  endif
	enddo
C	
	do i = 1,3
		do j = 1,NDOFEL
			Bl(i,j) = (1/NArea)*Bl(i,j)
		enddo
	enddo
c
c Calculation of initial stiffness and tangent stifness matrices
c
        pb=zero
		S=zero
        Fbar=zero
        BB=zero
		GT=zero
		GC=zero
        GCGT=zero
        BBS=zero
        BBSBT=zero
        GB=zero
        pbs=zero
        lam=zero
        mu=zero
        mu1=zero
        dd1=zero
        CMAT=zero
        F=zero
        detF=zero
        CC=zero
        INVCC=zero
        detCC=zero
        dudx = zero
        dudy = zero
        dvdx = zero
        dvdy = zero        
c
		do j = 1,(NDOFEL*half)
			dudx = dudx+Bl(1,(2*j-1))*UNL(2*j-1)
            dudy = dudy+Bl(2,2*j)*UNL(2*j-1)
            dvdx = dvdx+Bl(1,(2*j-1))*UNL(2*j)
            dvdy = dvdy+Bl(2,2*j)*UNL(2*j)
		enddo
c  
		lam=(EMU*E)/((1+EMU)*(1-2*EMU))	
        mu=E/(2*(1+EMU))	

		F(1,1)=dudx+1
        F(1,2)=dudy
        F(2,1)=dvdx
        F(2,2)=dvdy+1
        F(3,3)=1

        CC(1,1)=F(1,1)*F(1,1)+F(2,1)*F(2,1)
        CC(1,2)=F(1,2)*F(1,1)+F(2,1)*F(2,2)
        CC(2,1)=CC(1,2)
        CC(2,2)=F(1,2)*F(1,2)+F(2,2)*F(2,2)
        CC(3,3)=1

        detCC=(CC(1,1)*CC(2,2))-(CC(1,2)*CC(2,1))

        INVCC(1,1)=CC(2,2)/detCC
        INVCC(1,2)=-CC(1,2)/detCC
        INVCC(2,1)=-CC(2,1)/detCC
        INVCC(2,2)=CC(1,1)/detCC
        INVCC(3,3)=1
        
		mu1=mu-lam*half*log(detCC)

	   	DO ii=1,3
	  	DO jj=1,3
        	if (ii .EQ. jj)  then
	   		pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*(INVCC(ii,jj)-1)
        	else
            pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*INVCC(ii,jj)
          	endif
	  	ENDDO
        ENDDO
        
C Define CMAT matrix for Neohookian material model  
		
		DO x1=1,3
	  	DO j1=1,3
		DO k1=1,3
	  	DO l1=1,3
        dd1(x1,j1,k1,l1)=lam*INVCC(x1,j1)*INVCC(k1,l1)+
     1mu1*(INVCC(x1,k1)*INVCC(j1,l1)+INVCC(x1,l1)*INVCC(j1,k1))
	  	ENDDO
        ENDDO
        ENDDO
        ENDDO
		
	  	CMAT(1,1)=dd1(1,1,1,1)
	  	CMAT(1,2)=dd1(1,1,2,2)
	  	CMAT(1,3)=dd1(1,1,1,2)
	  	CMAT(2,1)=dd1(2,2,1,1)
      	CMAT(2,2)=dd1(2,2,2,2)
	  	CMAT(2,3)=dd1(2,2,1,2)  
      	CMAT(3,1)=dd1(1,2,1,1)
      	CMAT(3,2)=dd1(1,2,2,2)
	  	CMAT(3,3)=dd1(1,2,1,2)
c        
		S(1,1)=pb(1,1)
 		S(1,2)=pb(1,2)
        S(2,1)=pb(2,1)
        S(2,2)=pb(2,2)
        S(3,3)=pb(1,1)
        S(3,4)=pb(1,2)
        S(4,3)=pb(2,1)
        S(4,4)=pb(2,2)

        Fbar(1,1)=F(1,1)
 		Fbar(1,3)=F(2,1)
        Fbar(2,2)=F(1,2)
        Fbar(2,4)=F(2,2)
        Fbar(3,1)=F(1,2)
        Fbar(3,2)=F(1,1)
        Fbar(3,3)=F(2,2)
        Fbar(3,4)=F(2,1)
c
        DO i=1,(NDOFEL*half)
		BB((2*i-1),1)=Bl(1,(2*i-1))
        BB((2*i-1),2)=Bl(2,(2*i))
        BB((2*i),3)=Bl(1,(2*i-1))
        BB((2*i),4)=Bl(2,(2*i))
        ENDDO
c		
		DO i=1,3
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	GT(i,j)=GT(i,j)+Fbar(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
        pbs(1,1)=pb(1,1)
        pbs(2,1)=pb(2,2)
        pbs(3,1)=pb(1,2)
c
		DO i=1,NDOFEL
	  	DO j=1,3
        DO kk1=1,3
	   	GC(i,j)=GC(i,j)+GT(kk1,i)*CMAT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,3
	   	GCGT(i,j)=GCGT(i,j)+GC(i,kk1)*GT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c
	  	DO i=1,NDOFEL
        DO j=1,4
        DO kk1=1,4
	   	BBS(i,j)=BBS(i,j)+BB(i,kk1)*S(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	BBSBT(i,j)=BBSBT(i,j)+BBS(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
		DO i=1,NDOFEL
        DO kk2=1,3
	   	GB(i,1)=GB(i,1)+GT(kk2,i)*pbs(kk2,1)
	  	ENDDO
        ENDDO
c
      	DO i=1,NDOFEL
        KSMATRX(i,1)=KSMATRX(i,1)+GB(i,1)
        DO j=1,NDOFEL
		KTMATRX(i,j)=KTMATRX(i,j)+GCGT(i,j)+BBSBT(i,j)
        ENDDO
		ENDDO
C===============================C
C        calculate  RHS         C
C===============================C
  	  DO K1=1, NDOFEL              !***!   
	    RHS(K1, 1)=RHS(K1, 1)-NArea*KSMATRX(K1,1)
	  ENDDO
c
  	  DO K1=1,NDOFEL              !***!  DEFINING K MATRIX 
	     DO K2=1,NDOFEL              
	     AMATRX(K1, K2)=AMATRX(K1, K2)+NArea*KTMATRX(K1,K2)
	     ENDDO
	  ENDDO       
C---------- Output the saved results           
c         write (*,*) 'Node number',JELEM
c         write (*,*) 'Bmat matrix', GT
c         write (*,*) 'Bgeo matrix transpose', BB 
c         write (*,*) 'F matrix', F
c		  write (*,*) 'dd1 matrix', dd1
c         write (*,*) 'pb matrix', pb
c         write (*,*) 'Cmat matrix', CMAT
c         write (*,*) 'KSMATRX matrix', KSMATRX
c         write (*,*) 'KTMATRX matrix', KTMATRX
c         write (*,*) 'Displacement values', U
c         pause
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
!This is for a single node having twelve surrounding nodes including the node itself
!UEL actually calculates nodal stiffness; not the element stiffness
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
		parameter(NEG=-1.D0,zero=0.D0, half=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
		INTEGER :: numnode,numelem,Ncen
		REAL(8) :: node(numnode,2)
		INTEGER :: element(numelem,3),Nodnum(NDOFEL*half)  
C
		REAL(8) NN(6),NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2)
        REAL(8) N3(2),NArea,A11,B1(3,6),UNL(NDOFEL),detCC
		REAL(8) Bl(3,NDOFEL),pb(3,3),S(4,4),GT(3,NDOFEL),mu1
        REAL(8) lam,mu,F(3,3),CC(3,3),INVCC(3,3),dd1(3,3,3,3)
        REAL(8) CMAT(3,3),Fbar(3,4),BB(NDOFEL,4),pbs(3,1)
        REAL(8) GC(NDOFEL,3),GCGT(NDOFEL,NDOFEL),BBS(NDOFEL,4)
        REAL(8) BBSBT(NDOFEL,NDOFEL),GB(NDOFEL,1)
        REAL(8) KSMATRX(NDOFEL,1),KTMATRX(NDOFEL,NDOFEL) 
C
        INTEGER i,j,kk1,kk2,i2,nn1,nn2,nn3,x1,j1,k1,l1
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
		AMATRX=zero
		Bl=zero
        KSMATRX=zero
        KTMATRX=zero
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
    	do I=1,NDOFEL
         UNL(I) = U(I)
         RHS(I,1)=ZERO
        enddo
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,(NDOFEL*half)
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=ZERO
C---------- Assemble Constitutive matrix---D matrix-----C-------C-------C-------C-------
C---------- Pass material properties into variables-----C-------C-------C-------C-------
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition------C-------C-------C-------C-------C-------C-------
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data
! Loop for all elements		
	do ie = 1,numelem
	  A11=ZERO
	  N1=ZERO
	  N2=ZERO
	  N3=ZERO
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
	  call get_NArea_NBmat(N1,N2,N3,A11,B1)
	  do i=1,3
		    do j = 1,(NDOFEL*half)
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			  Bl(ii,(2*(j-1)+k)) = Bl(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
	  enddo
		NArea=NArea+A11/3
	  else
	  endif
	enddo
C	
	do i = 1,3
		do j = 1,NDOFEL
			Bl(i,j) = (1/NArea)*Bl(i,j)
		enddo
	enddo
c
c Calculation of initial stiffness and tangent stifness matrices
c
        pb=zero
		S=zero
        Fbar=zero
        BB=zero
		GT=zero
		GC=zero
        GCGT=zero
        BBS=zero
        BBSBT=zero
        GB=zero
        pbs=zero
        lam=zero
        mu=zero
        mu1=zero
        dd1=zero
        CMAT=zero
        F=zero
        detF=zero
        CC=zero
        INVCC=zero
        detCC=zero
        dudx = zero
        dudy = zero
        dvdx = zero
        dvdy = zero        
c
		do j = 1,(NDOFEL*half)
			dudx = dudx+Bl(1,(2*j-1))*UNL(2*j-1)
            dudy = dudy+Bl(2,2*j)*UNL(2*j-1)
            dvdx = dvdx+Bl(1,(2*j-1))*UNL(2*j)
            dvdy = dvdy+Bl(2,2*j)*UNL(2*j)
		enddo
c  
		lam=(EMU*E)/((1+EMU)*(1-2*EMU))	
        mu=E/(2*(1+EMU))	

		F(1,1)=dudx+1
        F(1,2)=dudy
        F(2,1)=dvdx
        F(2,2)=dvdy+1
        F(3,3)=1

        CC(1,1)=F(1,1)*F(1,1)+F(2,1)*F(2,1)
        CC(1,2)=F(1,2)*F(1,1)+F(2,1)*F(2,2)
        CC(2,1)=CC(1,2)
        CC(2,2)=F(1,2)*F(1,2)+F(2,2)*F(2,2)
        CC(3,3)=1

        detCC=(CC(1,1)*CC(2,2))-(CC(1,2)*CC(2,1))

        INVCC(1,1)=CC(2,2)/detCC
        INVCC(1,2)=-CC(1,2)/detCC
        INVCC(2,1)=-CC(2,1)/detCC
        INVCC(2,2)=CC(1,1)/detCC
        INVCC(3,3)=1
        
		mu1=mu-lam*half*log(detCC)

	   	DO ii=1,3
	  	DO jj=1,3
        	if (ii .EQ. jj)  then
	   		pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*(INVCC(ii,jj)-1)
        	else
            pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*INVCC(ii,jj)
          	endif
	  	ENDDO
        ENDDO
        
C Define CMAT matrix for Neohookian material model  
		
		DO x1=1,3
	  	DO j1=1,3
		DO k1=1,3
	  	DO l1=1,3
        dd1(x1,j1,k1,l1)=lam*INVCC(x1,j1)*INVCC(k1,l1)+
     1mu1*(INVCC(x1,k1)*INVCC(j1,l1)+INVCC(x1,l1)*INVCC(j1,k1))
	  	ENDDO
        ENDDO
        ENDDO
        ENDDO
		
	  	CMAT(1,1)=dd1(1,1,1,1)
	  	CMAT(1,2)=dd1(1,1,2,2)
	  	CMAT(1,3)=dd1(1,1,1,2)
	  	CMAT(2,1)=dd1(2,2,1,1)
      	CMAT(2,2)=dd1(2,2,2,2)
	  	CMAT(2,3)=dd1(2,2,1,2)  
      	CMAT(3,1)=dd1(1,2,1,1)
      	CMAT(3,2)=dd1(1,2,2,2)
	  	CMAT(3,3)=dd1(1,2,1,2)
c        
		S(1,1)=pb(1,1)
 		S(1,2)=pb(1,2)
        S(2,1)=pb(2,1)
        S(2,2)=pb(2,2)
        S(3,3)=pb(1,1)
        S(3,4)=pb(1,2)
        S(4,3)=pb(2,1)
        S(4,4)=pb(2,2)

        Fbar(1,1)=F(1,1)
 		Fbar(1,3)=F(2,1)
        Fbar(2,2)=F(1,2)
        Fbar(2,4)=F(2,2)
        Fbar(3,1)=F(1,2)
        Fbar(3,2)=F(1,1)
        Fbar(3,3)=F(2,2)
        Fbar(3,4)=F(2,1)
c
        DO i=1,(NDOFEL*half)
		BB((2*i-1),1)=Bl(1,(2*i-1))
        BB((2*i-1),2)=Bl(2,(2*i))
        BB((2*i),3)=Bl(1,(2*i-1))
        BB((2*i),4)=Bl(2,(2*i))
        ENDDO
c		
		DO i=1,3
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	GT(i,j)=GT(i,j)+Fbar(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
        pbs(1,1)=pb(1,1)
        pbs(2,1)=pb(2,2)
        pbs(3,1)=pb(1,2)
c
		DO i=1,NDOFEL
	  	DO j=1,3
        DO kk1=1,3
	   	GC(i,j)=GC(i,j)+GT(kk1,i)*CMAT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,3
	   	GCGT(i,j)=GCGT(i,j)+GC(i,kk1)*GT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c
	  	DO i=1,NDOFEL
        DO j=1,4
        DO kk1=1,4
	   	BBS(i,j)=BBS(i,j)+BB(i,kk1)*S(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	BBSBT(i,j)=BBSBT(i,j)+BBS(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
		DO i=1,NDOFEL
        DO kk2=1,3
	   	GB(i,1)=GB(i,1)+GT(kk2,i)*pbs(kk2,1)
	  	ENDDO
        ENDDO
c
      	DO i=1,NDOFEL
        KSMATRX(i,1)=KSMATRX(i,1)+GB(i,1)
        DO j=1,NDOFEL
		KTMATRX(i,j)=KTMATRX(i,j)+GCGT(i,j)+BBSBT(i,j)
        ENDDO
		ENDDO
C===============================C
C        calculate  RHS         C
C===============================C
  	  DO K1=1, NDOFEL              !***!   
	    RHS(K1, 1)=RHS(K1, 1)-NArea*KSMATRX(K1,1)
	  ENDDO
c
  	  DO K1=1,NDOFEL              !***!  DEFINING K MATRIX 
	     DO K2=1,NDOFEL              
	     AMATRX(K1, K2)=AMATRX(K1, K2)+NArea*KTMATRX(K1,K2)
	     ENDDO
	  ENDDO       
C---------- Output the saved results           
c         write (*,*) 'Node number',JELEM
c         write (*,*) 'Bmat matrix', GT
c         write (*,*) 'Bgeo matrix transpose', BB 
c         write (*,*) 'F matrix', F
c		  write (*,*) 'dd1 matrix', dd1
c         write (*,*) 'pb matrix', pb
c         write (*,*) 'Cmat matrix', CMAT
c         write (*,*) 'KSMATRX matrix', KSMATRX
c         write (*,*) 'KTMATRX matrix', KTMATRX
c         write (*,*) 'Displacement values', U
c         pause
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
!This is for a single node having thirteen surrounding nodes including the node itself
!UEL actually calculates nodal stiffness; not the element stiffness
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
		parameter(NEG=-1.D0,zero=0.D0, half=0.5D0,ONE=1.D0,TWO=2.D0,SIX=6.D0)
		INTEGER :: numnode,numelem,Ncen
		REAL(8) :: node(numnode,2)
		INTEGER :: element(numelem,3),Nodnum(NDOFEL*half)  
C
		REAL(8) NN(6),NSTRAIN(3,1),NSTRESS(3,1),N1(2),N2(2)
        REAL(8) N3(2),NArea,A11,B1(3,6),UNL(NDOFEL),detCC
		REAL(8) Bl(3,NDOFEL),pb(3,3),S(4,4),GT(3,NDOFEL),mu1
        REAL(8) lam,mu,F(3,3),CC(3,3),INVCC(3,3),dd1(3,3,3,3)
        REAL(8) CMAT(3,3),Fbar(3,4),BB(NDOFEL,4),pbs(3,1)
        REAL(8) GC(NDOFEL,3),GCGT(NDOFEL,NDOFEL),BBS(NDOFEL,4)
        REAL(8) BBSBT(NDOFEL,NDOFEL),GB(NDOFEL,1)
        REAL(8) KSMATRX(NDOFEL,1),KTMATRX(NDOFEL,NDOFEL) 
C
        INTEGER i,j,kk1,kk2,i2,nn1,nn2,nn3,x1,j1,k1,l1
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
C---------- Initialize matrices
		AMATRX=zero
		Bl=zero
        KSMATRX=zero
        KTMATRX=zero
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
    	do I=1,NDOFEL
         UNL(I) = U(I)
         RHS(I,1)=ZERO
        enddo
C---------- Define coordinates for N1-N3--------C-------C-------C-------C-------C-------
! Get the node numbers for surrounding nodes	
	Ncen=JELEM
!	write (*,*) Ncen
	do i = 1,(NDOFEL*half)
		do j = 1,numnode
		if (COORDS(1,i) .EQ. node(j,1) .AND. COORDS(2,i) .EQ. node(j,2)) then
			Nodnum(i)=j
		endif
		enddo
	enddo	
	NArea=ZERO
C---------- Assemble Constitutive matrix---D matrix-----C-------C-------C-------C-------
C---------- Pass material properties into variables-----C-------C-------C-------C-------
	THICK=PROPS(1)
	E=PROPS(2)
	EMU=PROPS(3)
C---------- Plane Stress condition------C-------C-------C-------C-------C-------C-------
C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------C-------
! Construct B matrix using element order data
! Loop for all elements		
	do ie = 1,numelem
	  A11=ZERO
	  N1=ZERO
	  N2=ZERO
	  N3=ZERO
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
	  call get_NArea_NBmat(N1,N2,N3,A11,B1)
	  do i=1,3
		    do j = 1,(NDOFEL*half)
		    if (element(ie,i) .EQ. Nodnum(j)) then
			    do ii=1,3
			    do k=1,2
			  Bl(ii,(2*(j-1)+k)) = Bl(ii,(2*(j-1)+k))+B1(ii,(2*(i-1)+k))
			    enddo
			    enddo
			else
		    endif
		    enddo
	  enddo
		NArea=NArea+A11/3
	  else
	  endif
	enddo
C	
	do i = 1,3
		do j = 1,NDOFEL
			Bl(i,j) = (1/NArea)*Bl(i,j)
		enddo
	enddo
c
c Calculation of initial stiffness and tangent stifness matrices
c
        pb=zero
		S=zero
        Fbar=zero
        BB=zero
		GT=zero
		GC=zero
        GCGT=zero
        BBS=zero
        BBSBT=zero
        GB=zero
        pbs=zero
        lam=zero
        mu=zero
        mu1=zero
        dd1=zero
        CMAT=zero
        F=zero
        detF=zero
        CC=zero
        INVCC=zero
        detCC=zero
        dudx = zero
        dudy = zero
        dvdx = zero
        dvdy = zero        
c
		do j = 1,(NDOFEL*half)
			dudx = dudx+Bl(1,(2*j-1))*UNL(2*j-1)
            dudy = dudy+Bl(2,2*j)*UNL(2*j-1)
            dvdx = dvdx+Bl(1,(2*j-1))*UNL(2*j)
            dvdy = dvdy+Bl(2,2*j)*UNL(2*j)
		enddo
c  
		lam=(EMU*E)/((1+EMU)*(1-2*EMU))	
        mu=E/(2*(1+EMU))	

		F(1,1)=dudx+1
        F(1,2)=dudy
        F(2,1)=dvdx
        F(2,2)=dvdy+1
        F(3,3)=1

        CC(1,1)=F(1,1)*F(1,1)+F(2,1)*F(2,1)
        CC(1,2)=F(1,2)*F(1,1)+F(2,1)*F(2,2)
        CC(2,1)=CC(1,2)
        CC(2,2)=F(1,2)*F(1,2)+F(2,2)*F(2,2)
        CC(3,3)=1

        detCC=(CC(1,1)*CC(2,2))-(CC(1,2)*CC(2,1))

        INVCC(1,1)=CC(2,2)/detCC
        INVCC(1,2)=-CC(1,2)/detCC
        INVCC(2,1)=-CC(2,1)/detCC
        INVCC(2,2)=CC(1,1)/detCC
        INVCC(3,3)=1
        
		mu1=mu-lam*half*log(detCC)

	   	DO ii=1,3
	  	DO jj=1,3
        	if (ii .EQ. jj)  then
	   		pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*(INVCC(ii,jj)-1)
        	else
            pb(ii,jj)=lam*half*log(detCC)*INVCC(ii,jj)-mu*INVCC(ii,jj)
          	endif
	  	ENDDO
        ENDDO
        
C Define CMAT matrix for Neohookian material model  
		
		DO x1=1,3
	  	DO j1=1,3
		DO k1=1,3
	  	DO l1=1,3
        dd1(x1,j1,k1,l1)=lam*INVCC(x1,j1)*INVCC(k1,l1)+
     1mu1*(INVCC(x1,k1)*INVCC(j1,l1)+INVCC(x1,l1)*INVCC(j1,k1))
	  	ENDDO
        ENDDO
        ENDDO
        ENDDO
		
	  	CMAT(1,1)=dd1(1,1,1,1)
	  	CMAT(1,2)=dd1(1,1,2,2)
	  	CMAT(1,3)=dd1(1,1,1,2)
	  	CMAT(2,1)=dd1(2,2,1,1)
      	CMAT(2,2)=dd1(2,2,2,2)
	  	CMAT(2,3)=dd1(2,2,1,2)  
      	CMAT(3,1)=dd1(1,2,1,1)
      	CMAT(3,2)=dd1(1,2,2,2)
	  	CMAT(3,3)=dd1(1,2,1,2)
c        
		S(1,1)=pb(1,1)
 		S(1,2)=pb(1,2)
        S(2,1)=pb(2,1)
        S(2,2)=pb(2,2)
        S(3,3)=pb(1,1)
        S(3,4)=pb(1,2)
        S(4,3)=pb(2,1)
        S(4,4)=pb(2,2)

        Fbar(1,1)=F(1,1)
 		Fbar(1,3)=F(2,1)
        Fbar(2,2)=F(1,2)
        Fbar(2,4)=F(2,2)
        Fbar(3,1)=F(1,2)
        Fbar(3,2)=F(1,1)
        Fbar(3,3)=F(2,2)
        Fbar(3,4)=F(2,1)
c
        DO i=1,(NDOFEL*half)
		BB((2*i-1),1)=Bl(1,(2*i-1))
        BB((2*i-1),2)=Bl(2,(2*i))
        BB((2*i),3)=Bl(1,(2*i-1))
        BB((2*i),4)=Bl(2,(2*i))
        ENDDO
c		
		DO i=1,3
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	GT(i,j)=GT(i,j)+Fbar(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
        pbs(1,1)=pb(1,1)
        pbs(2,1)=pb(2,2)
        pbs(3,1)=pb(1,2)
c
		DO i=1,NDOFEL
	  	DO j=1,3
        DO kk1=1,3
	   	GC(i,j)=GC(i,j)+GT(kk1,i)*CMAT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,3
	   	GCGT(i,j)=GCGT(i,j)+GC(i,kk1)*GT(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c
	  	DO i=1,NDOFEL
        DO j=1,4
        DO kk1=1,4
	   	BBS(i,j)=BBS(i,j)+BB(i,kk1)*S(kk1,j)
	  	ENDDO
	  	ENDDO
        ENDDO
c	
		DO i=1,NDOFEL
	  	DO j=1,NDOFEL
        DO kk1=1,4
	   	BBSBT(i,j)=BBSBT(i,j)+BBS(i,kk1)*BB(j,kk1)
	  	ENDDO
	  	ENDDO
        ENDDO
c
		DO i=1,NDOFEL
        DO kk2=1,3
	   	GB(i,1)=GB(i,1)+GT(kk2,i)*pbs(kk2,1)
	  	ENDDO
        ENDDO
c
      	DO i=1,NDOFEL
        KSMATRX(i,1)=KSMATRX(i,1)+GB(i,1)
        DO j=1,NDOFEL
		KTMATRX(i,j)=KTMATRX(i,j)+GCGT(i,j)+BBSBT(i,j)
        ENDDO
		ENDDO
C===============================C
C        calculate  RHS         C
C===============================C
  	  DO K1=1, NDOFEL              !***!   
	    RHS(K1, 1)=RHS(K1, 1)-NArea*KSMATRX(K1,1)
	  ENDDO
c
  	  DO K1=1,NDOFEL              !***!  DEFINING K MATRIX 
	     DO K2=1,NDOFEL              
	     AMATRX(K1, K2)=AMATRX(K1, K2)+NArea*KTMATRX(K1,K2)
	     ENDDO
	  ENDDO       
C---------- Output the saved results           
c         write (*,*) 'Node number',JELEM
c         write (*,*) 'Bmat matrix', GT
c         write (*,*) 'Bgeo matrix transpose', BB 
c         write (*,*) 'F matrix', F
c		  write (*,*) 'dd1 matrix', dd1
c         write (*,*) 'pb matrix', pb
c         write (*,*) 'Cmat matrix', CMAT
c         write (*,*) 'KSMATRX matrix', KSMATRX
c         write (*,*) 'KTMATRX matrix', KTMATRX
c         write (*,*) 'Displacement values', U
c         pause
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
	2b(3,6),A(3,3),DetA,C1
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
	  C1 = abs(DetA)
	  area = C1/2
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