! This program performs numerical calculation of a scalar or vector field described by a partial differential equation of the form:
! D*(ddU/ddx + ddU/ddy + ddU/ddz ) + C*(dU/dt + Vx*dU/dx + Vy*dU/dy + Vz*dU/dz ) = -F
! where U is a scalar or vector; D,C are constant coefficients; Vx, Vy, Vz are constant velocities in x, y, z directions
! F - vector or scalar field of external influences.
! Dirichlet, Neumann (zero only) and Robin (open boundary option only) boundary conditions are provided.
! Initial conditions are provided for the scalar field.
! Movement in 3D space is provided for external field sources.

! Two types of finite-difference approximation are used: for calculation by the method of Successive Over-Relaxation (SOR, matrix-free approximation) 
! and by the method of BiConjugate Gradient STABilized  with rewriting (BiCGSTABwr, a sparse matrix is ​​formed). Fortran code is located in the file solvers.f90

! The initial data is read by the vxc2data subroutine from the working file "in.vxc", which is on created in the "VoxCad" graphic editor.
! The Fortran code is located in the files vxc2data.f90 and m_vxc2data.f90

! Fortran codes created by J.Sochor   ( https://github.com/JNSresearcher )

PROGRAM fields_3d

USE m_fparser
USE m_vxc2data
IMPLICIT NONE

! --------------internal working variables--------------!
! variables to measure calculation time
INTEGER:: k0,k1,k2
REAL(8)   T, Tcalc, Tsavedata

INTEGER nCells,  &                  ! total number of cells
        Ncells0, &                  ! number of cells with scalar field sources
        NcellsX, NcellsY, NcellsZ   ! number of cells with vector field sources along the axes "x", "y", "z"

INTEGER m1, m2,  i,j,k, L, m, n, nn,ios
INTEGER  kim ,kjm, kkm, kip,kjp,kkp, kkpp, kdz
INTEGER  nim,nijm,njm,njkm,nkm,nikm,nijkm, nip,njp,nkp
REAL(8) sxp,sxm,sx,dsx,   syp,sym,sy,dsy,    szp,szm,sz,dsz, a, Ynod(15)  
INTEGER  ncolon(15)

! arrays for new cells of field sources as they move
INTEGER, ALLOCATABLE :: new_nodes0(:), new_nodesX(:), new_nodesY(:), new_nodesZ(:) 
INTEGER  nij, Lnew, jnew, inew, flag_shift, flag_move,  movestop(3)

! arrays of parameters of difference schemes
REAL(8), ALLOCATABLE :: QxschP(:), QyschP(:), QzschP(:), QxschM(:), QyschM(:),QzschM(:),QschP(:)
REAL(8) WschP

! an array of column references for each row and a list of columns for BCG method
!the compressed sparse row  (CSR) format storage
INTEGER, ALLOCATABLE::  irow(:), jcol(:)   
REAL(8), ALLOCATABLE :: valA(:)  
INTEGER  num_nz    ! number of non-zero elements

! variables to control the output of calculation results
INTEGER n_new, Ntime, Npoint, Nout, Nprint, Nprint_display, Nout_display
CHARACTER (LEN = 4)  ch_sd

REAL(8),ALLOCATABLE ::  Ux(:),    Uy(:),    Uz(:),     &    ! arrays for vector fields (dimension nCells)
                        Fx(:),    Fy(:),    Fz(:),     &    ! arrays for  sources  of vector field 
                        Fxbuf(:), Fybuf(:), Fzbuf(:),  &    ! buffer arrays for storing inertial sources of vector field
                        U0(:),                         &    ! arrays for scalar field
                        F0(:),                         &    ! array for sources of scalar field 
                        F0buf(:)                            ! buffer array for storing inertial sources of scalar field
!---------------------------------------!
! input data from SUBROUTINE vxc2data   !
!---------------------------------------!
REAL(8) :: delta(3)         ! array grid spacing along X, Y and Z   
REAL(8) :: dt,   &          ! time step
           Time, &          ! stop time
           dtt              ! jump duration
INTEGER :: sdx,sdy,sdz      ! number of cells along X,  Y and  Z   
INTEGER :: nsub,     &      ! number of physical domains
           nsub_air, &      ! number of environment domains
           nsub_glob        ! total number of domains
INTEGER :: numfun,   &      ! number of functions for calculating external sources
           numMech,  &      ! number of functions for calculating movements of external sources
           num_nodV, &      ! sign of the presence of a vector field
           num_nod0         ! sign of the presence of a scalar field
INTEGER :: size_PHYS_C,  &  ! number of domains with diffusion properties
           size_PHYS_U0, &  ! number of domains with initial potential  (for scalar field)
           size_PHYS_Ubnd   ! number of domains with constant potential (for scalar field)

! solv - character string corresponding to the name of methods: 'BCG' or 'SOR'
CHARACTER(LEN=3) :: solv 
REAL(8)      :: BND(3,2)  ! values ​​of Dirichlet boundary conditions

! character string for boundary conditions on 6 faces: N - Neumann D - Dirichlet A - Absorption
CHARACTER(6) :: bound     

CHARACTER(16):: files      ! name for output files
REAL(8) ::  tolerance      ! convergence criterion
INTEGER ::  iter, itmax    ! actual number of iterations,  maximum number of iterations   
!-------------------------------!
! call to data entry subroutine !
!-------------------------------!
CALL vxc2data ( delta,    dt,        Time,     dtt,                   &
                sdx,    sdy,      sdz,       nsub,     nsub_air, nsub_glob, size_PHYS_C, size_PHYS_U0, size_PHYS_Ubnd,&
                numfun, numMech,  num_nodV, num_nod0,     &
                solv,   files,    tolerance, itmax,    bound,  BND )
!----------------------!
! START of calculation !
!----------------------!
nCells = sdx*sdy*sdz; kdz=sdx*sdy
sz = 1.d0/(delta(3)*delta(3))
sy = 1.d0/(delta(2)*delta(2))
sx = 1.d0/(delta(1)*delta(1)) 
dsx = 0.5d0/delta(1) 
dsy = 0.5d0/delta(2)
dsz = 0.5d0/delta(3)
!---------------------------------------------------------------------------!
! Calculation of parameters for the SOR solver  (Successive OverRelaxation) !
! and formation of a sparse matrix for the BCG solver                       !
! (BiConjugate Gradient stabilized method with rewrite)                     !
!---------------------------------------------------------------------------!
IF      (solv == 'SOR') THEN
    CALL init_SOR                   
ELSEIF  (solv == 'BCG') THEN
    CALL gen_sparse_matrix        ! generation of sparse matrix 
ELSE
   PRINT*, 'ERROR in name solv: "',solv, '", use only SOR or BCG'
   STOP 
ENDIF
!-------------------------------!
! preparation of time intervals !
!-------------------------------!
T=0.
Ntime=0                              ! counter of each point with step dt
Nout_display=nint(Time/(100.0*DT))   ! jump size for display of symbol  ">"
Nprint_display=Ntime + Nout_display  ! counter of points for display symbol ">" with jump Nout_display
Npoint=0                             ! counter of calculated points with step dtt
Nout=nint(DTT/DT)                    ! jump size for output to files
Nprint = Ntime + Nout                ! point counter with step Nout for output to files with step Nout
!-----------------------------------!
! preparing arrays of scalar field  !
!-----------------------------------!
IF (num_nod0 /=0) THEN
    ALLOCATE (F0(nCells), U0(nCells), source=0.d0)
    
    IF ( size_PHYS_C /= 0 ) ALLOCATE (F0buf(nCells), source=0.d0)
    IF (solv == 'SOR' .and. index (bound, 'D') /=0  ) THEN
        CALL boundary(U0,BND)                                    ! Dirichlet conditions for SOR
    ELSEIF (solv == 'BCG' .and.index (bound, 'D') /=0  ) THEN      
        CALL boundary(F0,BND)                                    ! Dirichlet conditions for BCG
    ENDIF
    
    IF (size_PHYS_U0 /=0 ) THEN  ! initial conditions within domains - U0 
        DO i=1, size_PHYS_U0
            n = PHYS_U0(i)%siznod
            a = PHYS_U0(i)%valdom
            DO j=1,n
                k = PHYS_U0(i)%nod(j)
                U0(k) = a
            ENDDO
        ENDDO
    ENDIF

    IF (size_PHYS_Ubnd /=0 ) THEN  ! constant field within domains - Ubnd
        DO i=1, size_PHYS_Ubnd
            n = PHYS_Ubnd(i)%siznod
            a = PHYS_Ubnd(i)%valdom
            DO j=1,n
                k = PHYS_Ubnd(i)%nod(j)
                IF (solv == 'BCG' ) THEN
                    F0(k) = a
                ELSE
                    U0(k) = a
                ENDIF
            ENDDO
        ENDDO
    ENDIF
endif
!-----------------------------------!
! preparing arrays of vector fields !
!-----------------------------------!
IF (num_nodV /=0) THEN
    ALLOCATE (Fx(nCells), Fy(nCells), Fz(nCells), &
              Ux(nCells), Uy(nCells), Uz(nCells), source=0.d0)
    IF ( size_PHYS_C /= 0 ) THEN
        ALLOCATE (Fxbuf(nCells),Fybuf(nCells),Fzbuf(nCells), source=0.d0)
    ENDIF
    IF     (solv == 'SOR' .and. index (bound, 'D') /=0  ) THEN
        CALL boundary(Ux,BND);  CALL boundary(Uy,BND);  CALL boundary(Uz,BND)
    ELSEIF (solv == 'BCG' .and. index (bound, 'D') /=0  ) THEN
        CALL boundary(Fx,BND);  CALL boundary(Fy,BND);  CALL boundary(Fz,BND)
    ENDIF
ENDIF
!---------------------------------------------------------------!
! preparation of arrays for calculating the movement of sources !
!---------------------------------------------------------------!
flag_move = 0 ! global flag of motion in input data

IF (numfun /=0 ) THEN
    DO i=1,numfun
    ! motion search processing
        fun_nod(i)%Distance = 0.d0 
    
        ! unlock
        IF (fun_nod(i)% num_Vmech(1) == 0  .and. fun_nod(i)%move(1) /=0) THEN ! it is not function
        !  distance traveled in time dt with speed Vx in fractions of cell size along dX
            fun_nod(i)%shift(1) = fun_nod(i)%vel_Vmech(1)*dt/delta(1)   !Vx*dt/dX
            flag_move = 1
        ELSEIF (  fun_nod(i)%move(1) /=0) THEN 
            flag_move = 1
        ENDIF
        
        IF (fun_nod(i)% num_Vmech(2) == 0  .and. fun_nod(i)%move(2) /=0) THEN ! it is not function
            fun_nod(i)%shift(2) = fun_nod(i)%vel_Vmech(2)*dt/delta(2)
            flag_move = 1; 
        ELSEIF (  fun_nod(i)%move(2) /=0) THEN 
            flag_move = 1
        ENDIF
    
        IF (fun_nod(i)% num_Vmech(3) == 0 .and. fun_nod(i)%move(3) /=0 ) THEN ! it is not function, but not a shift either
            fun_nod(i)%shift(3) = fun_nod(i)%vel_Vmech(3)*dt/delta(3)
            flag_move = 1; 
        ELSEIF (  fun_nod(i)%move(3) /=0) THEN 
            flag_move = 1
        ENDIF
    ENDDO

! initial allocation of array memory for source nodes, regardless of whether there is movement or not
    Ncells0 = 0; NcellsX = 0; NcellsY = 0; NcellsZ = 0; 
    DO n=1,numfun
        IF     (Fun(n)%ex == 'X') THEN
            NcellsX = NcellsX + fun_nod(n)%numnod_Fx
        ELSEIF (Fun(n)%ex == 'Y') THEN
            NcellsY = NcellsY + fun_nod(n)%numnod_Fy
        ELSEIF (Fun(n)%ex == 'Z') THEN
            NcellsZ = NcellsZ + fun_nod(n)%numnod_Fz
        ELSEIF (Fun(n)%ex == '0') THEN
            Ncells0 = Ncells0 + fun_nod(n)%numnod_F0
        ENDIF 
    ENDDO
    IF (NcellsX /= 0) ALLOCATE (new_nodesX(NcellsX), source=0)
    IF (NcellsY /= 0) ALLOCATE (new_nodesY(NcellsY), source=0)
    IF (NcellsZ /= 0) ALLOCATE (new_nodesZ(NcellsZ), source=0)
    IF (Ncells0 /= 0) ALLOCATE (new_nodes0(Ncells0), source=0)
ENDIF !numfun

! filling source node arrays if there is no movement
IF (flag_move == 0 .and. numfun /=0 ) THEN
    Ncells0 = 0; NcellsX = 0; NcellsY = 0; NcellsZ = 0; ! counters
    DO n=1,numfun
        IF (Fun(n)%ex == 'X') THEN
            DO k = 1,fun_nod(n)%numnod_Fx
                m = fun_nod(n)%nods_Fx(k)
                NcellsX = NcellsX + 1
                new_nodesX(NcellsX) = m
            ENDDO
        ELSEIF (Fun(n)%ex == 'Y') THEN
            DO k = 1,fun_nod(n)%numnod_Fy
                m = fun_nod(n)%nods_Fy(k)
                NcellsY = NcellsY + 1
                new_nodesY(NcellsY) = m
            ENDDO
        ELSEIF (Fun(n)%ex == 'Z') THEN
            DO k = 1,fun_nod(n)%numnod_Fz
                m = fun_nod(n)%nods_Fz(k)
                NcellsZ = NcellsZ + 1
                new_nodesZ(NcellsZ) = m
            ENDDO
        ELSEIF (Fun(n)%ex == '0') THEN
            DO k = 1, fun_nod(n)%numnod_F0
                m = fun_nod(n)%nods_F0(k)
                F0(m) = a
                Ncells0 = Ncells0 + 1
                new_nodes0(Ncells0) = m
            ENDDO
        ELSE
            STOP 'err calc filling  array src, if not move'
        ENDIF  
    ENDDO
ENDIF

! CALL system ( 'mkdir out_'//trim(files), i )
! PRINT *, 'create dir "out_'//trim(files)//'": status=', i

CALL system ( 'mkdir out', i )
PRINT *, 'create dir "out": status=', i
!!!!!!  START time counter
CALL system_clock(k1,k0)
Tsavedata = 0.0
movestop = 1  ! flag of local start-stop movement in all coordinates
write(*,'( "start solver: ", a)') solv

2000 CONTINUE   
!-------------------------------------!
! Calculation of independent sources  !
!-------------------------------------!
IF (numfun /= 0) THEN
    CALL initf (numfun) 
    DO i=1,numfun
        DO m=1,Fun(i)%args
            IF ( trim(Fun(i)%namex(m)) == 'T') THEN
                Fun(i)%velx(m)=T
            ENDIF
        ENDDO
        CALL parsef (i, trim(fun(i)%eqn), fun(i)%namex(1:fun(i)%args) )
        fun(i)%vely = evalf (i, fun(i)%velx(1:fun(i)%args))
    ENDDO
ENDIF
!-------------------------------------------------!
! Calculation of the speed of movement of sources !
!-------------------------------------------------!
IF (numMech /= 0) THEN
    CALL initf (numMech) 
    DO i=1,numMech
        DO m=1,Vmech(i)%args
            IF ( trim(Vmech(i)%namex(m)) == 'T') THEN
                Vmech(i)%velx(m)=T
            ENDIF
        ENDDO
        CALL parsef (i, trim(Vmech(i)%eqn), Vmech(i)%namex(1:Vmech(i)%args) )
        Vmech(i)%vely = evalf (i, Vmech(i)%velx(1:Vmech(i)%args))
    ENDDO
ENDIF
!--------------------------------------------------!
! Distribution of independent sources into domains !
!--------------------------------------------------!
IF (num_nodV /=0) THEN

    IF (flag_move == 1) THEN 
        ! there are moving sources ! 
        IF (size_PHYS_C /= 0  ) THEN      ! only inertial sources 
            Fxbuf =0.d0; Fybuf =0.d0; Fzbuf =0.d0
            
            ! recording inertial sources to the buffer
            DO m=1,size_PHYS_C 
                n = PHYS_C(m)%siznod
                Fxbuf( PHYS_C(m)%nod(1:n) ) =  Fx( PHYS_C(m)%nod(1:n) )
                Fybuf( PHYS_C(m)%nod(1:n) ) =  Fy( PHYS_C(m)%nod(1:n) )
                Fzbuf( PHYS_C(m)%nod(1:n) ) =  Fz( PHYS_C(m)%nod(1:n) )
            ENDDO
            ! reset all sources
            Fx = 0.d0; Fy = 0.d0; Fz = 0.d0; 
            ! restore of inertial sources
            Fx = Fxbuf; Fy = Fybuf; Fz = Fzbuf; ! only C
        ELSE
            Fx = 0.d0; Fy = 0.d0; Fz = 0.d0; 
        ENDIF
        
        ! calculation of new coordinates of nodes of independent sources 
        IF (numfun /= 0) THEN
            NcellsX=0; NcellsY=0; NcellsZ=0
            
            DO n=1,numfun
                CALL motion_calc    ! motion calculation
            
                a = Fun(n)%vely
                IF (Fun(n)%ex == 'X') THEN
                    DO k = 1,fun_nod(n)%numnod_Fx
                        m = fun_nod(n)%nods_Fx(k)
                        CALL new_m !(m)
                        Fx(m) = a
                        
                        NcellsX = NcellsX + 1
                        new_nodesX(NcellsX) = m
                    ENDDO
                ELSEIF (Fun(n)%ex == 'Y') THEN
                    DO k = 1,fun_nod(n)%numnod_Fy
                        m = fun_nod(n)%nods_Fy(k)
                        CALL new_m !(m)
                        Fy(m) = a
                        
                        NcellsY = NcellsY + 1
                        new_nodesY(NcellsY) = m
                    ENDDO
                ELSEIF (Fun(n)%ex == 'Z') THEN
                    DO k = 1,fun_nod(n)%numnod_Fz
                        m = fun_nod(n)%nods_Fz(k)
                        CALL new_m !(m)
                        Fz(m) = a
                        
                        NcellsZ = NcellsZ + 1
                        new_nodesZ(NcellsZ) = m
                    ENDDO
                ELSEIF (Fun(n)%ex == '0') THEN
                ELSE
                    STOP 'err calc Fx or Fy or Fz '
                ENDIF                    
            ENDDO
        ENDIF
    ELSE
        !--------------------!
        ! no moving sources  !
        !--------------------!
        IF (numfun /= 0) THEN
            DO n=1,numfun
                a = Fun(n)%vely
                IF (Fun(n)%ex == 'X') THEN
                    DO k = 1,fun_nod(n)%numnod_Fx
                        m = fun_nod(n)%nods_Fx(k)
                        Fx(m) = a
                    ENDDO
                ELSEIF (Fun(n)%ex == 'Y') THEN
                    DO k = 1,fun_nod(n)%numnod_Fy
                        m = fun_nod(n)%nods_Fy(k)
                        Fy(m) = a
                    ENDDO
                ELSEIF (Fun(n)%ex == 'Z') THEN
                    DO k = 1,fun_nod(n)%numnod_Fz
                        m = fun_nod(n)%nods_Fz(k)
                        Fz(m) = a
                    ENDDO
                ELSEIF (Fun(n)%ex == '0') THEN
                ELSE
                    STOP 'err calc Fx or Fy or Fz '
                ENDIF  
            ENDDO
        ENDIF
    ENDIF

    IF ( index (bound, 'N') /=0 .or. index (bound, 'A') /=0 ) THEN
        CALL boundary(Ux, BND); CALL boundary(Uy, BND);  CALL boundary(Uz, BND)
    ENDIF

    IF (size_PHYS_C /=0) THEN
    !---------------------------------!
    ! calculation of inertial sources !
    !---------------------------------!
        DO m=1,size_PHYS_C 
            a = PHYS_C(m)%valdom
            DO j=1, PHYS_C(m)%siznod
                k = PHYS_C(m)%nod(j)
                Fx(k) = a * Ux(k) +  Fx(k) 
                Fy(k) = a * Uy(k) +  Fy(k)   ! Jb = 2С/dt * U_из + Ib  
                Fz(k) = a * Uz(k) +  Fz(k)
            ENDDO
        ENDDO
    ENDIF
                                  !---------------------------!
                                  ! SOLVERS for vector fields !
                                  !---------------------------!
    IF     (solv == 'BCG') THEN
        IF (numfun ==0 .and.  (size_PHYS_U0 /=0 .or. size_PHYS_Ubnd /=0 .or.  count(BND == 0.d0) /=6 ) ) THEN
            ! pure boundary value problem (RHS of linear equations = 0)
            CALL sprsBCGstabwr2(valA, irow, jcol, nCells, Fx, Ux,  tolerance, itmax, iter)
            CALL sprsBCGstabwr2(valA, irow, jcol, nCells, Fy, Uy,  tolerance, itmax, iter)
            CALL sprsBCGstabwr2(valA, irow, jcol, nCells, Fz, Uz,  tolerance, itmax, iter)
        ELSE
            ! RHS of linear equations /= 0
            CALL sprsBCGstabwr (valA, irow, jcol, nCells, Fx, Ux,  tolerance, itmax, iter)
            CALL sprsBCGstabwr (valA, irow, jcol, nCells, Fy, Uy,  tolerance, itmax, iter)
            CALL sprsBCGstabwr (valA, irow, jcol, nCells, Fz, Uz,  tolerance, itmax, iter)
        ENDIF
    ELSEIf (solv == 'SOR') THEN
        CALL SOR(sdx,sdy,sdz, geoPHYS, geoUbnd, QschP,WschP, QxschM, QyschM,QzschM, QxschP, QyschP,QzschP, Fx, Ux, tolerance, itmax )
        CALL SOR(sdx,sdy,sdz, geoPHYS, geoUbnd, QschP,WschP, QxschM, QyschM,QzschM, QxschP, QyschP,QzschP, Fy, Uy, tolerance, itmax )
        CALL SOR(sdx,sdy,sdz, geoPHYS, geoUbnd, QschP,WschP, QxschM, QyschM,QzschM, QxschP, QyschP,QzschP, Fz, Uz, tolerance, itmax )
    ENDIF

    IF (size_PHYS_C /=0) THEN
    !---------------------------------------------!
    ! calculation of sources of inertial branches !
    !---------------------------------------------!
        DO m=1,size_PHYS_C 
            a = PHYS_C(m)%valdom
            DO j=1, PHYS_C(m)%siznod
                k = PHYS_C(m)%nod(j)
                Fx(k) = a * Ux(k) -  Fx(k)  
                Fy(k) = a * Uy(k) -  Fy(k)    ! Ib = 2С/dt * Ub - Jb 
                Fz(k) = a * Uz(k) -  Fz(k) 
            ENDDO
        ENDDO
    ENDIF
ENDIF ! num_nodV /=0

IF (num_nod0 /=0) THEN
    IF (flag_move == 1) THEN  ! exist shift 
            !  delete inerc
        IF (size_PHYS_C /= 0  ) THEN  ! only C to buffer
            F0buf =0.d0
            DO m=1,size_PHYS_C 
                n = PHYS_C(m)%siznod
                DO j=1,n
                    k=PHYS_C(m)%nod(j)
                    F0buf( k ) =  F0( k ) 
                ENDDO
                ! F0buf( PHYS_C(m)%nod(1:n) ) =  F0( PHYS_C(m)%nod(1:n) ) 
            ENDDO
            F0 = 0.d0
            F0 =  F0buf ! only C
        ELSE
            F0 = 0.d0
        ENDIF

        IF (numfun /= 0) THEN
            Ncells0 = 0
            DO n=1,numfun
                CALL motion_calc    ! motion calculation
                a = Fun(n)%vely
                IF (Fun(n)%ex == '0') THEN
                    DO k = 1, fun_nod(n)%numnod_F0
                        m = fun_nod(n)%nods_F0(k)
                        CALL new_m !(m)
                        F0(m) = a
                        Ncells0 = Ncells0 + 1
                        new_nodes0(Ncells0) = m
                    ENDDO
                ELSEIF (Fun(n)%ex == 'X') THEN
                ELSEIF (Fun(n)%ex == 'Y') THEN
                ELSEIF (Fun(n)%ex == 'Z') THEN
                ELSE
                    STOP 'err calc F0 '
                ENDIF  
            ENDDO
        ENDIF
    ELSE  ! no  delete inerc
        IF (numfun /= 0) THEN
            DO n=1,numfun
                a = Fun(n)%vely
                IF (Fun(n)%ex == '0') THEN
                    DO k = 1, fun_nod(n)%numnod_F0
                        m = fun_nod(n)%nods_F0(k)
                        F0(m) = a
                    ENDDO
                ELSEIF (Fun(n)%ex == 'X') THEN
                ELSEIF (Fun(n)%ex == 'Y') THEN
                ELSEIF (Fun(n)%ex == 'Z') THEN
                ELSE
                    STOP 'err calc F0 '
                ENDIF 
            ENDDO
        ENDIF
    ENDIF
    
    IF ( index (bound, 'N') /=0 .or. index (bound, 'A') /=0) CALL boundary(U0, BND )
    
    IF (size_PHYS_C /=0) THEN
       !---------------------------------!
       ! calculation of inertial sources !
       !---------------------------------! 
        DO m=1,size_PHYS_C 
            a = PHYS_C(m)%valdom
            DO j=1, PHYS_C(m)%siznod
                k = PHYS_C(m)%nod(j)
                F0(k) = a * U0(k) +  F0(k) !   i + U0*2C/dt + i0 =  U*2C/dt 
            ENDDO
        ENDDO
        
        IF     (solv == 'BCG') THEN
            IF (   index (bound, 'D') /=0 ) CALL boundary(F0, BND ) 
            IF (size_PHYS_Ubnd /=0 ) THEN  ! initial U0
                DO i=1, size_PHYS_Ubnd
                    n = PHYS_Ubnd(i)%siznod
                    a = PHYS_Ubnd(i)%valdom
                    DO j=1,n
                        k = PHYS_Ubnd(i)%nod(j)
                        F0(k) = a
                    ENDDO
                ENDDO
            ENDIF
        ENDIF
    ENDIF
                                  !---------------------------!
                                  ! SOLVERS for scalar fields !
                                  !---------------------------!
    IF     (solv == 'BCG') THEN  
        IF (numfun ==0 .and.  (size_PHYS_U0 /=0 .or. size_PHYS_Ubnd /=0 .or.  count(BND == 0.d0) /=6 ) ) THEN
            ! only boundary value problem (RHS of linear equations = 0)
            CALL sprsBCGstabwr2(valA, irow, jcol, nCells, F0, U0,  tolerance, itmax, iter)  ! pure bnd
        ELSE
            ! RHS of linear equations /= 0
            CALL sprsBCGstabwr (valA, irow, jcol, nCells, F0, U0,  tolerance, itmax, iter)
        ENDIF
    ELSEIF (solv == 'SOR') THEN
        CALL SOR(sdx,sdy,sdz, geoPHYS, geoUbnd, QschP,WschP, QxschM, QyschM,QzschM, QxschP, QyschP,QzschP, F0, U0, tolerance, itmax )
    ENDIF

    IF (size_PHYS_C /=0) THEN
    !---------------------------------------------!
    ! calculation of sources of inertial branches !
    !---------------------------------------------!
        DO m=1,size_PHYS_C 
            a = PHYS_C(m)%valdom
            DO j=1, PHYS_C(m)%siznod
                k = PHYS_C(m)%nod(j)
                F0(k) = a * U0(k) -  F0(k) 
            ENDDO
        ENDDO
    ENDIF
ENDIF !num_nod0 /=0

!===========================================================================================
IF (Ntime >= Nprint .and. Ntime /=0 ) THEN  ! output with step dtt skip 1st point
    Nprint = Ntime  + Nout
    Npoint= Npoint + 1                      ! point counter with step dtt
    ios=0
    CALL writeVtk_field (Npoint, num_nod0, num_nodV, sdx, sdy, sdz, delta,  U0, F0, Ux,Uy,Uz, Fx,Fy,Fz, geoPHYS, typPHYS, size_PHYS_C, files)
    ! IF (flag_move == 1) THEN
    IF (numfun /= 0) CALL writeVtk_src ( Npoint, numfun,  Ncells0, NcellsX, NcellsY, NcellsZ, new_nodes0, new_nodesX, new_nodesY,  new_nodesZ,  &
                          num_nod0, num_nodV,  sdx, sdy, sdz, delta, files)
    ! ENDIF
ENDIF

IF (Ntime >= Nprint_display) THEN
    Nprint_display = Ntime + Nout_display

    write(*,'( a,$ )') '>'
ENDIF

Ntime = Ntime + 1   ! this is every point
T = T + DT

IF (T < Time) GOTO 2000
PRINT*,'|'

CALL system_clock(k2);
Tcalc =  REAL(k2-k1)/REAL(k0)
   
write(*,'( "solve complet. Tcalc= ", g10.3)') Tcalc

CONTAINS

SUBROUTINE gen_sparse_matrix  
! GENERATION of compressed sparse matrix discretizing the Laplace operator
! d^2 / dx^2 + d^2 / dy^2 + d^2 / dz^2  is created from a 7-point
!  stencil on an Nx by Ny by Nz grid.

! This is an approximate version! The boundary conditions at corners, edges, and faces are not fully taken into account.
num_nz = 0
DO k = 1, sdz;  DO j = 1, sdy;  DO i = 1, sdx
    IF (i==1 .or.j==1 .or.k==1 .or.i==sdx .or.j==sdy .or.k==sdz  ) THEN
!boundary X
        IF      (i==1   ) THEN
            IF     ( bound(1:1) == 'N' .or. bound(1:1) == 'A' ) THEN 
                num_nz = num_nz + 2
            ELSEIF ( bound(1:1) == 'D') THEN
                num_nz = num_nz + 1
            ENDIF
        ELSEIF (i==sdx ) THEN
            IF     ( bound(2:2) == 'N' .or. bound(1:1) == 'A' ) THEN
                num_nz = num_nz + 2
            ELSEIF ( bound(2:2) == 'D') THEN
                num_nz = num_nz + 1
            ENDIF
!boundary Y  
        ELSEIF (j==1   )  THEN 
            IF     ( bound(3:3) == 'N' .or. bound(1:1) == 'A' ) THEN
                num_nz = num_nz + 2
            ELSEIF ( bound(3:3) == 'D') THEN 
                num_nz = num_nz + 1
            ENDIF
        ELSEIF (j==sdy  ) THEN
            IF     ( bound(4:4) == 'N' .or. bound(1:1) == 'A' ) THEN
                num_nz = num_nz + 2
            ELSEIF ( bound(4:4) == 'D') THEN
                num_nz = num_nz + 1
            ENDIF
!boundary Z
        ELSEIF (k==1   ) THEN  
            IF     ( bound(5:5) == 'N' .or. bound(1:1) == 'A' ) THEN
                num_nz = num_nz + 2
            ELSEIF ( bound(5:5) == 'D') THEN
                num_nz = num_nz + 1
            ENDIF
        ELSEIF (k==sdz ) THEN
            IF     ( bound(6:6) == 'N' .or. bound(1:1) == 'A' ) THEN
                num_nz = num_nz + 2
            ELSEIF ( bound(6:6) == 'D') THEN
                num_nz = num_nz + 1
            ENDIF
        ENDIF
    ELSEIF (geoUbnd(i,j,k) /= 0) THEN
        num_nz = num_nz + 1
    ELSE
        num_nz = num_nz + 7
    ENDIF
ENDDO; ENDDO; ENDDO;

! memory allocation for CSR storage format arrays
ALLOCATE (valA(num_nz), source=0.d0)
ALLOCATE (jcol(num_nz), irow(nCells+1), source=0)

nn = 0; L=0; m=0
DO k=1,sdz; DO j=1,sdy;   DO i=1,sdx
    n = geoPHYS(i,j,k)   
    nn =  nn + 1    
    kim = nn-1;  kjm = nn-sdx;  kkm = nn-kdz;
    kip = nn+1;  kjp = nn+sdx;  kkp = nn+kdz; 
        
    irow(nn)=L+1
    IF (i==1 .or.j==1 .or.k==1 .or.i==sdx .or.j==sdy .or.k==sdz ) THEN
!boundary X
        IF      (i==1   ) THEN
            IF     ( bound(1:1) == 'N' .or. bound(1:1) == 'A' ) THEN
                L = L + 1; jcol(L) = nn;  valA(L) = 1.d0; 
                L = L + 1; jcol(L) = kip; valA(L) = -1.d0 
            ELSEIF ( bound(1:1) == 'D') THEN
                L = L + 1; jcol(L) = nn; valA(L) = 1.d0
            ENDIF
        ELSEIF (i==sdx ) THEN
            IF     ( bound(2:2) == 'N' .or. bound(2:2) == 'A' ) THEN 
                L = L + 1; jcol(L) = kim; valA(L) = -1.d0
                L = L + 1; jcol(L) = nn;  valA(L) = 1.d0;  
            ELSEIF ( bound(2:2) == 'D') THEN 
                L = L + 1;  jcol(L) = nn; valA(L) = 1.d0 
            ENDIF
!boundary Y  
        ELSEIF (j==1   )  THEN 
            IF     ( bound(3:3) == 'N' .or. bound(3:3) == 'A' ) THEN
                L = L + 1; jcol(L) = nn;  valA(L) = 1.d0;
                L = L + 1; jcol(L) = kjp; valA(L) = -1.d0 
            ELSEIF ( bound(3:3) == 'D') THEN 
                L = L + 1; jcol(L) = nn; valA(L) = 1.d0
            ENDIF
        ELSEIF (j==sdy  ) THEN
            IF     ( bound(4:4) == 'N' .or. bound(4:4) == 'A' ) THEN
                L = L + 1; jcol(L) = kjm; valA(L) = -1.d0 
                L = L + 1; jcol(L) = nn;  valA(L) = 1.d0; 
            ELSEIF ( bound(4:4) == 'D') THEN 
                L = L + 1; jcol(L) = nn; valA(L) = 1.d0
            ENDIF
!boundary Z
        ELSEIF (k==1   ) THEN  
            IF     ( bound(5:5) == 'N' .or. bound(5:5) == 'A' ) THEN 
                L = L + 1; jcol(L) = nn;  valA(L) = 1.d0; 
                L = L + 1; jcol(L) = kkp; valA(L) = -1.d0
            ELSEIF ( bound(5:5) == 'D') THEN 
                L = L + 1; jcol(L) = nn; valA(L) = 1.d0;
            ENDIF
        ELSEIF (k==sdz ) THEN
            IF     ( bound(6:6) == 'N' .or. bound(6:6) == 'A' ) THEN 
                L = L + 1; jcol(L) = kkm; valA(L) = -1.d0 
                L = L + 1; jcol(L) = nn;  valA(L) = 1.d0; 
            ELSEIF ( bound(6:6) == 'D') THEN ! V(n) = B0(3,2)
                L = L + 1; jcol(L) = nn; valA(L) = 1.d0
            ENDIF
        ENDIF
    ELSEIF (geoUbnd(i,j,k) /= 0) THEN
        L = L + 1; jcol(L) = nn; valA(L) = 1.d0
    ELSE
        nim = geoPHYS(i-1,j,k);  nijm = geoPHYS(i-1,j-1,k);
        njm = geoPHYS(i,j-1,k);  njkm = geoPHYS(i,j-1,k-1);
        nkm = geoPHYS(i,j,k-1);  nikm = geoPHYS(i-1,j,k-1); 
        nijkm = geoPHYS(i-1,j-1,k-1); 

        sxm= sx * 0.25d0*(valPHYS(nim,1)  + valPHYS(nijm,1)  + valPHYS(nijkm,1) + valPHYS(nikm,1) ) + valPHYS(n,3)/(2.d0*delta(1)) !  Vex
        sxp= sx * 0.25d0*(valPHYS(n,1)    + valPHYS(njm,1)   + valPHYS(njkm,1)  + valPHYS(nkm,1) )  - valPHYS(n,3)/(2.d0*delta(1)) 
        sym= sy * 0.25d0*(valPHYS(nijm,1) + valPHYS(nijkm,1) + valPHYS(njkm,1)  + valPHYS(njm,1) )  + valPHYS(n,4)/(2.d0*delta(2)) !  Vey 
        syp= sy * 0.25d0*(valPHYS(n,1)    + valPHYS(nim,1)   + valPHYS(nikm,1)  + valPHYS(nkm,1) )  - valPHYS(n,4)/(2.d0*delta(2)) 
        szm= sz * 0.25d0*(valPHYS(nkm,1)  + valPHYS(nikm,1)  + valPHYS(nijkm,1) + valPHYS(njkm,1) ) + valPHYS(n,5)/(2.d0*delta(3)) !  Vez
        szp= sz * 0.25d0*(valPHYS(n,1)    + valPHYS(nim,1)   + valPHYS(nijm,1)  + valPHYS(njkm,1))  - valPHYS(n,5)/(2.d0*delta(3))  
        a = (sxm+sxp + sym+syp + szm+szp) + 2.d0*valPHYS(n,2)/dt 
        jcol(L+1:L+7) = [ kkm,  kjm,  kim, nn,  kip,  kjp,  kkp ]
        valA(L+1:L+7) = [-szm, -sym, -sxm,  a, -sxp, -syp, -szp ]
        L = L + 7
    ENDIF  

ENDDO ;  ENDDO; ENDDO
irow(nn+1)=L+1

PRINT '(a,i3,a,i3,a,i3,     a,i9, a,g12.5,a/)', 'sparse matrix generation completed on grid (', sdx,' x ',sdy,' x ',sdz, &
                 ' ), Non zero elem= ',num_nz, ' Density of matrix:', 100.0* REAL(num_nz)/REAL(nCells)/REAL(nCells),'%'

END SUBROUTINE gen_sparse_matrix

SUBROUTINE boundary(V, B0) 
! This is an approximate version! The boundary conditions at corners, edges, and faces are not fully taken into account.
REAL(8)  :: V(nCells), B0(3,2)
    n=0 
    DO k = 1,sdz;  DO j = 1,sdy;   DO i = 1,sdx
        n=n+1
        kjm = n-sdx;  kkm = n-kdz;
        kjp = n+sdx;  kkp = n+kdz; 
!boundary X
        IF      (i==1   ) THEN
            IF     ( bound(1:1) == 'N') THEN 
                V(n) = V(n+1)
            ELSEIF ( bound(1:1) == 'D') THEN
                V(n) = B0(1,1)
            ELSEIF ( bound(1:1) == 'A') THEN
                V(n) = 0.9d0*V(n+1)
            ENDIF
        ELSEIF (i==sdx ) THEN
            IF     ( bound(2:2) == 'N') THEN
                V(n) = V(n-1)
            ELSEIF ( bound(2:2) == 'D') THEN
                V(n) = B0(1,2)
            ELSEIF ( bound(2:2) == 'A') THEN
                V(n) = 0.9d0*V(n-1)
            ENDIF
!boundary Y  
        ELSEIF (j==1   )  THEN 
            IF     ( bound(3:3) == 'N') THEN
                V(n) = V(kjp)
            ELSEIF ( bound(3:3) == 'D') THEN
                V(n) = B0(2,1)
            ELSEIF ( bound(3:3) == 'A') THEN
                V(n) = 0.9d0*V(kjp)
            ENDIF
        ELSEIF (j==sdy  ) THEN
            IF     ( bound(4:4) == 'N') THEN
                V(n) = V(kjm)
            ELSEIF ( bound(4:4) == 'D') THEN
                V(n) = B0(2,2)
            ELSEIF ( bound(4:4) == 'A') THEN
                V(n) = 0.9d0*V(kjm)
            ENDIF
!boundary Z
        ELSEIF (k==1   ) THEN  
            IF     ( bound(5:5) == 'N') THEN
                V(n) = V(kkp) 
            ELSEIF ( bound(5:5) == 'D') THEN
                V(n) = B0(3,1)
            ELSEIF ( bound(5:5) == 'A') THEN
                V(n) = 0.9d0*V(kkp)
            ENDIF
        ELSEIF (k==sdz ) THEN
            IF     ( bound(6:6) == 'N') THEN
                V(n) = V(kkm)  
            ELSEIF ( bound(6:6) == 'D') THEN
                V(n) = B0(3,2)
            ELSEIF ( bound(6:6) == 'A') THEN
                V(n) = 0.9d0*V(kkm)
            ENDIF
        ENDIF
     ENDDO; ENDDO;ENDDO
END SUBROUTINE boundary


SUBROUTINE init_SOR
    ALLOCATE ( QxschP(nCells), QyschP(nCells),QzschP(nCells),QxschM(nCells), QyschM(nCells),QzschM(nCells),QschP(nCells), source=0.d0 )
    IF     (sdx <= sdy .and. sdx <= sdz ) THEN
        i=sdy; j=sdz
    ELSEIF (sdy <= sdx .and. sdy <= sdz ) THEN
        i=sdx; j=sdz
    ELSEIF (sdz <= sdx .and. sdz <= sdy ) THEN
        i=sdx; j=sdy
    ELSE
        STOP 'no find sdx ,sdy , sdz'
    ENDIF

    a = (COS( 3.1415926_8/i ) +  COS( 3.1415926_8/j ))/2.0_8 
    a = 2.0_8/(1.0_8 + sqrt(1.0_8 - a*a)   )
    WschP = a   
    
! filling array
DO k=2,sdz-1 ; DO j=2,sdy-1;  DO i=2,sdx-1
    n=geoPHYS(i,j,k);
    m   = i   + sdx*(j-1)   + sdx*sdy*(k-1)
    
    nim = geoPHYS(i-1,j,k);  nijm = geoPHYS(i-1,j-1,k);
    njm = geoPHYS(i,j-1,k);  njkm = geoPHYS(i,j-1,k-1);
    nkm = geoPHYS(i,j,k-1);  nikm = geoPHYS(i-1,j,k-1); 
    nijkm = geoPHYS(i-1,j-1,k-1); 

QxschM(m) = sx * 0.25d0*(valPHYS(nim,1)  + valPHYS(nijm,1)  + valPHYS(nijkm,1) + valPHYS(nikm,1) ) + valPHYS(n,3)/(2.d0*delta(1)) !  Vex
QxschP(m) = sx * 0.25d0*(valPHYS(n,1)    + valPHYS(njm,1)   + valPHYS(njkm,1)  + valPHYS(nkm,1) )  - valPHYS(n,3)/(2.d0*delta(1)) 
QyschM(m) = sy * 0.25d0*(valPHYS(nijm,1) + valPHYS(nijkm,1) + valPHYS(njkm,1)  + valPHYS(njm,1) )  + valPHYS(n,4)/(2.d0*delta(2)) !  Vey 
QyschP(m) = sy * 0.25d0*(valPHYS(n,1)    + valPHYS(nim,1)   + valPHYS(nikm,1)  + valPHYS(nkm,1) )  - valPHYS(n,4)/(2.d0*delta(2)) 
QzschM(m) = sz * 0.25d0*(valPHYS(nkm,1)  + valPHYS(nikm,1)  + valPHYS(nijkm,1) + valPHYS(njkm,1) ) + valPHYS(n,5)/(2.d0*delta(3)) !  Vez
QzschP(m) = sz * 0.25d0*(valPHYS(n,1)    + valPHYS(nim,1)   + valPHYS(nijm,1)  + valPHYS(njkm,1) ) - valPHYS(n,5)/(2.d0*delta(3))  

    QschP(m) = QxschP(m)+QxschM(m) + QyschP(m)+QyschM(m) + QzschP(m)+QzschM(m) + 2.d0*valPHYS(n,2)/dt ! c/dt

ENDDO; ENDDO;  ENDDO

PRINT '(a,i3,a,i3,a,i3,     a,i9,      a,g12.5)', 'SOR init completed on grid (', sdx,' x ',sdy,' x ',sdz, &
                 ' ), itmax= ',itmax, ' tolerance=', tolerance
END SUBROUTINE init_SOR

SUBROUTINE motion_calc 
    DO i=1,3
        IF (fun_nod(n)% num_Vmech(i) == 0 ) THEN     ! это не фун 
            fun_nod(n)%Distance(i) = fun_nod(n)%Distance(i) + movestop(1)*fun_nod(n)%shift(i)
            fun_nod(n)%length(i) = nint(fun_nod(n)%Distance(i))
        ELSE                                         ! это фун  
            fun_nod(n)%Distance(i) = fun_nod(n)%Distance(i) +  Vmech(fun_nod(n)%num_Vmech(i))%vely*dt/delta(i)
            fun_nod(n)%length(i) = nint(fun_nod(n)%Distance(i))
        ENDIF
    ENDDO
END SUBROUTINE motion_calc 

SUBROUTINE new_m 
    L = ceiling(REAL(m)/( REAL(sdx*sdy) ) ) 
    Lnew = L + fun_nod(n)%length(3)
    
    IF  ( Lnew > sdz-2  ) THEN 
        movestop(3) =0; Lnew = sdz-2
    ELSEIF(Lnew < 2    ) THEN
        movestop(3) =0;  Lnew = 2
    ELSEIF ( movestop(3) == 0 .and. (Lnew < sdz-2 .or. Lnew > 2)  ) THEN
             movestop(3) = 1 
    ENDIF
    IF (L == 1) THEN
        nij = m
    ELSE
        nij = m - (L-1)*sdx*sdy
    ENDIF
    j = ceiling( REAL(nij) / REAL(sdx) )  
    jnew = j + fun_nod(n)%length(2)
!-------------------------------------------------------------------------------------------
! comment/uncomment to check for out of bounds along the y-axis
!  1 variant
    IF  ( jnew > sdy-2  ) THEN 
        movestop(2) =0; jnew = sdy-2
    ELSEIF(jnew < 2    ) THEN
        movestop(2) =0;  jnew = 2
    ELSEIF ( movestop(2) == 0 .and. (jnew < sdy-2 .or. jnew > 2)  ) THEN
             movestop(2) = 1 
    ENDIF
!======================================================
!  2 variant
    ! IF  ( jnew > sdy  ) THEN 
        ! movestop(2) =0; jnew = sdy
    ! ELSEIF(jnew < 0    ) THEN
        ! movestop(2) =0;  jnew = 1
    ! ELSEIF ( movestop(2) == 0 .and. (jnew < sdy .or. jnew > 0 )  ) THEN
             ! movestop(2) = 1 
    ! ENDIF
!-----------------------------------------------------------------------------------------
    i = nij - (j - 1) * sdx  
    inew = i + fun_nod(n)%length(1)

    IF  ( inew > sdx-2  ) THEN  !
        movestop(1) =0;  inew = sdx-2
    ELSEIF(inew < 2    ) THEN !
        movestop(1) =0; inew = 2
    ELSEIF ( movestop(1) == 0 .and. (inew < sdx-2 .or. inew > 2)  ) THEN
             movestop(1) = 1 
    ENDIF

    m  = inew + sdx*(jnew-1) + sdx*sdy*(Lnew-1)
END SUBROUTINE new_m

END program fields_3d


