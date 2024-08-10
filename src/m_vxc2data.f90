! module with declarations of global variables, arrays and structures

! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )

MODULE m_vxc2data
implicit none

! structures for calculating functions
TYPE tFun
    CHARACTER (LEN=1)               :: ex           ! source direction identifier: X, Y, Z or 0
    REAL(8)                         :: vely         ! function value
    CHARACTER (LEN=50)              :: eqn          ! function expression
    INTEGER                         :: args,nomsch  ! number of arguments and domain number
    CHARACTER (LEN=8),ALLOCATABLE   :: namex(:)     ! argument names
    REAL(8),ALLOCATABLE             :: velx(:)      ! argument values
END TYPE tFun
    TYPE (tFun), ALLOCATABLE :: Fun(:), Vmech(:)    !  structures of functions for sources and functions for sources movement
! U0bnd(:) in future - for bound scalar field 

! structure for cells in which functions are located
TYPE tfun_nod
!  number of cells for vector and scalar functions
    INTEGER             :: numnod_Fx,  numnod_Fy,  numnod_Fz,  numnod_F0
    
! arrays with cells numbers for vector and scalar functions
    INTEGER,ALLOCATABLE :: nods_Fx(:), nods_Fy(:), nods_Fz(:), nods_F0(:)

! speed, distance and relative shift of the source movement along the X Y Z axes
    REAL(8)             :: vel_Vmech(3),  Distance(3),  shift(3)  

    INTEGER             :: num_Vmech(3),  &    ! X Y Z axis motion function number
                           move(3)             ! sign of motion in input data
                           
    INTEGER             :: length(3)           ! distance in integers: length = nint(Distance)
END TYPE tfun_nod
    TYPE (tfun_nod), ALLOCATABLE :: fun_nod(:)  

! structure for domain parameters
TYPE tPHYS 
    INTEGER, ALLOCATABLE  :: nod(:)     ! arrays of nodes
    INTEGER :: numdom,       siznod     ! size
    REAL(8) :: valdom                   ! value
END TYPE tPHYS
TYPE (tPHYS), ALLOCATABLE :: PHYS_C(:),    & ! inertial domains
                             PHYS_Ubnd(:), & ! domains with a constant scalar field
                             PHYS_U0(:)      ! domains with initial conditions for scalar field
! INTEGER :: size_PHYS_C,  size_PHYS_Ubnd, size_PHYS_U0

INTEGER(1),        ALLOCATABLE  :: geoPHYS(:,:,:),   &  ! array with physical domain numbers
                                   geoUbnd(:,:,:)       ! array of domains  with constant scalar field

! valPHYS - array of domain parameters
! valPHYS(:,1) - "diffusion" parameter D
! valPHYS(:,2) - "inertial"  parameter C
! valPHYS(:,3) - domain speed along the X axis
! valPHYS(:,4) - domain speed along the Y axis
! valPHYS(:,5) - domain speed along the Z axis
REAL(8), ALLOCATABLE  :: valPHYS(:,:) 

CHARACTER (LEN=6), ALLOCATABLE  :: typPHYS(:), &  ! array of symbolic domain types
                                   namePHYS(:)    ! array of domain names

END MODULE m_vxc2data
