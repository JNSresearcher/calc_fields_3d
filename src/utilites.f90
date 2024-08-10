! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )

SUBROUTINE writeVtk_src ( Npoint, numfun,  Ncells0, NcellsX, NcellsY, NcellsZ, new_nodes0, new_nodesX, new_nodesY,  new_nodesZ,  &
                          num_nod0, num_nodV,  sdx, sdy, sdz, delta, files)
                          
! This program generates a .vtk file for sources (scalar or vector, moving or stationary)

USE m_vxc2data
IMPLICIT NONE

INTEGER,      INTENT(IN):: Npoint, numfun, num_nod0, num_nodV, sdx,sdy,sdz 
INTEGER,      INTENT(IN):: Ncells0, NcellsX, NcellsY, NcellsZ, &
                           new_nodes0(Ncells0),new_nodesX(NcellsX),new_nodesY(NcellsY),new_nodesZ(NcellsZ) !(nsub)
REAL(8),      INTENT(IN):: delta(3)
CHARACTER(*), INTENT(IN):: files

CHARACTER (LEN = 1)  ci
CHARACTER (LEN = 4)  ch_sd
CHARACTER (LEN = 30) buf1, buf2
INTEGER:: i,j,k, l, m, n,numcells, nn0, nnX, nnY, nnZ, nij,  ios, n_new !nim,njm,nkm,nip,njp,nkp,
REAL(8):: s,sxm,sym,szm

numcells = Ncells0 + NcellsX + NcellsY + NcellsZ

WRITE(ch_sd,'(i4)')Npoint
ios=0

! for ifort
! open(newunit=n_new,file='out/src_'//trim(files)//trim(adjustl(ch_sd))//'.vtk',form="UNFORMATTED", buffered='yes', blocksize=524288, buffercount=1,  &
! access="STREAM",  convert="big_endian", iostat=ios)  

OPEN(newunit=n_new,file='out/src_'//trim(files)//trim(adjustl(ch_sd))//'.vtk',form="UNFORMATTED",&
access="STREAM",  convert="big_endian", iostat=ios)   

IF (ios/=0) THEN
    PRINT*,'error! Could not open the outfile result.';   STOP
ENDIF

! header output
WRITE (n_new   ) "# vtk DataFile Version 3.0"//new_line(ci)//"out data result"//new_line(ci)//"BINARY" //new_line(ci)
WRITE(n_new  )"DATASET UNSTRUCTURED_GRID"//new_line(ci)
WRITE(buf1,'(i8)')   numcells*8  !
WRITE(n_new  ) "POINTS "//trim(adjustl(buf1))//" double" //new_line(ci)

! output of coordinates
IF (num_nod0 /=0) THEN   ! сначала скал источники
    nn0 = 0              ! счетчик узлов скал источ
    DO n=1,numfun
        IF (Fun(n)%ex == '0') THEN
            DO l = 1, fun_nod(n)%numnod_F0
                nn0 = nn0 + 1
                m = new_nodes0(nn0)
            
                CALL find_coord
                CALL write_coord
            ENDDO
        ENDIF
    ENDDO
ENDIF ! if num_nod0

IF (num_nodV /=0) THEN ! затем вект источники
    nnX = 0; nnY = 0;nnZ = 0;   ! счетчик узлов вект источ
    DO n=1,numfun
        IF (Fun(n)%ex == 'X') THEN
            DO l = 1, fun_nod(n)%numnod_Fx
                nnX = nnX + 1
                m = new_nodesX(nnX)
                CALL find_coord
                CALL write_coord
            ENDDO
        ELSEIF (Fun(n)%ex == 'Y') THEN
            DO l = 1,fun_nod(n)%numnod_Fy
                nnY = nnY + 1
                m = new_nodesY(nnY)
                CALL find_coord
                CALL write_coord
            ENDDO
        ELSEIF (Fun(n)%ex == 'Z') THEN
            DO l = 1,fun_nod(n)%numnod_Fz
                nnZ = nnZ + 1
                m = new_nodesZ(nnZ)
                CALL find_coord
                CALL write_coord
            ENDDO
        ENDIF
    ENDDO
ENDIF ! if num_nodV

WRITE(n_new)  new_line(ci)

WRITE(buf1,'(i8)')   numcells 
WRITE(buf2,'(i8)') 9*numcells  
WRITE(n_new) "CELLS "//trim(adjustl(buf1))//" "//trim(adjustl(buf2))//new_line(ci)
DO i=0, numcells - 1   !(Ncells0 + NcellsX + NcellsY + NcellsZ)-1
    WRITE(n_new) 8, 8*i + [0, 1, 2, 3, 4, 5, 6, 7]
ENDDO
WRITE(n_new)   new_line(ci)

WRITE(n_new ) "CELL_TYPES "//trim(adjustl(buf1))//new_line(ci)
DO i=1,numcells  !(Ncells0 + NcellsX + NcellsY + NcellsZ)
    WRITE(n_new) 11
ENDDO
WRITE(n_new)  new_line(ci)

WRITE(n_new) "CELL_DATA "//trim(adjustl(buf1))//new_line(ci)

IF (num_nod0 /=0) THEN 
!         scalar field output
    WRITE(n_new)"SCALARS "//'Scalar_field '//"double" //" 1"//new_line(ci)//"LOOKUP_TABLE default"//new_line(ci)
    DO n=1,numfun
        s = Fun(n)%vely
        IF (Fun(n)%ex == '0') WRITE(n_new) (s, k = 1, fun_nod(n)%numnod_F0)
    ENDDO
ENDIF

IF (num_nodV /=0) THEN 
!          vector field output : grad(U0)
    WRITE(n_new) "VECTORS "//"Vector_field_SRC"//" double"//new_line(ci)
    DO n=1,numfun
        s = Fun(n)%vely
        IF     (Fun(n)%ex == 'X') THEN
            sxm = s; sym=0.d0;; szm=0.d0; 
            WRITE(n_new) ( [sxm, sym, szm], k = 1, fun_nod(n)%numnod_Fx)
        ELSEIF (Fun(n)%ex == 'Y') THEN
            sxm = 0.d0; sym=s;; szm=0.d0; 
            WRITE(n_new) ( [sxm, sym, szm], k = 1, fun_nod(n)%numnod_Fy)
        ELSEIF (Fun(n)%ex == 'Z') THEN
            sxm = 0.d0; sym=0.d0; szm = s
            WRITE(n_new) ( [sxm, sym, szm], k = 1, fun_nod(n)%numnod_Fz)
        ENDIF
    ENDDO
ENDIF

WRITE(n_new)  new_line(ci)
   
CLOSE(n_new)

CONTAINS

    SUBROUTINE find_coord
        k = ceiling(REAL(m)/( REAL(sdx*sdy) ) ) 
        IF (k == 1) THEN
            nij = m
        ELSE
            nij = m - (k-1)*sdx*sdy
        ENDIF
        j = ceiling( REAL(nij) / REAL(sdx) ) 
        i = nij - (j - 1) * sdx  
    END SUBROUTINE find_coord

    SUBROUTINE write_coord
        ! point 0 x = i; y=j ; z=k
        sxm = REAL(i,8)*delta(1) - delta(1)
        sym = REAL(j,8)*delta(2) - delta(2) 
        szm = REAL(k,8)*delta(3) - delta(3) 
        WRITE(n_new ) sxm, sym, szm
        
        ! point 1 x = i+1; y=j ; z=k
        sxm = REAL(i+1,8)*delta(1) - delta(1)
        sym = REAL(j,8)*delta(2) - delta(2) 
        szm = REAL(k,8)*delta(3) - delta(3) 
        WRITE(n_new) sxm, sym, szm
        
        ! point 2 x = i; y=j+1 ; z=k
        sxm = REAL(i,8)*delta(1) - delta(1)
        sym = REAL(j+1,8)*delta(2) - delta(2) 
        szm = REAL(k,8)*delta(3) - delta(3) 
        WRITE(n_new) sxm, sym, szm
        
        ! point 3 x = i+1; y=j+1 ; z=k
        sxm = REAL(i+1,8)*delta(1) - delta(1)
        sym = REAL(j+1,8)*delta(2) - delta(2) 
        szm = REAL(k,8)*delta(3) - delta(3) 
        WRITE(n_new) sxm, sym, szm
        
        ! point 4 x = i; y=j ; z=k+1
        sxm = REAL(i,8)*delta(1) - delta(1)
        sym = REAL(j,8)*delta(2) - delta(2) 
        szm = REAL(k+1,8)*delta(3) - delta(3) 
        WRITE(n_new ) sxm, sym, szm
        
        ! point 5 x = i+1; y=j ; z=k+1
        sxm = REAL(i+1,8)*delta(1) - delta(1)
        sym = REAL(j,8)*delta(2) - delta(2) 
        szm = REAL(k+1,8)*delta(3) - delta(3) 
        WRITE(n_new ) sxm, sym, szm
        
        ! point 6 x = i; y=j+1 ; z=k+1
        sxm = REAL(i,8)*delta(1) - delta(1)
        sym = REAL(j+1,8)*delta(2) - delta(2) 
        szm = REAL(k+1,8)*delta(3) - delta(3) 
        WRITE(n_new) sxm, sym, szm
        
        ! point 7 x = i+1; y=j+1 ; z=k+1
        sxm = REAL(i+1,8)*delta(1) - delta(1)
        sym = REAL(j+1,8)*delta(2) - delta(2) 
        szm = REAL(k+1,8)*delta(3) - delta(3) 
        WRITE(n_new) sxm, sym, szm

    END SUBROUTINE write_coord
END SUBROUTINE writeVtk_src


SUBROUTINE writeVtk_field (Npoint, num_nod0, num_nodV, sdx, sdy, sdz, delta, V0, SRC0, Ux,Uy,Uz, SRCx,SRCy,SRCz, geoPHYS, typPHYS, size_PHYS_C, files)
! This program generates a .vtk file for fields 

IMPLICIT NONE

INTEGER,      INTENT(IN):: Npoint, num_nod0, num_nodV, sdx,sdy,sdz, size_PHYS_C
REAL(8),      INTENT(IN):: V0(*), SRC0(*), Ux(*),Uy(*),Uz(*), SRCx(*),SRCy(*),SRCz(*),  delta(3)
CHARACTER(*), INTENT(IN):: files
INTEGER(1),   INTENT(IN):: geoPHYS(sdx,sdy,sdz)
CHARACTER(6), INTENT(IN):: typPHYS(*)

CHARACTER (LEN = 1)  ci
CHARACTER (LEN = 4)  ch_sd
CHARACTER (LEN = 30) buf1
INTEGER:: i,j,k,m, n,kdz,  nim,njm,nkm,nip,njp,nkp, ios, n_new
REAL(8):: s,sxm,sym,szm
    
kdz=sdx*sdy

WRITE(ch_sd,'(i4)')Npoint
ios=0

! for ifort
! OPEN(newunit=n_new,file='out/Field_'//trim(files)//trim(adjustl(ch_sd))//'.vtk',form="UNFORMATTED", buffered='yes', blocksize=524288, buffercount=1,  &
! access="STREAM",  convert="big_endian", iostat=ios)  

! for gfortran
OPEN(newunit=n_new,file='out/Field_'//trim(files)//trim(adjustl(ch_sd))//'.vtk',form="UNFORMATTED",  &
access="STREAM",  convert="big_endian", iostat=ios)  

IF (ios/=0) THEN
    print *,'error! Could not open the outfile result.';   stop
ENDIF

! header output
WRITE (n_new) "# vtk DataFile Version 3.0"//new_line(ci)//"out data result"//new_line(ci)//"BINARY"//new_line(ci)
WRITE(buf1,'(i8," ",i8," ",i8)') sdx, sdy, sdz  
WRITE(n_new)"DATASET STRUCTURED_GRID"//new_line(ci)//"DIMENSIONS "//trim(adjustl(buf1))//new_line(ci)
WRITE(buf1,'(i8)')   (sdx) * (sdy) * (sdz)
WRITE(n_new) "POINTS "//trim(adjustl(buf1))//" float"//new_line(ci)

! output of coordinates
DO k=1,sdz
    szm = REAL(k,8)*delta(3) - delta(3) 
    DO j=1,sdy
        sym = REAL(j,8)*delta(2) - delta(2) 
        DO i=1,sdx
            sxm = REAL(i,8)*delta(1) - delta(1) 
            WRITE(n_new) REAL(sxm,4), REAL(sym,4), REAL(szm,4)
        ENDDO
    ENDDO
ENDDO
WRITE(n_new)  new_line(ci)

WRITE(n_new) "POINT_DATA "//trim(adjustl(buf1))//new_line(ci)

IF (num_nod0 /=0) THEN
   !        scalar field output
    WRITE(n_new)"SCALARS "//'Scalar_field_U0 '//"float" //" 1"//new_line(ci)//"LOOKUP_TABLE default"//new_line(ci)
    m=0
    DO k = 1,sdz;  DO j = 1,sdy;  DO i = 1,sdx
        m=m+1
        WRITE(n_new) REAL( V0(m) )
    ENDDO;  ENDDO;  ENDDO;
    WRITE(n_new)  new_line(ci)

    !     vector field output : grad(U0)
    WRITE(n_new) "VECTORS "//"Vector_field_Grad(U0)"//" float"//new_line(ci)
    DO k = 1,sdz;  DO j = 1,sdy;  DO i = 1,sdx;
        m   = i   + sdx*(j-1)  + kdz*(k-1)
        nim = m-1;  njm = m-sdx;  nkm = m-kdz;
        nip = m+1;  njp = m+sdx;  nkp = m+kdz; 
        IF (k==1 .or.k==sdz .or.j==1 .or.j==sdy .or.i==1 .or.i==sdx  )  THEN !.or. SRC0(m) /=0
            sxm=0.d0; sym=0.d0; szm=0.d0
        ELSE
            sxm = -0.5d0*(V0(nip) - V0(nim)) / delta(1)   ! Ux = dU0/dx
            sym = -0.5d0*(V0(njp) - V0(njm)) / delta(2)   ! Uy = dU0/dy 
            szm = -0.5d0*(V0(nkp) - V0(nkm)) / delta(3)   ! Uz = dU0/dz
        ENDIF
        WRITE(n_new) REAL(sxm,4), REAL(sym,4), REAL(szm,4)
    ENDDO;  ENDDO;  ENDDO
    WRITE(n_new)  new_line(ci)
   
    !     scalar field output
    WRITE(n_new)"SCALARS "//'Scalar_field_SOURCE '//"float" //" 1"//new_line(ci)//"LOOKUP_TABLE default"//new_line(ci)
    m=0
    DO k = 1,sdz;  DO j = 1,sdy;  DO i = 1,sdx
        m=m+1
        WRITE(n_new) REAL( SRC0(m) )
    ENDDO;  ENDDO;  ENDDO
    WRITE(n_new)  new_line(ci)
ENDIF !num_nod0

IF (num_nodV /=0) THEN
    !  vector field output :  curl(U)
    WRITE(n_new) "VECTORS "//"Vector_field_curl(U)"//" float"//new_line(ci)

    m=0
    DO k=1,sdz               ! dUz/dy - dUy/dz = curl_x(V)
    DO j=1,sdy               ! dUx/dz - dUz/dx = curl_y(V)
    DO i=1,sdx               ! dUy/dx - dUx/dy = curl_z(V)
        m=m+1
        nim = m - 1;  njm = m - sdx;  nkm = m - kdz;
        nip = m + 1;  njp = m + sdx;  nkp = m + kdz; 
        IF (i==1) nim = m; IF (i==sdx) nip  =m; 
        IF (j==1) njm = m; IF (j==sdy) njp = m; 
        IF (k==1) nkm = m; IF (k==sdz) nkp = m; 
        sxm = 0.5d0*( ( Uz(njp) - Uz(njm)) /delta(2) - ( Uy(nkp) - Uy(nkm) ) /delta(3) )
        sym = 0.5d0*( ( Ux(nkp) - Ux(nkm)) /delta(3) - ( Uz(nip) - Uz(nim) ) /delta(1) )
        szm = 0.5d0*( ( Uy(nip) - Uy(nim)) /delta(1) - ( Ux(njp) - Ux(njm) ) /delta(2) )
        WRITE(n_new) REAL(sxm,4), REAL(sym,4), REAL(szm,4)
    ENDDO
    ENDDO
    ENDDO
    WRITE(n_new)  new_line(ci)

    WRITE(n_new) "VECTORS "//"Vector_field_U"//" float"//new_line(ci)
    m = 0
    DO k = 1,sdz;  DO j = 1,sdy;  DO i = 1,sdx
        m = m + 1
        WRITE(n_new) REAL( Ux(m),4 ), REAL( Uy(m),4  ), REAL(Uz(m),4)
    ENDDO;  ENDDO;  ENDDO
    WRITE(n_new)  new_line(ci)

    IF (size_PHYS_C /= 0 ) THEN
        ! vector field output :  eddy current
        WRITE(n_new) "VECTORS "//"Vector_field_eddy"//" float"//new_line(ci)
        m=0; s=0.d0
        DO k = 1,sdz;  DO j = 1,sdy;  DO i = 1,sdx
            m = m+1
            n = geoPHYS(i,j,k)
            IF (index(typPHYS(n), 'C' ) /= 0) THEN ! is domain with C
                WRITE(n_new) REAL( -SRCx(m),4 ), REAL( -SRCy(m),4  ), REAL( -SRCz(m),4 )
            ELSE
                WRITE(n_new) REAL( s,4 ), REAL( s,4  ), REAL( s,4 )
            ENDIF
        ENDDO;  ENDDO;  ENDDO
        WRITE(n_new)  new_line(ci)
    
        ! vector field output :  source
        WRITE(n_new) "VECTORS "//"Vector_field_SOURCE"//" float"//new_line(ci)
        m=0; s=0.d0
        DO k = 1,sdz;  DO j = 1,sdy;  DO i = 1,sdx
            m = m+1
            n = geoPHYS(i,j,k)
            IF (index(typPHYS(n), 'C' ) == 0) THEN ! is domain with C
                WRITE(n_new) REAL( SRCx(m),4 ), REAL( SRCy(m),4  ), REAL( SRCz(m),4 )
            ELSE
                WRITE(n_new) REAL( s,4 ), REAL( s,4  ), REAL( s,4 )
            ENDIF
        ENDDO;  ENDDO;  ENDDO
        WRITE(n_new)  new_line(ci)
    ELSE
         ! vector field output :  only source 
        WRITE(n_new) "VECTORS "//"Vector_field_SOURCE"//" float"//new_line(ci)
        m=0
        DO k = 1,sdz;  DO j = 1,sdy;  DO i = 1,sdx
            m = m+1
            WRITE(n_new) REAL( SRCx(m),4 ), REAL( SRCy(m),4  ), REAL( SRCz(m),4 )
        ENDDO;  ENDDO;  ENDDO
        WRITE(n_new)  new_line(ci)
    ENDIF 
ENDIF !num_nodV /=0

CLOSE(n_new)
END SUBROUTINE writeVtk_field

SUBROUTINE string2words(Ligne,words,num_words)
! converting string to array of words

CHARACTER (LEN = *), INTENT(IN) :: Ligne
CHARACTER (LEN = *), INTENT(OUT):: words(*) 
INTEGER Start
INTEGER Finish
INTEGER num_words
INTEGER length 

    length= len_trim(Ligne)
    Start = 0;Finish=0;num_words=0
    DO j = 1, length
        IF (Ligne(j:j) == ' ' .or. Ligne(j:j)== '	') THEN
            IF (Start > 0) THEN
                num_words = num_words + 1
                words(num_words) = Ligne(Start : Finish)
                Start = 0
            ENDIF
        ELSE IF (Start == 0) THEN
            Start = j;  Finish = j
        ELSE
            Finish = Finish + 1
        ENDIF
    ENDDO
    IF (Start > 0) THEN
        num_words = num_words + 1
        words(num_words) = Ligne(Start : Finish)
        Start = 0
    ENDIF

END SUBROUTINE string2words


SUBROUTINE UPP (st1, st2)
! converting string to upper case

    CHARACTER (LEN=*),  INTENT(in) :: st1
    CHARACTER (LEN=*), INTENT(out) :: st2
    CHARACTER (LEN=*),   PARAMETER :: lsym = 'abcdefghijklmnopqrstuvwxyz'
    CHARACTER (LEN=*),   PARAMETER :: usym = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    st2 = st1
    DO j=1,LEN_TRIM(st1)
        k = INDEX(lsym,st1(j:j))
        IF (k > 0) st2(j:j) = usym(k:k)
    END DO
END SUBROUTINE UPP

REAL(8) function numeric(SA)
! Converting a string to a number using decimal prefixes.
! (Spice-like format used in electrical circuit simulation programs).

! Accepted decimal prefix formats in this program  :
! 'f',   femto    1e-15
! 'p',   pico     1e-12
! 'n',   nano     1e-9
! 'u',   micro    1e-6
! 'm',   milli    1e-3
! 'c'    centi    1e-2
! 'h'    hecto    1e+2
! 'k',   kilo     1e+3
! 'meg'  mega     1e+6
! 'g',   giga     1e+9
! 't',   tera     1e+12
! 'pet'  peta     1e+15

! some features:
! - Instead of a decimal point, on can use prefix: 1.3k = 1k3
! - can put exponential format at the end 
! - spaces inside are ignored
! - prefixes are automatiCALLy converted to uppercase during conversion

! Fortran code created by J.Sochor  ( https://github.com/JNSresearcher/convert_prefix ) 

IMPLICIT NONE

    CHARACTER (LEN=*),   PARAMETER :: lsym = 'abcdefghijklmnopqrstuvwxyz'
    CHARACTER (LEN=*),   PARAMETER :: usym = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!                                            1   2   3    4    5   6   7    8    9   10
    CHARACTER*1 :: ch, prefixes(10)=(/'M','K','U','N','P','G','T','F','C','H'/) 
    CHARACTER (LEN = *) sa
    CHARACTER*20 inter
    REAL(8) multiplier   ! factor
    INTEGER i,j,k,l,m,n,lm,lp
    
    DO j=1,LEN_TRIM(SA)
        k = INDEX(lsym,SA(j:j))
        IF (k > 0) SA(j:j) = usym(k:k)
    END DO
    
    multiplier = 1.0_8;lp=0;lm=0
    m=index(SA,',')
    IF (m /= 0) SA(m:m)='.'

    DO i=1,9
        l=index(SA,prefixes(i))
        IF (l /= 0) THEN
            m=1
            IF ( i == 1 ) THEN
                lm=index(SA, 'MEG')
                IF (lm /= 0) THEN
                    multiplier=1.E+06_8
                ELSE
                    multiplier=1.E-3_8
                ENDIF
                exit
            ELSEIF ( i == 2 ) THEN
                multiplier=1.E+03_8
                exit
            ELSEIF ( i == 3 ) THEN
                multiplier=1.E-06_8
                exit
            ELSEIF ( i == 4 ) THEN
                multiplier=1.E-09_8
                exit
            ELSEIF ( i == 5 ) THEN
                lm=index(SA, 'PET')
                IF (lm /= 0) THEN
                    multiplier=1.E+15_8
                ELSE
                    multiplier=1.E-12_8
                ENDIF
                exit
            ELSEIF ( i == 6 ) THEN
                multiplier=1.E+09_8
                exit
            ELSEIF ( i == 7 ) THEN
                multiplier=1.E+12_8
                exit
            ELSEIF ( i == 8 ) THEN
                multiplier=1.E-15_8
                exit
            ELSEIF ( i == 9 ) THEN
                multiplier=1.E-2_8
                exit
            ELSEIF ( i == 10 ) THEN
                multiplier=1.E+2_8
                exit
            ENDIF
        ENDIF
    ENDDO

    k=index(SA,'.');j=len_trim(SA)
    IF ( (m /= 0).AND.(k==0) )  THEN
        SA(l:l)='.'
        IF (lm /= 0) THEN
            SA(lm+1:lm+2)='  '
            sa=sa(1:lm)//sa(lm+3:)
        ENDIF
    ELSE
        IF (lm /= 0) THEN
            SA(lm:lm+2)='   '
            sa=sa(1:lm)//sa(lm+3:)
        ENDIF
    ENDIF
    
    j=len_trim(SA)
    i=index(SA,'E')
    IF (i == 0) THEN
        DO i=1,j
            ch=SA(i:i)
            n=ichar(ch)
            IF ((n.GT.57).or.(n.LT.48)) THEN
                IF     (ch == '.') THEN
                    cycle
                ELSEIF (ch == '-') THEN
                    cycle
                ELSE
                    SA(i:i) = ' '
                ENDIF
            ELSE
            ENDIF
        ENDDO
    ENDIF

    WRITE (inter,'(A20)') SA(1:j)
    !print*,SA(1:j)
    READ (inter,'(BN,G20.0)') numeric
    numeric = multiplier * numeric

END FUNCTION numeric


