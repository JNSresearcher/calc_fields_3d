! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )

SUBROUTINE SOR (sdx,sdy,sdz, geoPHYS, geoUbnd, Q,W,Qmx,Qmy,Qmz,Qpx,Qpy,Qpz, rhs, X, tolerance, itmax )
! Successive OverRelaxation https://en.wikipedia.org/wiki/Successive_over-relaxation

IMPLICIT NONE

INTEGER,   INTENT(IN)    :: itmax, sdx,sdy,sdz
REAL(8),   INTENT(IN)    :: rhs(sdx*sdy*sdz), tolerance
REAL(8),   INTENT(IN)    :: Q(*),Qmx(*),Qmy(*),Qmz(*),Qpx(*),Qpy(*),Qpz(*)
REAL(8),   INTENT(IN)    :: W
INTEGER(1),INTENT(IN)    :: geoPHYS(sdx,sdy,sdz), geoUbnd(sdx,sdy,sdz)
REAL(8),   INTENT(INOUT) :: X(sdx*sdy*sdz)
REAL(8) RMIN, RESID, a
INTEGER i,j,k, L, m, n,  kdz, nim, njm, nkm, nip, njp, nkp

    kdz=sdx*sdy
    RMIN=1.1d0*tolerance; 
    L=0
    DO WHILE( RMIN >= tolerance  .and. L < itmax )
        L = L + 1
        RMIN = 0.d0;
 
        DO k=2,sdz-1;  DO j=2,sdy-1;  DO i=2,sdx-1
            IF (geoUbnd(i,j,k) == 0 ) THEN
                m   = i   + sdx*(j-1)   + sdx*sdy*(k-1)
              
                n=geoPHYS(i,j,k);
                a = W/Q(m);

                nim = m-1;  njm = m-sdx;  nkm = m-kdz;
                nip = m+1;  njp = m+sdx;  nkp = m+kdz;  
    
                RESID =  a*(   ( Qpx(m)* X(nip) + Qmx(m)*X(nim)  )      &
                            +  ( Qpy(m)* X(njp) + Qmy(m)*X(njm)  )      &
                            +  ( Qpz(m)* X(nkp) + Qmz(m)*X(nkm)  )      &
                            -    Q(m)  * X(m)   + rhs(m)    )  
                RMIN = RMIN + ABS(RESID)
                X(m) = X(m) + RESID
            ENDIF    
        ENDDO;   ENDDO;  ENDDO
    ENDDO
    
    IF (L < itmax) THEN
    ELSE
        PRINT '(a,i5)', 'No iter ',L
        ! write(*,'((a),I6,a,es14.7,a,es14.7,a,es14.7)') 'Iter=', L, ' err_x=',RMINx, ' err_y=',RMINy, ' err_z=',RMINz
    ENDIF

END SUBROUTINE SOR

SUBROUTINE sprsBCGstabWR (valA, irow, jcol, n, b, x, tolerance, itmax, iter)
! Biconjugate Gradient Stabilized with rewrite method 
! This is the version if the right side (vector B) is not 0
! discussed in detail in my project https://github.com/JNSresearcher/SOLVERS_BCGSTAB_GMRES

INTEGER jcol(*), irow(n+1)
REAL(8) valA(*)
INTEGER iter, itmax,  n 
REAL(8):: alpha, beta, omega, tolerance, rr0, rr0_new, Bnorm, &
                    R(n), R0(n), P(n), AP(n), S(n), AS(n), B(n), X(n)

    iter=0 
    R = sprsAx(X)
    DO j=1,n
        R(j) = B(j) - R(j)
    ENDDO 
    R0 = R
    P = R

    Bnorm=norm2(b)
    IF (Bnorm == 0.d0) RETURN
    DO
        IF (iter > itmax) THEN
            PRINT*, norm2(R)
            EXIT
        ENDIF
        iter = iter + 1
        AP = sprsAx(P)
        rr0 = DOT_PRODUCT (R,R0)
        alpha = rr0/dot_product(AP,R0) 
        S = R - alpha*AP
        IF (norm2(S) < tolerance) THEN
            X = X + alpha*P;
            EXIT
        ENDIF
        AS = sprsAx(S)
        omega = dot_product(AS,S)/dot_product(AS,AS)
        X = X + alpha*P + omega*S
        R = S - omega*AS
        IF ( (norm2(R)/Bnorm) < tolerance) EXIT
        rr0_new = dot_product(R,R0)
        beta = (alpha/omega)*rr0_new/rr0
        P = R + beta*(P - omega*AP)
        IF ( (abs(rr0_new)/Bnorm) < tolerance) THEN
            R0 = R; P = R
        ENDIF
    ENDDO

CONTAINS

FUNCTION sprsAx  (V) 
    REAL(8)  :: sprsAx(n), V(n)
    INTEGER i1,i2
    DO  i=1,n
        i1=irow(i); i2=irow(i+1)-1
        sprsAx(i) =  dot_product(valA(i1:i2), V(jcol(i1:i2) ) )
    ENDDO 
END FUNCTION sprsAx 

END SUBROUTINE sprsBCGstabWR

SUBROUTINE sprsBCGstabWR2 (valA, irow, jcol, n, b, x, tolerance, itmax, iter)
! Biconjugate Gradient Stabilized with rewrite method. 
! This is the version if the right side (vector B) is 0
! discussed in detail in my project https://github.com/JNSresearcher/SOLVERS_BCGSTAB_GMRES

INTEGER jcol(*), irow(n+1)
REAL(8) valA(*)
INTEGER iter, itmax,  n 
REAL(8):: alpha, beta, omega, tolerance, rr0, rr0_new, &  !  Bnorm,
                    R(n), R0(n), P(n), AP(n), S(n), AS(n), B(n), X(n)
    iter=0 
    R = sprsAx(X)
    DO j=1,n
        R(j) = B(j) - R(j)
    ENDDO 
    R0 = R
    P = R
    DO
        IF (iter > itmax) THEN
            PRINT*, norm2(R)
            EXIT
        ENDIF
        iter = iter + 1
        AP = sprsAx(P)
        rr0 = DOT_PRODUCT (R,R0)
        alpha = rr0/dot_product(AP,R0) 
        S = R - alpha*AP
        IF (norm2(S) < tolerance) THEN
            X = X + alpha*P;
            EXIT
        ENDIF
        AS = sprsAx(S)
        omega = dot_product(AS,S)/dot_product(AS,AS)
        X = X + alpha*P + omega*S
        R = S - omega*AS
        IF ( norm2(R) < tolerance) EXIT
        rr0_new = dot_product(R,R0)
        beta = (alpha/omega)*rr0_new/rr0
        P = R + beta*(P - omega*AP)
        IF (abs(rr0_new) < tolerance) THEN
            R0 = R; P = R
        ENDIF
    ENDDO

CONTAINS

FUNCTION sprsAx  (V) 
    REAL(8)  :: sprsAx(n), V(n)
    INTEGER i1,i2
    DO  i=1,n
        i1=irow(i); i2=irow(i+1)-1
        sprsAx(i) =  dot_product(valA(i1:i2), V(jcol(i1:i2) ) )
    ENDDO 
END FUNCTION sprsAx 

END SUBROUTINE sprsBCGstabWR2