      SUBROUTINE ZLAQPS( M, N, OFFSET, NB, KB, A, LDA, JPVT, TAU, VN1,
     $                   VN2, AUXV, F, LDF )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            KB, LDA, LDF, M, N, NB, OFFSET
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      DOUBLE PRECISION   VN1( * ), VN2( * )
      COMPLEX*16         A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLAQPS computes a step of QR factorization with column pivoting
*  of a complex M-by-N matrix A by using Blas-3.  It tries to factorize
*  NB columns from A starting from the row OFFSET+1, and updates all
*  of the matrix with Blas-3 xGEMM.
*
*  In some cases, due to catastrophic cancellations, it cannot
*  factorize NB columns.  Hence, the actual number of factorized
*  columns is returned in KB.
*
*  Block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A. N >= 0
*
*  OFFSET  (input) INTEGER
*          The number of rows of A that have been factorized in
*          previous steps.
*
*  NB      (input) INTEGER
*          The number of columns to factorize.
*
*  KB      (output) INTEGER
*          The number of columns actually factorized.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, block A(OFFSET+1:M,1:KB) is the triangular
*          factor obtained and block A(1:OFFSET,1:N) has been
*          accordingly pivoted, but no factorized.
*          The rest of the matrix, block A(OFFSET+1:M,KB+1:N) has
*          been updated.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  JPVT    (input/output) INTEGER array, dimension (N)
*          JPVT(I) = K <==> Column K of the full matrix A has been
*          permuted into position I in AP.
*
*  TAU     (output) COMPLEX*16 array, dimension (KB)
*          The scalar factors of the elementary reflectors.
*
*  VN1     (input/output) DOUBLE PRECISION array, dimension (N)
*          The vector with the partial column norms.
*
*  VN2     (input/output) DOUBLE PRECISION array, dimension (N)
*          The vector with the exact column norms.
*
*  AUXV    (input/output) COMPLEX*16 array, dimension (NB)
*          Auxiliar vector.
*
*  F       (input/output) COMPLEX*16 array, dimension (LDF,NB)
*          Matrix F' = L*Y'*A.
*
*  LDF     (input) INTEGER
*          The leading dimension of the array F. LDF >= max(1,N).
*
*  Further Details
*  ===============
*
*  Based on contributions by
*    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
*    X. Sun, Computer Science Dept., Duke University, USA
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0,
     $                   CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            ITEMP, J, K, LASTRK, LSTICC, PVT, RK
      DOUBLE PRECISION   TEMP, TEMP2
      COMPLEX*16         AKK
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMM, ZGEMV, ZLARFG, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, MAX, MIN, NINT, SQRT
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DZNRM2
      EXTERNAL           IDAMAX, DZNRM2
*     ..
*     .. Executable Statements ..
*
      LASTRK = MIN( M, N+OFFSET )
      LSTICC = 0
      K = 0
*
*     Beginning of while loop.
*
   10 CONTINUE
      IF( ( K.LT.NB ) .AND. ( LSTICC.EQ.0 ) ) THEN
         K = K + 1
         RK = OFFSET + K
*
*        Determine ith pivot column and swap if necessary
*
         PVT = ( K-1 ) + IDAMAX( N-K+1, VN1( K ), 1 )
         IF( PVT.NE.K ) THEN
            CALL ZSWAP( M, A( 1, PVT ), 1, A( 1, K ), 1 )
            CALL ZSWAP( K-1, F( PVT, 1 ), LDF, F( K, 1 ), LDF )
            ITEMP = JPVT( PVT )
            JPVT( PVT ) = JPVT( K )
            JPVT( K ) = ITEMP
            VN1( PVT ) = VN1( K )
            VN2( PVT ) = VN2( K )
         END IF
*
*        Apply previous Householder reflectors to column K:
*        A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)'.
*
         IF( K.GT.1 ) THEN
*CC            CALL ZGEMM( 'No transpose', 'Conjugate transpose',
*CC     $                  M-RK+1, 1, K-1, -CONE, A( RK, 1 ), LDA,
*CC     $                  F( K, 1 ), LDF, CONE, A( RK, K ), LDA )
            DO 20 J = 1, K - 1
               F( K, J ) = DCONJG( F( K, J ) )
   20       CONTINUE
            CALL ZGEMV( 'No transpose', M-RK+1, K-1, -CONE, A( RK, 1 ),
     $                  LDA, F( K, 1 ), LDF, CONE, A( RK, K ), 1 )
            DO 30 J = 1, K - 1
               F( K, J ) = DCONJG( F( K, J ) )
   30       CONTINUE
         END IF
*
*        Generate elementary reflector H(k).
*
         IF( RK.LT.M ) THEN
            CALL ZLARFG( M-RK+1, A( RK, K ), A( RK+1, K ), 1, TAU( K ) )
         ELSE
            CALL ZLARFG( 1, A( RK, K ), A( RK, K ), 1, TAU( K ) )
         END IF
*
         AKK = A( RK, K )
         A( RK, K ) = CONE
*
*        Compute Kth column of F:
*
*        Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)'*A(RK:M,K).
*
         IF( K.LT.N ) THEN
            CALL ZGEMV( 'Conjugate transpose', M-RK+1, N-K, TAU( K ),
     $                  A( RK, K+1 ), LDA, A( RK, K ), 1, CZERO,
     $                  F( K+1, K ), 1 )
         END IF
*
*        Padding F(1:K,K) with zeros.
*
         DO 40 J = 1, K
            F( J, K ) = CZERO
   40    CONTINUE
*
*        Incremental updating of F:
*        F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)'
*                    *A(RK:M,K).
*
         IF( K.GT.1 ) THEN
            CALL ZGEMV( 'Conjugate transpose', M-RK+1, K-1, -TAU( K ),
     $                  A( RK, 1 ), LDA, A( RK, K ), 1, CZERO,
     $                  AUXV( 1 ), 1 )
*
            CALL ZGEMV( 'No transpose', N, K-1, CONE, F( 1, 1 ), LDF,
     $                  AUXV( 1 ), 1, CONE, F( 1, K ), 1 )
         END IF
*
*        Update the current row of A:
*        A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)'.
*
         IF( K.LT.N ) THEN
            CALL ZGEMM( 'No transpose', 'Conjugate transpose', 1, N-K,
     $                  K, -CONE, A( RK, 1 ), LDA, F( K+1, 1 ), LDF,
     $                  CONE, A( RK, K+1 ), LDA )
         END IF
*
*        Update partial column norms.
*
         IF( RK.LT.LASTRK ) THEN
            DO 50 J = K + 1, N
               IF( VN1( J ).NE.ZERO ) THEN
                  TEMP = ABS( A( RK, J ) ) / VN1( J )
                  TEMP = MAX( ZERO, ( ONE+TEMP )*( ONE-TEMP ) )
                  TEMP2 = ONE + 0.05D0*TEMP*( VN1( J ) / VN2( J ) )**2
                  IF( TEMP2.EQ.ONE ) THEN
                     VN2( J ) = DBLE( LSTICC )
                     LSTICC = J
                  ELSE
                     VN1( J ) = VN1( J )*SQRT( TEMP )
                  END IF
               END IF
   50       CONTINUE
         END IF
*
         A( RK, K ) = AKK
*
*        End of while loop.
*
         GO TO 10
      END IF
      KB = K
      RK = OFFSET + KB
*
*     Apply the block reflector to the rest of the matrix:
*     A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) -
*                         A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)'.
*
      IF( KB.LT.MIN( N, M-OFFSET ) ) THEN
         CALL ZGEMM( 'No transpose', 'Conjugate transpose', M-RK, N-KB,
     $               KB, -CONE, A( RK+1, 1 ), LDA, F( KB+1, 1 ), LDF,
     $               CONE, A( RK+1, KB+1 ), LDA )
      END IF
*
*     Recomputation of difficult columns.
*
   60 CONTINUE
      IF( LSTICC.GT.0 ) THEN
         ITEMP = NINT( VN2( LSTICC ) )
         VN1( LSTICC ) = DZNRM2( M-RK, A( RK+1, LSTICC ), 1 )
         VN2( LSTICC ) = VN1( LSTICC )
         LSTICC = ITEMP
         GO TO 60
      END IF
*
      RETURN
*
*     End of ZLAQPS
*
      END
