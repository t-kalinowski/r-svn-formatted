      SUBROUTINE ZUNGR2( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZUNGR2 generates an m by n complex matrix Q with orthonormal rows,
*  which is defined as the last m rows of a product of k elementary
*  reflectors of order n
*
*        Q  =  H(1)' H(2)' . . . H(k)'
*
*  as returned by ZGERQF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. N >= M.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. M >= K >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the (m-k+i)-th row must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by ZGERQF in the last k rows of its array argument
*          A.
*          On exit, the m-by-n matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZGERQF.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (M)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, II, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLACGV, ZLARF, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGR2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 )
     $   RETURN
*
      IF( K.LT.M ) THEN
*
*        Initialise rows 1:m-k to rows of the unit matrix
*
         DO 20 J = 1, N
            DO 10 L = 1, M - K
               A( L, J ) = ZERO
   10       CONTINUE
            IF( J.GT.N-M .AND. J.LE.N-K )
     $         A( M-N+J, J ) = ONE
   20    CONTINUE
      END IF
*
      DO 40 I = 1, K
         II = M - K + I
*
*        Apply H(i)' to A(1:m-k+i,1:n-k+i) from the right
*
         CALL ZLACGV( N-M+II-1, A( II, 1 ), LDA )
         A( II, N-M+II ) = ONE
         CALL ZLARF( 'Right', II-1, N-M+II, A( II, 1 ), LDA,
     $               DCONJG( TAU( I ) ), A, LDA, WORK )
         CALL ZSCAL( N-M+II-1, -TAU( I ), A( II, 1 ), LDA )
         CALL ZLACGV( N-M+II-1, A( II, 1 ), LDA )
         A( II, N-M+II ) = ONE - DCONJG( TAU( I ) )
*
*        Set A(m-k+i,n-k+i+1:n) to zero
*
         DO 30 L = N - M + II + 1, N
            A( II, L ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of ZUNGR2
*
      END
