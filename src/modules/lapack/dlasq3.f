      SUBROUTINE DLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL,
     $                   ITER, NDIV )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            I0, ITER, N0, NDIV, NFAIL, PP
      DOUBLE PRECISION   DESIG, DMIN, QMAX, SIGMA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
*     ..
*
*  Purpose
*  =======
*  DLASQ3 checks for deflation, computes a shift (TAU) and calls dqds.
*  In case of failure it changes shifts, and tries again until output
*  is positive.
*
*  Arguments
*  =========
*
*  I0     (input) INTEGER
*         First index.
*
*  N0     (input) INTEGER
*         Last index.
*
*  Z      (input) DOUBLE PRECISION array, dimension ( 4*N )
*         Z holds the qd array.
*
*  PP     (input) INTEGER
*         PP=0 for ping, PP=1 for pong.
*
*  DMIN   (output) DOUBLE PRECISION
*         Minimum value of d.
*
*  SIGMA  (output) DOUBLE PRECISION
*         Sum of shifts used in current segment.
*
*  DESIG  (input/output) DOUBLE PRECISION
*         Lower order part of SIGMA
*
*  QMAX   (input) DOUBLE PRECISION
*         Maximum value of q.
*
*  NFAIL  (output) INTEGER
*         Number of times shift was too big.
*
*  ITER   (output) INTEGER
*         Number of iterations.
*
*  NDIV   (output) INTEGER
*         Number of divisions.
*
*  TTYPE  (output) INTEGER
*         Shift type.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   CBIAS
      PARAMETER          ( CBIAS = 1.50D0 )
      DOUBLE PRECISION   ZERO, QURTR, HALF, ONE, TWO, TEN
      PARAMETER          ( ZERO = 0.0D0, QURTR = 0.250D0, HALF = 0.5D0,
     $                   ONE = 1.0D0, TWO = 2.0D0, TEN = 10.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IPN4, J4, N0IN, NN, TTYPE
      DOUBLE PRECISION   DMIN1, DMIN2, DN, DN1, DN2, EPS, EPS2, S,
     $                   SFMIN, T, TAU, TEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASQ4, DLASQ5, DLASQ6
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Save statement ..
      SAVE               TTYPE, DMIN1, DMIN2, DN, DN1, DN2, TAU
*     ..
*     .. Data statements ..
      DATA               TTYPE / 0 /
      DATA               DMIN1 / ZERO / , DMIN2 / ZERO / , DN / ZERO / ,
     $                   DN1 / ZERO / , DN2 / ZERO / , TAU / ZERO /
*     ..
*     .. Executable Statements ..
*
      N0IN = N0
      EPS = DLAMCH( 'Precision' )*TEN
      SFMIN = DLAMCH( 'Safe minimum' )
      EPS2 = EPS**2
*
*     Check for deflation.
*
   10 CONTINUE
*
      IF( N0.LT.I0 )
     $   RETURN
      IF( N0.EQ.I0 )
     $   GO TO 20
      NN = 4*N0 + PP
      IF( N0.EQ.( I0+1 ) )
     $   GO TO 40
*
*     Check whether E(N0-1) is negligible, 1-by-1 case.
*
      IF( Z( NN-5 ).GT.EPS2*( SIGMA+Z( NN-3 ) ) .AND. Z( NN-2*PP-4 ).GT.
     $    EPS2*Z( NN-7 ) )GO TO 30
*
   20 CONTINUE
*
      Z( 4*N0-3 ) = Z( 4*N0+PP-3 ) + SIGMA
      N0 = N0 - 1
      GO TO 10
*
*     Check  whether E(N0-2) is negligible, 2-by-2 case.
*
   30 CONTINUE
*
      IF( Z( NN-9 ).GT.EPS2*SIGMA .AND. Z( NN-2*PP-8 ).GT.EPS2*
     $    Z( NN-11 ) )GO TO 50
*
   40 CONTINUE
*
      IF( Z( NN-3 ).GT.Z( NN-7 ) ) THEN
         S = Z( NN-3 )
         Z( NN-3 ) = Z( NN-7 )
         Z( NN-7 ) = S
      END IF
      IF( Z( NN-5 ).GT.Z( NN-3 )*EPS2 ) THEN
         T = HALF*( ( Z( NN-7 )-Z( NN-3 ) )+Z( NN-5 ) )
         S = Z( NN-3 )*( Z( NN-5 ) / T )
         IF( S.LE.T ) THEN
            S = Z( NN-3 )*( Z( NN-5 ) / ( T*( ONE+SQRT( ONE+S /
     $          T ) ) ) )
         ELSE
            S = Z( NN-3 )*( Z( NN-5 ) / ( T+SQRT( T )*SQRT( T+S ) ) )
         END IF
         T = Z( NN-7 ) + ( S+Z( NN-5 ) )
         Z( NN-3 ) = Z( NN-3 )*( Z( NN-7 ) / T )
         Z( NN-7 ) = T
      END IF
      Z( 4*N0-7 ) = Z( NN-7 ) + SIGMA
      Z( 4*N0-3 ) = Z( NN-3 ) + SIGMA
      N0 = N0 - 2
      GO TO 10
*
   50 CONTINUE
*
*     Reverse the qd-array, if warranted.
*
      IF( DMIN.LE.ZERO .OR. N0.LT.N0IN ) THEN
         IF( CBIAS*Z( 4*I0+PP-3 ).LT.Z( 4*N0+PP-3 ) ) THEN
            IPN4 = 4*( I0+N0 )
            DO 60 J4 = 4*I0, 2*( I0+N0-1 ), 4
               TEMP = Z( J4-3 )
               Z( J4-3 ) = Z( IPN4-J4-3 )
               Z( IPN4-J4-3 ) = TEMP
               TEMP = Z( J4-2 )
               Z( J4-2 ) = Z( IPN4-J4-2 )
               Z( IPN4-J4-2 ) = TEMP
               TEMP = Z( J4-1 )
               Z( J4-1 ) = Z( IPN4-J4-5 )
               Z( IPN4-J4-5 ) = TEMP
               TEMP = Z( J4 )
               Z( J4 ) = Z( IPN4-J4-4 )
               Z( IPN4-J4-4 ) = TEMP
   60       CONTINUE
            IF( N0-I0.LE.4 ) THEN
               Z( 4*N0+PP-1 ) = Z( 4*I0+PP-1 )
               Z( 4*N0-PP ) = Z( 4*I0-PP )
            END IF
            DMIN2 = MIN( DMIN2, Z( 4*N0+PP-1 ) )
            Z( 4*N0+PP-1 ) = MIN( Z( 4*N0+PP-1 ), Z( 4*I0+PP-1 ),
     $                       Z( 4*I0+PP+3 ) )
            Z( 4*I0-PP ) = MIN( Z( 4*N0-PP ), Z( 4*I0-PP ),
     $                     Z( 4*I0-PP+4 ) )
            QMAX = MAX( QMAX, Z( 4*I0+PP-3 ), Z( 4*I0+PP+1 ) )
            DMIN = -ZERO
         END IF
      END IF
*
   70 CONTINUE
*
      IF( DMIN.LT.ZERO .OR. SFMIN*QMAX.LE.
     $    MIN( Z( 4*N0+PP-1 ), Z( 4*N0+PP-9 ), DMIN2+Z( 4*N0-PP ) ) )
     $     THEN
*
*        Choose a shift.
*
         CALL DLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, DN1,
     $                DN2, TAU, TTYPE )
*
*        Call dqds until DMIN > 0.
*
   80    CONTINUE
*
         CALL DLASQ5( I0, N0, Z, PP, TAU, DMIN, DMIN1, DMIN2, DN, DN1,
     $                DN2 )
*
         ITER = ITER + 1
         NDIV = NDIV + ( N0-I0+2 )
*
*        Check for NaN: "DMIN.NE.DMIN"
*
         IF( DMIN.NE.DMIN ) THEN
            Z( 4*N0+PP-1 ) = ZERO
            TAU = ZERO
            GO TO 70
         END IF
*
*        Check for convergence hidden by negative DN.
*
         IF( DMIN.LT.ZERO .AND. DMIN1.GT.ZERO .AND.
     $       Z( 4*( N0-1 )-PP ).LT.EPS*( SIGMA+DN1 ) .AND. ABS( DN ).LT.
     $       EPS*SIGMA ) THEN
            Z( 4*( N0-1 )-PP+2 ) = ZERO
            DMIN = ABS( DMIN )
         END IF
*
         IF( DMIN.LT.ZERO ) THEN
*
*           Failure. Select new TAU and try again.
*
            NFAIL = NFAIL + 1
*
*           Failed twice. Play it safe.
*
            IF( TTYPE.LT.-22 ) THEN
               TAU = ZERO
               GO TO 80
            END IF
*
            IF( DMIN1.GT.ZERO ) THEN
*
*              Late failure. Gives excellent shift.
*
               TAU = ( TAU+DMIN )*( ONE-TWO*EPS )
               TTYPE = TTYPE - 11
            ELSE
*
*              Early failure. Divide by 4.
*
               TAU = QURTR*TAU
               TTYPE = TTYPE - 12
            END IF
            GO TO 80
         END IF
      ELSE
         CALL DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, DN1, DN2 )
         ITER = ITER + 1
         NDIV = NDIV + ( N0-I0 )
         TAU = ZERO
      END IF
*
      IF( TAU.LT.SIGMA ) THEN
         DESIG = DESIG + TAU
         T = SIGMA + DESIG
         DESIG = DESIG - ( T-SIGMA )
      ELSE
         T = SIGMA + TAU
         DESIG = SIGMA - ( T-TAU ) + DESIG
      END IF
      SIGMA = T
*
      RETURN
*
*     End of DLASQ3
*
      END
