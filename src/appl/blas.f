c
c     takes the sum of the absolute values.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision function dasum(n,dx,incx)
      double precision dx(*),dtemp
      integer i,incx,m,mp1,n,nincx
c
      dasum = 0.0d0
      dtemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return

      if(incx.ne.1) then
c     
c     code for increment not equal to 1
c     
	 nincx = n*incx
	 do 10 i = 1,nincx,incx
	    dtemp = dtemp + dabs(dx(i))
 10	 continue
	 dasum = dtemp
	 
      else
c
c     code for increment equal to 1
c
c
c     clean-up loop
c
	 m = mod(n,6)
	 if( m .ne. 0 ) then
	    do 30 i = 1,m
	       dtemp = dtemp + dabs(dx(i))
 30	    continue
	    if( n .lt. 6 ) go to 60
	 endif
	 mp1 = m + 1
	 do 50 i = mp1,n,6
	    dtemp = dtemp + dabs(dx(i))	   + dabs(dx(i + 1)) +
     &			    dabs(dx(i + 2))+ dabs(dx(i + 3)) +
     &			    dabs(dx(i + 4))+ dabs(dx(i + 5))
 50	 continue
 60	 dasum = dtemp
      endif
      return
      end


c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      subroutine daxpy(n,da,dx,incx,dy,incy)
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c	 code for unequal increments or equal increments
c	   not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
	dy(iy) = dy(iy) + da*dx(ix)
	ix = ix + incx
	iy = iy + incy
   10 continue
      return
c
c	 code for both increments equal to 1
c
c
c	 clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
	dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
	dy(i) = dy(i) + da*dx(i)
	dy(i + 1) = dy(i + 1) + da*dx(i + 1)
	dy(i + 2) = dy(i + 2) + da*dx(i + 2)
	dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end


c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      subroutine  dcopy(n,dx,incx,dy,incy)
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c	 code for unequal increments or equal increments
c	   not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
	dy(iy) = dx(ix)
	ix = ix + incx
	iy = iy + incy
   10 continue
      return
c
c	 code for both increments equal to 1
c
c
c	 clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
	dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
	dy(i) = dx(i)
	dy(i + 1) = dx(i + 1)
	dy(i + 2) = dx(i + 2)
	dy(i + 3) = dx(i + 3)
	dy(i + 4) = dx(i + 4)
	dy(i + 5) = dx(i + 5)
	dy(i + 6) = dx(i + 6)
   50 continue
      return
      end


c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision function ddot(n,dx,incx,dy,incy)
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c	 code for unequal increments or equal increments
c	   not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
	dtemp = dtemp + dx(ix)*dy(iy)
	ix = ix + incx
	iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c	 code for both increments equal to 1
c
c	 clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
	dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
	dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *	 dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end


c
c     smach computes machine parameters of floating point
c     arithmetic for use in testing only.  not required by
c     linpack proper.
c
c     if trouble with automatic computation of these quantities,
c     they can be set by direct assignment statements.
c     assume the computer has
c
c	 b = base of arithmetic
c	 t = number of base  b	digits
c	 l = smallest possible exponent
c	 u = largest possible exponent
c
c     then
c
c	 eps = b**(1-t)
c	 tiny = 100.0*b**(-l+t)
c	 huge = 0.01*b**(u-t)
c
c     dmach same as smach except t, l, u apply to
c     double precision.
c
c     cmach same as smach except if complex division
c     is done by
c
c	 1/(x+i*y) = (x-i*y)/(x**2+y**2)
c
c     then
c
c	 tiny = sqrt(tiny)
c	 huge = sqrt(huge)
c
c
c     job is 1, 2 or 3 for epsilon, tiny and huge, respectively.
c
      double precision function dmach(job)
      integer job
      double precision eps,tiny,huge,s
c
      eps = 1.0d0
   10 eps = eps/2.0d0
      s = 1.0d0 + eps
      if (s .gt. 1.0d0) go to 10
      eps = 2.0d0*eps
c
      s = 1.0d0
   20 tiny = s
      s = s/16.0d0
      if (s*1.0 .ne. 0.0d0) go to 20
      tiny = (tiny/eps)*100.0
      huge = 1.0d0/tiny
c
      dmach = 0.0d0
      if (job .eq. 1) dmach = eps
      if (job .eq. 2) dmach = tiny
      if (job .eq. 3) dmach = huge
      return
      end


*
*  dnrm2() returns the euclidean norm of a vector via the function
*  name, so that
*
*     dnrm2 := sqrt( x'*x )
*
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to DLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*
      double precision function dnrm2 ( n, x, incx )
*     .. scalar arguments ..
      integer				incx, n
*     .. array arguments ..
      double precision			x( * )
*     .. parameters ..
      double precision	    one		, zero
      parameter		  ( one = 1.0d+0, zero = 0.0d+0 )
*     .. local scalars ..
      integer		    ix
      double precision	    absxi, norm, scale, ssq
*     .. intrinsic functions ..
      intrinsic		    abs, sqrt
*     ..
*     .. executable statements ..
      if( n.lt.1 .or. incx.lt.1 )then
	 norm  = zero
      else if( n.eq.1 )then
	 norm  = abs( x( 1 ) )
      else
	 scale = zero
	 ssq   = one
*	 the following loop is equivalent to this call to the lapack
*	 auxiliary routine:
*	 call dlassq( n, x, incx, scale, ssq )
*
	 do 10, ix = 1, 1 + ( n - 1 )*incx, incx
	    if( x( ix ).ne.zero )then
	       absxi = abs( x( ix ) )
	       if( scale.lt.absxi )then
		  ssq	= one	+ ssq*( scale/absxi )**2
		  scale = absxi
	       else
		  ssq	= ssq	+     ( absxi/scale )**2
	       end if
	    end if
   10	 continue
	 norm  = scale * sqrt( ssq )
      end if
*
      dnrm2 = norm
      return
      end
*     end of dnrm2.


c
c     applies a plane rotation.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      subroutine  drot (n,dx,incx,dy,incy,c,s)
      double precision dx(*),dy(*),dtemp,c,s
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c	code for unequal increments or equal increments not equal
c	  to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
	dtemp = c*dx(ix) + s*dy(iy)
	dy(iy) = c*dy(iy) - s*dx(ix)
	dx(ix) = dtemp
	ix = ix + incx
	iy = iy + incy
   10 continue
      return
c
c	code for both increments equal to 1
c
   20 do 30 i = 1,n
	dtemp = c*dx(i) + s*dy(i)
	dy(i) = c*dy(i) - s*dx(i)
	dx(i) = dtemp
   30 continue
      return
      end


c
c     construct givens plane rotation.
c     jack dongarra, linpack, 3/11/78.
c
      subroutine drotg(da,db,c,s)
      double precision da,db,c,s,roe,scale,r,z
c
      roe = db
      if( dabs(da) .gt. dabs(db) ) roe = da
      scale = dabs(da) + dabs(db)
      if( scale .ne. 0.0d0 ) go to 10
	 c = 1.0d0
	 s = 0.0d0
	 r = 0.0d0
	 z = 0.0d0
	 go to 20
   10 r = scale*dsqrt((da/scale)**2 + (db/scale)**2)
      r = dsign(1.0d0,roe)*r
      c = da/r
      s = db/r
      z = 1.0d0
      if( dabs(da) .gt. dabs(db) ) z = s
      if( dabs(db) .ge. dabs(da) .and. c .ne. 0.0d0 ) z = 1.0d0/c
   20 da = r
      db = z
      return
      end


c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      subroutine  dscal(n,da,dx,incx)
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c	 code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
	dx(i) = da*dx(i)
   10 continue
      return
c
c	 code for increment equal to 1
c
c
c	 clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
	dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
	dx(i) = da*dx(i)
	dx(i + 1) = da*dx(i + 1)
	dx(i + 2) = da*dx(i + 2)
	dx(i + 3) = da*dx(i + 3)
	dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end


c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      subroutine  dswap (n,dx,incx,dy,incy)
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c	code for unequal increments or equal increments not equal
c	  to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
	dtemp = dx(ix)
	dx(ix) = dy(iy)
	dy(iy) = dtemp
	ix = ix + incx
	iy = iy + incy
   10 continue
      return
c
c	code for both increments equal to 1
c
c	clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
	dtemp = dx(i)
	dx(i) = dy(i)
	dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
	dtemp = dx(i)
	dx(i) = dy(i)
	dy(i) = dtemp
	dtemp = dx(i + 1)
	dx(i + 1) = dy(i + 1)
	dy(i + 1) = dtemp
	dtemp = dx(i + 2)
	dx(i + 2) = dy(i + 2)
	dy(i + 2) = dtemp
   50 continue
      return
      end


c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      integer function idamax(n,dx,incx)
      double precision dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c	 code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
	 if(dabs(dx(ix)).gt.dmax)then
	    idamax = i
	    dmax = dabs(dx(ix))
	 endif
	 ix = ix + incx
   10 continue
      return
c
c	 code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
	 if(dabs(dx(i)).le.dmax) go to 30
	 idamax = i
	 dmax = dabs(dx(i))
   30 continue
      return
      end
