C
C     Modified from the SMART package by J.H. Friedman, 10/10/84
C     Main change is to add spline smoothing modified from BRUTO,
C     calling code written for smooth.spline in S.
C
C     B.D. Ripley (ripley@stats.ox.ac.uk)  1994-7.
C
C
      subroutine BDRsmart(m,mu,p,q,n,w,x,y,ww,smod,nsmod,
     & sp,nsp,dp,ndp,edf)
      implicit double precision (a-h, o-z)
      integer p,q
      double precision x(p,n),y(q,n),w(n),ww(q),smod(nsmod),
     & sp(nsp),edf(m)
      double precision dp(ndp)
      smod(1)=m
      smod(2)=p
      smod(3)=q
      smod(4)=n
      call BDRsmart1(m,mu,p,q,n,w,x,y,ww,smod(6),smod(q+6),smod(q+7),
     & smod(q+7+p*m),smod(q+7+m*(p+q)),smod(q+7+m*(p+q+n)),  
     & smod(q+7+m*(p+q+2*n)),sp,sp(q*n+1),sp(n*(q+15)+1),  
     & sp(n*(q+15)+q+1),dp,smod(5),edf)
      return
      end

      subroutine BDRsmart1(m,mu,p,q,n,w,x,y,ww,yb,ys,a,b,f,t,
     & asr,r,sc,bt,g,dp,flm,edf)
      implicit double precision (a-h, o-z)
      integer p,q
      double precision w(n),x(p,n),y(1),ww(q),a(p,m),f(n,m),t(n,m),
     & r(q,n)
      double precision sc(n,15),g(p,3),yb(q),b(q,m),bt(q),edf(m)
      double precision dp(1)
      common /BDRparms/ ifl,lf,span,alpha,big
      common /BDRz01/ conv,maxit,mitone,cutmin,fdel,cjeps,mitcj
      sw=0.0
      do 10161 j=1,n
10161   sw=sw+w(j)
      do 10201 j=1,n
        do 10201 i=1,q
10201     r(i,j)=y(q*(j-1)+i)
      do 10241 i=1,q
        s=0.0
        do 10251 j=1,n
10251     s=s+w(j)*r(i,j)
        yb(i)=s/sw
10241   continue
c     yb is vector of means
      do 10261 j=1,n
        do 10261 i=1,q
10261     r(i,j)=r(i,j)-yb(i)
      ys=0.
      do 10281 i=1,q
        s=0.
        do 10291 j=1,n
10291     s=s+w(j)*r(i,j)**2
10281   ys=ys+ww(i)*s/sw
      if(ys .gt. 0.0) goto 10311
c     ys is the overall standard deviation -- quit if zero
      return
10311 continue
      ys=sqrt(ys)
      s=1./ys
      do 10331 j=1,n
        do 10331 i=1,q
10331     r(i,j)=r(i,j)*s
c     r is now standardized residuals
c     subfit adds up to m  terms one at time; lm is the number fitted.
      call BDRsubfit(m,p,q,n,w,sw,x,r,ww,lm,a,b,f,t,asr,sc,bt,g,dp,edf)
      if(lf.le.0) go to 10350
      call BDRfulfit(lm,lf,p,q,n,w,sw,x,r,ww,a,b,f,t,asr,sc,bt,g,dp,edf)
10371 continue
      do 10381 l=1,lm
        sc(l,1)=l+.1
        s=0.0
        do 10391 i=1,q
10391     s=s+ww(i)*abs(b(i,l))
10381   sc(l,2)=-s
      call BDRsort(sc(1,2),sc,1,lm)
      do 10461 j=1,n
        do 10461 i=1,q
10461     r(i,j)=y(q*(j-1)+i)
      do 10521 i=1,q
        do 10531 j=1,n
          r(i,j)=r(i,j)-yb(i)
          s=0.
          do 10541 l=1,lm
10541       s=s+b(i,l)*f(j,l)
          r(i,j)=r(i,j)/ys-s
10531     continue
10521   continue
      if(lm.le.mu) goto 10350
      l=sc(lm,1)
      asr=0.0
      do 10561 j=1,n
        do 10561 i=1,q
          r(i,j)=r(i,j)+b(i,l)*f(j,l)
10561     asr=asr+w(j)*ww(i)*r(i,j)**2
      asr=asr/sw
      if(l .ge. lm) goto 10591
      do 10601 i=1,p
10601   a(i,l)=a(i,lm)
      do 10611 i=1,q
10611   b(i,l)=b(i,lm)
      do 10621 j=1,n
        f(j,l)=f(j,lm)
10621   t(j,l)=t(j,lm)
10591 continue
      lm=lm-1
      call BDRfulfit(lm,lf,p,q,n,w,sw,x,r,ww,a,b,f,t,asr,sc,bt,g,dp,edf)
      goto 10371
10350 continue
      flm=lm
      return
      end

      subroutine BDRsubfit(m,p,q,n,w,sw,x,r,ww,lm,a,b,f,t,asr,sc,
     & bt,g,dp,edf)
      implicit double precision (a-h, o-z)
      integer p,q
      double precision w(n),x(p,n),ww(q),a(p,m),f(n,m),t(n,m),r(q,n)
      double precision sc(n,15),g(p,3),b(q,m),bt(q),edf(m)
      double precision dp(1)
      common /BDRparms/ ifl,lf,span,alpha,big
      common /BDRz01/ conv,maxit,mitone,cutmin,fdel,cjeps,mitcj
      asr=big
      lm=0
      do 10641 l=1,m
        lm=lm+1
        asrold=asr
        call BDRnewb(lm,q,ww,b)
        call BDRonetrm(0,p,q,n,w,sw,x,r,ww,a(1,lm),b(1,lm),
     &    f(1,lm),t(1,lm),asr,sc,g,dp,edf)
        do 10651 j=1,n
          do 10651 i=1,q
10651       r(i,j)=r(i,j)-b(i,lm)*f(j,lm)
      if(lm.eq.1) goto 10641
      if(lf.le.0) goto 10690
      if(lm.eq.m) return
      iflsv=ifl
      ifl=0
      call BDRfulfit(lm,1,p,q,n,w,sw,x,r,ww,a,b,f,t,asr,sc,bt,g,dp,edf)
      ifl=iflsv
10690 continue
      if(asr.le.0.0.or.(asrold-asr)/asrold.lt.conv) return
10641 continue
      return
      end

      subroutine BDRfulfit(lm,lbf,p,q,n,w,sw,x,r,ww,a,b,f,t,
     & asr,sc,bt,g,dp,edf)
      implicit double precision (a-h, o-z)
      integer p,q
      double precision w(n),x(p,n),ww(q),a(p,lm),f(n,lm),t(n,lm),r(q,n)
      double precision sc(n,15),g(p,3),b(q,lm),bt(q), asr(lm), edf(lm)
      double precision dp(1)
      common /BDRparms/ ifl,lf,span,alpha,big
      common /BDRz01/ conv,maxit,mitone,cutmin,fdel,cjeps,mitcj
      if(lbf.le.0) return
      asri=asr(1)
      iter=0
      fsv=cutmin
      isv=mitone
      if(lbf .ge. 3) goto 10711
      cutmin=1.0
      mitone=lbf-1
10711 continue
      asrold=asri
      iter=iter+1
      do 10731 lp=1,lm
        do 10741 i=1,q
10741     bt(i)=b(i,lp)
        do 10751 i=1,p
10751     g(i,3)=a(i,lp)
        do 10761 j=1,n
          do 10761 i=1,q
10761       r(i,j)=r(i,j)+bt(i)*f(j,lp)
        call BDRonetrm(1,p,q,n,w,sw,x,r,ww,g(1,3),bt,sc(1,14),sc(1,15),
     &    asri,sc,g,dp,edf(lp))
        if(asri .lt. asrold) then
        do 10801 i=1,q
10801     b(i,lp)=bt(i)
        do 10811 i=1,p
10811     a(i,lp)=g(i,3)
        do 10821 j=1,n
          f(j,lp)=sc(j,14)
10821     t(j,lp)=sc(j,15)
        else
           asri=asrold
        endif
        do 10841 j=1,n
          do 10841 i=1,q
10841       r(i,j)=r(i,j)-b(i,lp)*f(j,lp)
10731   continue
      if((iter .le. maxit) .and. ((asri .gt. 0.0) .and. 
     & ((asrold-asri)/asrold .ge. conv))) goto 10711
      cutmin=fsv
      mitone=isv
      if(ifl .le. 0) goto 10871
      asr(1+lm) = asri
      asr(1) = asri
10871 continue
      return
      end

      subroutine BDRonetrm(jfl,p,q,n,w,sw,x,y,ww,a,b,f,t,asr,
     & sc,g,dp,edf)
      implicit double precision (a-h, o-z)
      integer p,q
      double precision w(n),x(p,n),y(q,n),ww(q),a(p),b(q),f(n),t(n)
      double precision sc(n,13),g(p,2)
      double precision dp(1)
      common /BDRparms/ ifl,lf,span,alpha,big
      common /BDRz01/ conv,maxit,mitone,cutmin,fdel,cjeps,mitcj
      iter=0
      asr=big
10901 continue
      iter=iter+1
      asrold=asr
      do 10911 j=1,n
        s=0.0
        do 10921 i=1,q
10921     s=s+ww(i)*b(i)*y(i,j)
10911   sc(j,13)=s
      call BDRoneone(max0(jfl,iter-1),p,n,w,sw,sc(1,13),x,a,f,t,
     & asr,sc,g,dp,edf)
      do 10931 i=1,q
        s=0.0
        do 10941 j=1,n
10941     s=s+w(j)*y(i,j)*f(j)
10931   b(i)=s/sw
      asr=0.0
      do 10951 i=1,q
        s=0.0
        do 10961 j=1,n
10961     s=s+w(j)*(y(i,j)-b(i)*f(j))**2
10951   asr=asr+ww(i)*s/sw
      if((q .ne. 1) .and. ((iter .le. maxit) .and. ((asr .gt. 0.0) 
     & .and. ((asrold-asr)/asrold .ge. conv)))) goto 10901
      return
      end

      subroutine BDRoneone(ist,p,n,w,sw,y,x,a,f,t,asr,sc,g,dp,edf)
      implicit double precision (a-h, o-z)
      integer p
      double precision w(n),y(n),x(p,n),a(p),f(n),t(n),sc(n,12),g(p,2)
      double precision dp(1)
      common /BDRparms/ ifl,lf,span,alpha,big
      common /BDRz01/ conv,maxit,mitone,cutmin,fdel,cjeps,mitcj
      sml=1.0/big
      if(ist .gt. 0) goto 10981
      if(p .gt. 1) goto 11001
      a(1)=1.0
      goto 11011
11001 continue
      do 11021 j=1,n
        sc(j,2)=1.0
11021   continue
      call BDRdir(p,n,w,sw,y,x,sc(1,2),a,dp)
11011 continue
10981 continue
      s=0.0
      do 11031 i=1,p
        g(i,1)=0.0
        s=s+a(i)**2
11031   continue
      s=1.0/sqrt(s)
      do 11041 i=1,p
        a(i)=a(i)*s
11041   continue
      iter=0
      asr=big
      cut=1.0
11051 continue
      iter=iter+1
      asrold=asr
11060 continue
      s=0.0
      do 11071 i=1,p
        g(i,2)=a(i)+g(i,1)
        s=s+g(i,2)**2
11071 continue
      s=1./sqrt(s)
      do 11081 i=1,p
        g(i,2)=g(i,2)*s
11081   continue
      do 11091 j=1,n
        sc(j,1)=j+.1
        s=0.
        do 11101 i=1,p
          s=s+g(i,2)*x(i,j)
11101     continue
        sc(j,11)=s
11091   continue
      call BDRsort(sc(1,11),sc,1,n)
      do 11111 j=1,n
        k=sc(j,1)
        sc(j,2)=y(k)
        sc(j,3)=max(w(k),sml)
11111   continue
      call BDRsupsmu(n,sc(1,11),sc(1,2),sc(1,3),1,span,alpha,
     & sc(1,12),sc(1,4), edf)
      s=0.0
      do 11121 j=1,n
        s=s+sc(j,3)*(sc(j,2)-sc(j,12))**2
11121   continue
      s=s/sw
      if(s .lt. asr) goto 11141
      cut=cut*0.5
      if(cut.lt.cutmin) goto 11052
      do 11151 i=1,p
        g(i,1)=g(i,1)*cut
11151   continue
      go to 11060
11141 continue
      asr=s
      cut=1.0
      do 11161 i=1,p
        a(i)=g(i,2)
11161   continue
      do 11171 j=1,n
        k=sc(j,1)
        t(k)=sc(j,11)
        f(k)=sc(j,12)
11171   continue
      if(asr.le.0.0.or.(asrold-asr)/asrold.lt.conv) goto 11052
      if(iter.gt.mitone.or.p.le.1) goto 11052
      call BDRder(n,sc(1,11),sc(1,12),sc(1,3),fdel,sc(1,4),sc(1,5))
      do 11181 j=1,n
        k=sc(j,1)
        sc(j,5)=y(j)-f(j)
        sc(k,6)=sc(j,4)
11181   continue
      call BDRdir(p,n,w,sw,sc(1,5),x,sc(1,6),g,dp)
      goto 11051
11052 continue
      s=0.0
      v=s
      do 11191 j=1,n
        s=s+w(j)*f(j)
11191   continue
      s=s/sw
      do 11201 j=1,n
        f(j)=f(j)-s
        v=v+w(j)*f(j)**2
11201   continue
      if(v .le. 0.0) goto 11221
      v=1.0/sqrt(v/sw)
      do 11231 j=1,n
11231   f(j)=f(j)*v
11221 continue
      return
      end


      subroutine BDRdir(p,n,w,sw,r,x,d,e,g)
      implicit double precision (a-h, o-z)
      integer p
      double precision w(n),r(n),x(p,n),d(n),e(p)
      double precision s,g(1)
      common /BDRz01/ conv,maxit,mitone,cutmin,fdel,cjeps,mitcj
      do 11241 i=1,p
      s=0.0
      do 11251 j=1,n
      s=s+w(j)*d(j)*x(i,j)
11251 continue
      e(i)=s/sw
11241 continue
      k=0
      m1=p*(p+1)/2
      m2=m1+p
      do 11261 j=1,p
      s=0.0
      do 11271 l=1,n
      s=s+w(l)*r(l)*(d(l)*x(j,l)-e(j))
11271 continue
      g(m1+j)=s/sw
      do 11281 i=1,j
      s=0.0
      do 11291 l=1,n
      s=s+w(l)*(d(l)*x(i,l)-e(i))*(d(l)*x(j,l)-e(j))
11291 continue
      k=k+1
      g(k)=s/sw
11281 continue
11261 continue
      call BDRconj(p,g,g(m1+1),g(m2+1),cjeps,mitcj,g(m2+p+1))
      do 11301 i=1,p
      e(i)=g(m2+i)
11301 continue
      return
      end

      subroutine BDRconj(p,g,c,x,eps,maxit,sc)
      implicit double precision (a-h, o-z)
      integer p
      double precision g(1),c(p),x(p),sc(p,4)
      double precision beta,h,s,alpha,t
      do 11311 i=1,p
      x(i)=0.0
      sc(i,2)=0.0
11311 continue
      nit=0
11321 continue
      nit=nit+1
      h=0.0
      beta=0.0
      do 11331 i=1,p
      sc(i,4)=x(i)
      s=g(i*(i-1)/2+i)*x(i)
      im1=i-1
      j=1
      goto 11343
11341 j=j+1
11343 if((j).gt.(im1)) goto 11342
      s=s+g(i*(i-1)/2+j)*x(j)
      goto 11341
11342 continue
      j=i+1
      goto 11353
11351 j=j+1
11353 if((j).gt.(p)) goto 11352
      s=s+g(j*(j-1)/2+i)*x(j)
      goto 11351
11352 continue
      sc(i,1)=s-c(i)
      h=h+sc(i,1)**2
11331 continue
      if(h.le.0.0) goto 11322
      do 11361 iter=1,p
      do 11371 i=1,p
      sc(i,2)=beta*sc(i,2)-sc(i,1)
11371 continue
      t=0.0
      do 11381 i=1,p
      s=g(i*(i-1)/2+i)*sc(i,2)
      im1=i-1
      j=1
      goto 11393
11391 j=j+1
11393 if((j).gt.(im1)) goto 11392
      s=s+g(i*(i-1)/2+j)*sc(j,2)
      goto 11391
11392 continue
      j=i+1
      goto 11403
11401 j=j+1
11403 if((j).gt.(p)) goto 11402
      s=s+g(j*(j-1)/2+i)*sc(j,2)
      goto 11401
11402 continue
      sc(i,3)=s
      t=t+s*sc(i,2)
11381 continue
      alpha=h/t
      s=0.0
      do 11411 i=1,p
      x(i)=x(i)+alpha*sc(i,2)
      sc(i,1)=sc(i,1)+alpha*sc(i,3)
      s=s+sc(i,1)**2
11411 continue
      if(s.le.0.0) goto 11362
      beta=s/h
      h=s
11361 continue
11362 continue
      s=0.0
      do 11421 i=1,p
      s=dmax1(s,dabs(x(i)-sc(i,4)))
11421 continue
      if((s .ge. eps) .and. (nit .lt. maxit)) goto 11321
11322 continue
      return
      end

      subroutine BDRder (n,x,s,w,fdel,d,sc)
      implicit double precision (a-h, o-z)
      double precision x(n),s(n),w(n),d(n),sc(n,3)
      integer bl,el,bc,ec,br,er
      if(x(n) .gt. x(1)) goto 11441
      do 11451 j=1,n
      d(j)=0.0
11451 continue
      return
11441 continue
      i=n/4
      j=3*i
      scale=x(j)-x(i)
11461 if(scale.gt.0.0) goto 11462
      if(j.lt.n) j=j+1
      if(i.gt.1) i=i-1
      scale=x(j)-x(i)
      goto 11461
11462 continue
      del=fdel*scale*2.0
      do 11471 j=1,n
      sc(j,1)=x(j)
      sc(j,2)=s(j)
      sc(j,3)=w(j)
11471 continue
      call BDRpool (n,sc,sc(1,2),sc(1,3),del)
      bc=0
      br=bc
      er=br
11481 continue
      br=er+1
      er=br
11491 if(er .ge. n) goto 11492
      if(sc(br,1) .ne. sc(er+1,1)) goto 11511
      er=er+1
      goto 11521
11511 continue
      goto 11492
11521 continue
11501 continue
      goto 11491
11492 continue
      if(br .ne. 1) goto 11541
      bl=br
      el=er
      goto 11481
11541 continue
      if(bc .ne. 0) goto 11561
      bc=br
      ec=er
      do 11571 j=bl,el
      d(j)=(sc(bc,2)-sc(bl,2))/(sc(bc,1)-sc(bl,1))
11571 continue
      goto 11481
11561 continue
      do 11581 j=bc,ec
      d(j)=(sc(br,2)-sc(bl,2))/(sc(br,1)-sc(bl,1))
11581 continue
      if(er .ne. n) goto 11601
      do 11611 j=br,er
      d(j)=(sc(br,2)-sc(bc,2))/(sc(br,1)-sc(bc,1))
11611 continue
      goto 11482
11601 continue
      bl=bc
      el=ec
      bc=br
      ec=er
      goto 11481
11482 continue
      return
      end

      subroutine BDRpool (n,x,y,w,del)
      implicit double precision (a-h, o-z)
      double precision x(n),y(n),w(n)
      integer bb,eb,br,er,bl,el
      bb=0
      eb=bb
11621 if(eb.ge.n) goto 11622
      bb=eb+1
      eb=bb
11631 if(eb .ge. n) goto 11632
      if(x(bb) .ne. x(eb+1)) goto 11651
      eb=eb+1
      goto 11661
11651 continue
      goto 11632
11661 continue
11641 continue
      goto 11631
11632 continue
      if(eb .ge. n) goto 11681
      if(x(eb+1)-x(eb) .ge. del) goto 11701
      br=eb+1
      er=br
11711 if(er .ge. n) goto 11712
      if(x(er+1) .ne. x(br)) goto 11731
      er=er+1
      goto 11741
11731 continue
      goto 11712
11741 continue
11721 continue
      goto 11711
11712 continue
      if(x(er+1)-x(er).lt.x(eb+1)-x(eb)) goto 11621
      eb=er
      pw=w(bb)+w(eb)
      px=(x(bb)*w(bb)+x(eb)*w(eb))/pw
      py=(y(bb)*w(bb)+y(eb)*w(eb))/pw
      do 11751 i=bb,eb
      x(i)=px
      y(i)=py
      w(i)=pw
11751 continue
11701 continue
11681 continue
11761 continue
      if(bb.le.1) goto 11762
      if(x(bb)-x(bb-1).ge.del) goto 11762
      bl=bb-1
      el=bl
11771 if(bl .le. 1) goto 11772
      if(x(bl-1) .ne. x(el)) goto 11791
      bl=bl-1
      goto 11801
11791 continue
      goto 11772
11801 continue
11781 continue
      goto 11771
11772 continue
      bb=bl
      pw=w(bb)+w(eb)
      px=(x(bb)*w(bb)+x(eb)*w(eb))/pw
      py=(y(bb)*w(bb)+y(eb)*w(eb))/pw
      do 11811 i=bb,eb
      x(i)=px
      y(i)=py
      w(i)=pw
11811 continue
      goto 11761
11762 continue
      goto 11621
11622 continue
      return
      end

      subroutine BDRnewb(lm,q,ww,b)
      implicit double precision (a-h, o-z)
      integer q
      double precision ww(q),b(q,lm)
      common /BDRparms/ ifl,lf,span,alpha,big
      sml=1.0/big
      if(q .ne. 1) goto 11831
      b(1,lm)=1.0
      return
11831 continue
      if(lm .ne. 1) goto 11851
      do 11861 i=1,q
      b(i,lm)=i
11861 continue
      return
11851 continue
      lm1=lm-1
      do 11871 i=1,q
      b(i,lm)=0.0
11871 continue
      t=0.0
      do 11881 i=1,q
      s=0.0
      do 11891 l=1,lm1
      s=s+abs(b(i,l))
11891 continue
      b(i,lm)=s
      t=t+s
11881 continue
      do 11901 i=1,q
      b(i,lm)=ww(i)*(t-b(i,lm))
11901 continue
      l1=1
      if(lm.gt.q) l1=lm-q+1
      do 11911 l=l1,lm1
      s=0.0
      t=s
      do 11921 i=1,q
      s=s+ww(i)*b(i,lm)*b(i,l)
      t=t+ww(i)*b(i,l)**2
11921 continue
      s=s/sqrt(t)
      do 11931 i=1,q
      b(i,lm)=b(i,lm)-s*b(i,l)
11931 continue
11911 continue
      do 11941 i=2,q
      if(abs(b(i-1,lm)-b(i,lm)).gt.sml) return
11941 continue
      do 11951 i=1,q
      b(i,lm)=i
11951 continue
      return
      end

      block data BDRbkppr
      implicit double precision (a-h, o-z)
      common /BDRparms/ ifl,lf,span,alpha,big
      common /BDRz01/ conv,maxit,mitone,cutmin,fdel,cjeps,mitcj
      common /BDRspsmooth/ df, gcvpen, ismethod
      double precision df, gcvpen
      data ifl,maxit,conv,mitone,cutmin,fdel,span,alpha,big,cjeps,
     &mitcj, lf  /6,20,.005,20,0.1,0.02,0.0,0.0,1.0e20,0.001,1,2/
      data df, gcvpen, ismethod /4d0, 1d0, 0/
      end

      subroutine BDRsetppr(span1, alpha1, optlevel, ism, df1, 
     & gcvpen1)
      implicit double precision (a-h, o-z)
      integer optlevel
      double precision alpha1,df, df1, gcvpen, gcvpen1
      common /BDRparms/ ifl,lf,span,alpha,big
      common /BDRspsmooth/ df, gcvpen, ismethod
      span = span1
      lf = optlevel
      alpha = alpha1
      ismethod = ism
      df = df1
      gcvpen = gcvpen1
      end

      subroutine BDRfsort(mu,n,f,t,sp)
      implicit double precision (a-h, o-z)
      double precision f(n,mu),t(n,mu),sp(n,2)
      do 12151 l=1,mu
      do 12161 j=1,n
      sp(j,1)=j+.1
      sp(j,2)=f(j,l)
12161 continue
      call BDRsort(t(1,l),sp,1,n)
      do 12171 j=1,n
      k=sp(j,1)
      f(j,l)=sp(k,2)
12171 continue
12151 continue
      return
      end
      

      subroutine BDRpred(np,x,smod,y,sc)
      implicit double precision (a-h, o-z)
      integer p,q
      double precision x(np,*),y(np,*),smod(*), sc(*)
      integer place,low,high
      m=smod(1)+.1
      p=smod(2)+.1
      q=smod(3)+.1
      n=smod(4)+.1
      mu=smod(5)+.1
      ys=smod(q+6)
      ja=q+6
      jb=ja+p*m
      jf=jb+m*q
      jt=jf+n*m
      call BDRfsort(mu,n,smod(jf+1),smod(jt+1),sc)
      do 100 inp = 1, np
      ja=q+6
      jb=ja+p*m
      jf=jb+m*q
      jt=jf+n*m      
      do 12181 i=1,q
      y(inp,i)=0.0
12181 continue
      do 12191 l=1,mu
      s=0.0
      do 12201 j=1,p
      s=s+smod(ja+j)*x(inp,j)
12201 continue
      if(s .gt. smod(jt+1)) goto 12221
      place=1
      go to 12230
12221 continue
      if(s .lt. smod(jt+n)) goto 12251
      place=n
      go to 12230
12251 continue
      low=0
      high=n+1
12261 if(low+1.ge.high) goto 12262
      place=(low+high)/2
      t=smod(jt+place)
      if(s.eq.t) goto 12230
      if(s .ge. t) goto 12281
      high=place
      goto 12291
12281 continue
      low=place
12291 continue
12271 continue
      goto 12261
12262 continue
      jfl=jf+low
      jfh=jf+high
      jtl=jt+low
      jth=jt+high
      t=smod(jfl)+(smod(jfh)-smod(jfl))*(s-smod(jtl))  /
     & (smod(jth)-smod(jtl))
      go to 12300
12230 continue
      t=smod(jf+place)
12300 continue
      do 12311 i=1,q
      y(inp,i)=y(inp,i)+smod(jb+i)*t
12311 continue
      ja=ja+p
      jb=jb+q
      jf=jf+n
      jt=jt+n
12191 continue
      do 12321 i=1,q
      y(inp,i)=ys*y(inp,i)+smod(i+5)
12321 continue
100   continue
      return
      end

      subroutine BDRsetsmu
      implicit double precision (a-h, o-z)
      common /BDRspsmooth/ df, gcvpen, ismethod
      ismethod = 0
      end

      subroutine BDRsupsmu (n,x,y,w,iper,span,alpha,smo,sc,edf)
c
c------------------------------------------------------------------
c
c super smoother (Friedman, 1984).
c
c version 10/10/84
c
c coded  and copywrite (c) 1984 by:
c
c                        Jerome H. Friedman
c                     department of statistics
c                               and
c                stanford linear accelerator center
c                        stanford university
c
c all rights reserved.
c
c
c input:
c    n : number of observations (x,y - pairs).
c    x(n) : ordered abscissa values.
c    y(n) : corresponding ordinate (response) values.
c    w(n) : weight for each (x,y) observation.
c    iper : periodic variable flag.
c       iper=1 => x is ordered interval variable.
c       iper=2 => x is a periodic variable with values
c                 in the range (0.0,1.0) and period 1.0.
c    span : smoother span (fraction of observations in window).
c           span=0.0 => automatic (variable) span selection.
c    alpha : controls high frequency (small span) penality
c            used with automatic span selection (bass tone control).
c            (alpha.le.0.0 or alpha.gt.10.0 => no effect.)
c output:
c   smo(n) : smoothed ordinate (response) values.
c scratch:
c   sc(n,7) : internal working storage.
c
c note:
c    for small samples (n < 40) or if there are substantial serial
c    correlations between observations close in x - value, then
c    a prespecified fixed span smoother (span > 0) should be
c    used. reasonable span values are 0.2 to 0.4.
c
c------------------------------------------------------------------
c
      implicit double precision (a-h, o-z)
      dimension x(n),y(n),w(n),smo(n),sc(n,7)
      double precision df, gcvpen
      common /BDRspans/ spans(3) /BDRconsts/ big,sml,eps
      common /BDRspsmooth/ df, gcvpen, ismethod
      if (x(n).gt.x(1)) go to 30
      sy=0.0
      sw=sy
      do 10 j=1,n
      sy=sy+w(j)*y(j)
      sw=sw+w(j)
 10   continue
      a=0.0
      if (sw.gt.0.0) a=sy/sw
      do 20 j=1,n
      smo(j)=a
 20   continue
      return
 30   continue
C     change by BDR
      if (ismethod .ne. 0) then
        call BDRspline(n, x, y, w, smo, edf)
      else
      i=n/4
      j=3*i
      scale=x(j)-x(i)
 40   if (scale.gt.0.0) go to 50
      if (j.lt.n) j=j+1
      if (i.gt.1) i=i-1
      scale=x(j)-x(i)
      go to 40
 50   vsmlsq=(eps*scale)**2
      jper=iper
      if (iper.eq.2.and.(x(1).lt.0.0.or.x(n).gt.1.0)) jper=1
      if (jper.lt.1.or.jper.gt.2) jper=1
      if (span.le.0.0) go to 60
      call BDRsmooth (n,x,y,w,span,jper,vsmlsq,smo,sc)
      return
 60   do 70 i=1,3
      call BDRsmooth (n,x,y,w,spans(i),jper,vsmlsq,sc(1,2*i-1),sc(1,7))
      call BDRsmooth (n,x,sc(1,7),w,spans(2),-jper,vsmlsq,sc(1,2*i),h)
 70   continue
      do 90 j=1,n
      resmin=big
      do 80 i=1,3
      if (sc(j,2*i).ge.resmin) go to 80
      resmin=sc(j,2*i)
      sc(j,7)=spans(i)
 80   continue
      if (alpha.gt.0.0.and.alpha.le.10.0.and.resmin.lt.sc(j,6)
     & .and.resmin.gt.0.0) sc(j,7)=sc(j,7)+(spans(3)-sc(j,7))*
     & max(sml,resmin/sc(j,6))**(10.0-alpha)
 90   continue
      call BDRsmooth (n,x,sc(1,7),w,spans(2),-jper,vsmlsq,sc(1,2),h)
      do 110 j=1,n
      if (sc(j,2).le.spans(1)) sc(j,2)=spans(1)
      if (sc(j,2).ge.spans(3)) sc(j,2)=spans(3)
      f=sc(j,2)-spans(2)
      if (f.ge.0.0) go to 100
      f=-f/(spans(2)-spans(1))
      sc(j,4)=(1.0-f)*sc(j,3)+f*sc(j,1)
      go to 110
 100  f=f/(spans(3)-spans(2))
      sc(j,4)=(1.0-f)*sc(j,3)+f*sc(j,5)
 110  continue
      call BDRsmooth (n,x,sc(1,4),w,spans(1),-jper,vsmlsq,smo,h)
      edf = 0
      endif
      end

      subroutine BDRsmooth (n,x,y,w,span,iper,vsmlsq,smo,acvr)
      implicit double precision (a-h, o-z)
      dimension x(n),y(n),w(n),smo(n),acvr(n)
      integer in,out
      double precision wt,fbo,fbw,xm,ym,tmp,var,cvar,a,h,sy
      xm=0.0
      ym=xm
      var=ym
      cvar=var
      fbw=cvar
      jper=iabs(iper)
      ibw=0.5*span*n+0.5
      if (ibw.lt.2) ibw=2
      it=2*ibw+1
      do 20 i=1,it
      j=i
      if (jper.eq.2) j=i-ibw-1
      xti=x(j)
      if (j.ge.1) go to 10
      j=n+j
      xti=x(j)-1.0
 10   wt=w(j)
      fbo=fbw
      fbw=fbw+wt
      if (fbw.gt.0.0) xm=(fbo*xm+wt*xti)/fbw
      if (fbw.gt.0.0) ym=(fbo*ym+wt*y(j))/fbw
      tmp=0.0
      if (fbo.gt.0.0) tmp=fbw*wt*(xti-xm)/fbo
      var=var+tmp*(xti-xm)
      cvar=cvar+tmp*(y(j)-ym)
 20   continue
      do 80 j=1,n
      out=j-ibw-1
      in=j+ibw
      if ((jper.ne.2).and.(out.lt.1.or.in.gt.n)) go to 60
      if (out.ge.1) go to 30
      out=n+out
      xto=x(out)-1.0
      xti=x(in)
      go to 50
 30   if (in.le.n) go to 40
      in=in-n
      xti=x(in)+1.0
      xto=x(out)
      go to 50
 40   xto=x(out)
      xti=x(in)
 50   wt=w(out)
      fbo=fbw
      fbw=fbw-wt
      tmp=0.0
      if (fbw.gt.0.0) tmp=fbo*wt*(xto-xm)/fbw
      var=var-tmp*(xto-xm)
      cvar=cvar-tmp*(y(out)-ym)
      if (fbw.gt.0.0) xm=(fbo*xm-wt*xto)/fbw
      if (fbw.gt.0.0) ym=(fbo*ym-wt*y(out))/fbw
      wt=w(in)
      fbo=fbw
      fbw=fbw+wt
      if (fbw.gt.0.0) xm=(fbo*xm+wt*xti)/fbw
      if (fbw.gt.0.0) ym=(fbo*ym+wt*y(in))/fbw
      tmp=0.0
      if (fbo.gt.0.0) tmp=fbw*wt*(xti-xm)/fbo
      var=var+tmp*(xti-xm)
      cvar=cvar+tmp*(y(in)-ym)
 60   a=0.0
      if (var.gt.vsmlsq) a=cvar/var
      smo(j)=a*(x(j)-xm)+ym
      if (iper.le.0) go to 80
      h=0.0
      if (fbw.gt.0.0) h=1.0/fbw
      if (var.gt.vsmlsq) h=h+(x(j)-xm)**2/var
      acvr(j)=0.0
      a=1.0-w(j)*h
      if (a.le.0.0) go to 70
      acvr(j)=abs(y(j)-smo(j))/a
      go to 80
 70   if (j.le.1) go to 80
      acvr(j)=acvr(j-1)
 80   continue
      j=1
 90   j0=j
      sy=smo(j)*w(j)
      fbw=w(j)
      if (j.ge.n) go to 110
 100  if (x(j+1).gt.x(j)) go to 110
      j=j+1
      sy=sy+w(j)*smo(j)
      fbw=fbw+w(j)
      if (j.lt.n) go to 100
 110  if (j.le.j0) go to 130
      a=0.0
      if (fbw.gt.0.0) a=sy/fbw
      do 120 i=j0,j
      smo(i)=a
 120  continue
 130  j=j+1
      if (j.le.n) go to 90
      return
      end

      block data BDRbksupsmu
      implicit double precision (a-h, o-z)
      common /BDRspans/ spans(3) /BDRconsts/ big,sml,eps
c
c---------------------------------------------------------------
c
c this sets the compile time (default) values for various
c internal parameters :
c
c spans : span values for the three running linear smoothers.
c spans(1) : tweeter span.
c spans(2) : midrange span.
c spans(3) : woofer span.
c (these span values should be changed only with care.)
c big : a large representable floating point number.
c sml : a small number. should be set so that (sml)**(10.0) does
c       not cause floating point underflow.
c eps : used to numerically stabilize slope calculations for
c       running linear fits.
c
c these parameter values can be changed by declaring the
c relevant labeled common in the main program and resetting
c them with executable statements.
c
c-----------------------------------------------------------------
c
      data spans,big,sml,eps /0.05,0.2,0.5,1.0e20,1.0e-7,1.0e-3/
      end


      subroutine BDRspline (n, x, y, w, smo, edf)
c
c------------------------------------------------------------------
c
c input:
c    n : number of observations.
c    x(n) : ordered abscissa values.
c    y(n) : corresponding ordinate (response) values.
c    w(n) : weight for each (x,y) observation.
c output:
c   smo(n) : smoothed ordinate (response) values.
c   edf : equivalent degrees of freedom
c
c------------------------------------------------------------------
c
      implicit double precision (a-h, o-z)
      dimension x(n),y(n),w(n),smo(n)
      double precision knot(29), coef(25), work((17+25)*25)
      double precision df, df1, lambda, gcvpen, lev(250), s
      double precision dx(250), dy(250), dw(250), dsmo(250)
      integer par(2)
      double precision param(3), crit
      common /BDRspsmooth/ df, gcvpen, ismethod
      do 10 i = 1,n
        dx(i) = (x(i)-x(1))/(x(n)-x(1))
        dy(i) = y(i)
 10     dw(i) = w(i)
      nk = min(n,15)
      knot(1) = dx(1)
      knot(2) = dx(1)
      knot(3) = dx(1)
      knot(4) = dx(1)
      knot(nk+1) = dx(n)
      knot(nk+2) = dx(n)
      knot(nk+3) = dx(n)
      knot(nk+4) = dx(n)
      do 40 i = 5, nk
        p = (n-1)*real(i-4)/real(nk-3)
        ip = int(p)
        p = p-ip
        knot(i) = (1-p)*dx(ip+1) + p*dx(ip+2)
 40     continue
c      call dblepr('knots', 5, knot, nk+4)
C	par = (icrit, ispar)
      if (iabs(ismethod) .eq. 1) then
        par(1) = 3
        df1 = df
      else
         par(1) = 1
         df1 = 0d0
      endif
      par(2) = 0
      param(1) = 0d0
      param(2) = 1.5d0
      param(3) = 1d-2
      isetup = 0
      ier = 1
      call qsbart(gcvpen,df1,dx,dy,dw,n,knot,nk,coef,dsmo,lev,crit,
     & par,lambda,param,isetup, work,4,1,ier)
      if(ier .gt. 0) call intpr('TROUBLE:',8, ier, 1)
      do 50 i = 1,n
 50     smo(i) = dsmo(i)
c      call dblepr('smoothed',8, dsmo, n)
      s = 0
      do 60 i = 1, n
 60      s = s + lev(i)
      edf = s
      if(ismethod.lt.0) then
        call dblepr('lambda', 6, lambda, 1)
        call dblepr('df', 2, s, 1)
      endif
      end


***********************************************************************


      subroutine BDRsort (v,a,ii,jj)
      implicit double precision (a-h, o-z)
c
c     puts into a the permutation vector which sorts v into
c     increasing order.  only elements from ii to jj are considered.
c     arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements
c
c     this is a modification of CACM algorithm #347 by R. C. Singleton,
c     which is a modified Hoare quicksort.
c
      dimension a(jj),v(1),iu(20),il(20)
      integer t,tt
      double precision v
      m=1
      i=ii
      j=jj
 10   if (i.ge.j) go to 80
 20   k=i
      ij=(j+i)/2
      t=a(ij)
      vt=v(ij)
      if (v(i).le.vt) go to 30
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      v(ij)=v(i)
      v(i)=vt
      vt=v(ij)
 30   l=j
      if (v(j).ge.vt) go to 50
      a(ij)=a(j)
      a(j)=t
      t=a(ij)
      v(ij)=v(j)
      v(j)=vt
      vt=v(ij)
      if (v(i).le.vt) go to 50
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      v(ij)=v(i)
      v(i)=vt
      vt=v(ij)
      go to 50
 40   a(l)=a(k)
      a(k)=tt
      v(l)=v(k)
      v(k)=vtt
 50   l=l-1
      if (v(l).gt.vt) go to 50
      tt=a(l)
      vtt=v(l)
 60   k=k+1
      if (v(k).lt.vt) go to 60
      if (k.le.l) go to 40
      if (l-i.le.j-k) go to 70
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 90
 70   il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 90
 80   m=m-1
      if (m.eq.0) return
      i=il(m)
      j=iu(m)
 90   if (j-i.gt.10) go to 20
      if (i.eq.ii) go to 10
      i=i-1
 100  i=i+1
      if (i.eq.j) go to 80
      t=a(i+1)
      vt=v(i+1)
      if (v(i).le.vt) go to 100
      k=i
 110  a(k+1)=a(k)
      v(k+1)=v(k)
      k=k-1
      if (vt.lt.v(k)) go to 110
      a(k+1)=t
      v(k+1)=vt
      go to 100
      end
