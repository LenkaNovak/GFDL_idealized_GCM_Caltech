      subroutine gaqd(nlat,theta,wts,dwork,ldwork,ierror)
c
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1998 by ucar                 .
c  .                                                             .
c  .       university corporation for atmospheric research       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                         spherepack3.0                       .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c
c ... file gaqd.f
c
c     subroutine gaqd computes the nlat gaussian colatitudes and weights
c     in double precision. the colatitudes are in radians and lie in the
c     in the interval (0,pi).
c
c     input parameters
c
c     nlat    the number of gaussian colatitudes in the interval (0,pi)
c             (between the two poles).  nlat must be greater than zero.
c
c     dwork   a double precision temporary work space.
c
c     ldwork  the length of the work space  in the routine calling gaqd
c             ldwork must be at least nlat+1.
c
c     output parameters
c
c     theta   a double precision vector of length nlat containing the
c             nlat gaussian colatitudes on the sphere in increasing radians
c             in the interval (0,pi).
c
c     wts     a double precision vector of length nlat containing the
c             nlat gaussian weights.
c
c     ierror = 0 no errors
c            = 1 if ldwork.lt.nlat+1
c            = 2 if nlat.le.0
c            = 3 if unable to compute gaussian points
c                (failure in in the eigenvalue routine)
c
c  *****************************************************************
c
      dimension dwork(nlat+1),theta(nlat),wts(nlat)
      double precision x,theta,wts,dwork
      ierror = 1
c
c     check work space length
c
      if (ldwork.lt.nlat+1) return
      ierror = 2
      if (nlat.le.0) return
      ierror = 0
c
c     compute weights and points analytically when nlat=1,2
c
      if (nlat.eq.1) then
      theta(1) = dacos(0.0d0)
      wts(1) = 2.d0
      return
      end if
      if (nlat.eq.2) then
      x = dsqrt(1.0d0/3.0d0)
      theta(1) = dacos(x)
      theta(2) = dacos(-x)
      wts(1) = 1.d0
      wts(2) = 1.d0
      return
      end if
c
      call gaqd1(nlat,theta,wts,dwork,ierror)
      if (ierror.ne.0) then
      ierror = 3
      return
      end if
c
c     a newton correction
c
      call gpnewt(nlat,theta,dwork)
c
c     compute weights
c
      call gwts(nlat,theta,wts,dwork)
      return
      end
c
      subroutine gaqd1(n,theta,w,e,ier)
      dimension theta(n),w(n),e(n)
      double precision theta,w,e,temp
c
c     set symmetric tridiagnonal matrix subdiagonal and diagonal
c     coefficients for matrix coming from coefficients in the
c     recursion formula for legendre polynomials
c     a(n)*p(n-1)+b(n)*p(n)+c(n)*p(n+1) = 0.
c
      w(1)=0.
      e(1) = 0.
      do 100 j=2,n
      e(j)= dble((j-1)**2)/dble((j+j-1)*(j+j-3))
      w(j) = 0.
  100 continue
c
c     compute zeros as the eigenvalues of a symmetric
c     tridiagonal matrix
c
      call dval(n,w,e,ier)
      if (ier.ne.0) write(*,11) ier
 11   format(' error in evals =',i5)
c
c     compute gauss legendre points
c
      do 101 j=1,n
      theta(j) = dacos(w(j))
  101 continue
c
c     reverse order of gaussian points to be
c     monotonic increasing in radians
c
      n2 = n/2
      do 102 i=1,n2
      temp = theta(i)
      theta(i) = theta(n-i+1)
      theta(n-i+1) = temp
  102 continue
      return
      end
      subroutine gwts(n,theta,wts,work)
c
c     computes gauss weights as described in swarztrauber
c     and spotz, generalized harmonic transforms
c
      double precision theta(n),wts(n),work(n+1)
c
      call gwts1(n,theta,wts,work,work(n/2+2))
      return
      end
      subroutine gwts1(n,theta,wts,dcp,cp)
      double precision theta(n),wts(n),cp((n-1)/2+1),
     1                 dcp(n/2+1),fn,sqnn,pb,dpb      

      fn = n
      sqnn = dsqrt((fn+fn-1)*(fn+fn+1))
      call lfc (0,n-1,cp)
      call lfc (0,n,dcp)
      do i=1,n
      call lft (0,n-1,theta(i),cp,pb)
      call dlft (0,n,theta(i),dcp,dpb)
      wts(i) = -sqnn*dsin(theta(i))/(fn*pb*dpb)
      end do
      return
      end
      subroutine gpnewt(n,theta,work)
c
c     one newton iteration to correct the gauss points
c
      double precision theta(n),work(n/2+1)
c
      call gpnewt1(n,theta,work)
      return
      end
      subroutine gpnewt1(n,theta,cp)
      double precision theta(n),cp(n/2+1),pb,dpb      
c
      call lfc (0,n,cp)
      do i=1,n
      call lft (0,n,theta(i),cp,pb)
      call dlft (0,n,theta(i),cp,dpb)
      theta(i) = theta(i)-pb/dpb
      end do
      return
      end
      subroutine dval(n,d,e2,ierr)
c
c     name changed from tqlrat to avoid name conflicts
c
      integer i,j,l,m,n,ii,l1,mml,ierr
      double precision d(n),e2(n)
      double precision b,c,f,g,h,p,r,s,t,dzeps,dpytha
c
c     this subroutine is a translation of the algol procedure tqlrat,
c     algorithm 464, comm. acm 16, 689(1973) by reinsch.
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the rational ql method.
c
c     on input
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e2 contains the squares of the subdiagonal elements of the
c          input matrix in its last n-1 positions.  e2(1) is arbitrary.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues.
c
c        e2 has been destroyed.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls dpytha for  sqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1987.
c     modified by c. moler to fix underflow/overflow difficulties,
c     especially on the vax and other machines where epslon(1.0e0)**2
c     nearly underflows.  see the loop involving statement 102 and
c     the two statements just before statement 200.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e2(i-1) = e2(i)
c
      f = 0.0d0
      t = 0.0d0
      e2(n) = 0.0d0
c
      do 290 l = 1, n
         j = 0
         h = dabs(d(l)) + dsqrt(e2(l))
         if (t .gt. h) go to 105
         t = h
         b = dzeps(t)
         c = b * b
         if (c .ne. 0.0d0) go to 105
c        spliting tolerance underflowed.  look for larger value.
         do 102 i = l, n
            h = dabs(d(i)) + dsqrt(e2(i))
            if (h .gt. t) t = h
  102    continue
         b = dzeps(t)
         c = b * b
c     .......... look for small squared sub-diagonal element ..........
  105    do 110 m = l, n
            if (e2(m) .le. c) go to 120
c     .......... e2(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         s = dsqrt(e2(l))
         g = d(l)
         p = (d(l1) - g) / (2.0e0 * s)
         r = dpytha(p,1.0d0)
         d(l) = s / (p + sign(r,p))
         h = g - d(l)
c
         do 140 i = l1, n
  140    d(i) = d(i) - h
c
         f = f + h
c     .......... rational ql transformation ..........
         g = d(m)
         if (g .eq. 0.0d0) g = b
         h = g
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            p = g * h
            r = p + e2(i)
            e2(i+1) = s * r
            s = e2(i) / r
            d(i+1) = h + s * (h + d(i))
            g = d(i) - e2(i) / g
c           avoid division by zero on next pass
            if (g .eq. 0.0e0) g = dzeps(d(i))
            h = g * (p / r)
  200    continue
c
         e2(l) = s * g
         d(l) = h
c     .......... guard against underflow in convergence test ..........
         if (h .eq. 0.0d0) go to 210
         if (dabs(e2(l)) .le. dabs(c/h)) go to 210
         e2(l) = h * e2(l)
         if (e2(l) .ne. 0.0d0) go to 130
  210    p = d(l) + f
c     .......... order eigenvalues ..........
         if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
c
  250    i = 1
  270    d(i) = p
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      double precision function dpytha(a,b)
      double precision a,b
c     dpytha is a double precision modification of pythag off eispack
c     for use by dimtql

c
c     finds sqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dabs(a)
      if (dabs(b).ge.dabs(a)) p = dabs(b)
      if (p .eq. 0.0d0) go to 20
      r = (dabs(a)/p)**2
      if (dabs(b).lt.dabs(a)) r = (dabs(b)/p)**2
   10 continue
	 t = 4.0d0 + r
	 if (t .eq. 4.0d0) go to 20
	 s = r/t
	 u = 1.0d0 + 2.0d0*s
	 p = u*p
	 r = (s/u)**2 * r
      go to 10
   20 dpytha = p
      return
      end
      double precision function dzeps (x)
      double precision  x
c
c     estimate unit roundoff in quantities of size x.
c
      double precision a,b,c,eps
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing floating point
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to 
c            the accuracy used in floating point variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying 
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger floating point number.
c     the developers of eispack would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     this version dated 4/6/83.
c
      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = abs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      dzeps = eps*dabs(x)
      return
      end
      subroutine lfc (m,n,cp)
c
      double precision cp,fnum,fden,fnmh,a1,b1,c1,cp2,fnnp1,fnmsq,fk,
     1       t1,t2,pm1,sc10,sc20,sc40
      dimension       cp(1)
      parameter (sc10=1024.d0)
      parameter (sc20=sc10*sc10)
      parameter (sc40=sc20*sc20)
c
      cp(1) = 0.
      ma = iabs(m)
      if(ma .gt. n) return
      if(n-1) 2,3,5
    2 cp(1) = dsqrt(2.d0)
      return
    3 if(ma .ne. 0) go to 4
      cp(1) = dsqrt(1.5d0)
      return
    4 cp(1) = dsqrt(.75d0)
      if(m .eq. -1) cp(1) = -cp(1)
      return
    5 if(mod(n+ma,2) .ne. 0) go to 10
      nmms2 = (n-ma)/2
      fnum = n+ma+1
      fnmh = n-ma+1
      pm1 = 1.d0
      go to 15
   10 nmms2 = (n-ma-1)/2
      fnum = n+ma+2
      fnmh = n-ma+2
      pm1 = -1.d0
c      t1 = 1.
c      t1 = 2.d0**(n-1)
c      t1 = 1.d0/t1
 15   t1 = 1.d0/sc20
      nex = 20
      fden = 2.d0
      if(nmms2 .lt. 1) go to 20
      do 18 i=1,nmms2
      t1 = fnum*t1/fden
      if(t1 .gt. sc20) then
      t1 = t1/sc40
      nex = nex+40
      end if
      fnum = fnum+2.
      fden = fden+2.
   18 continue
   20 t1 = t1/2.d0**(n-1-nex)
      if(mod(ma/2,2) .ne. 0) t1 = -t1
      t2 = 1. 
      if(ma .eq. 0) go to 26
      do 25 i=1,ma
      t2 = fnmh*t2/(fnmh+pm1)
      fnmh = fnmh+2.
   25 continue
   26 cp2 = t1*dsqrt((n+.5d0)*t2)
      fnnp1 = n*(n+1)
      fnmsq = fnnp1-2.d0*ma*ma
      l = (n+1)/2
      if(mod(n,2) .eq. 0 .and. mod(ma,2) .eq. 0) l = l+1
      cp(l) = cp2
      if(m .ge. 0) go to 29
      if(mod(ma,2) .ne. 0) cp(l) = -cp(l)
   29 if(l .le. 1) return
      fk = n
      a1 = (fk-2.)*(fk-1.)-fnnp1
      b1 = 2.*(fk*fk-fnmsq)
      cp(l-1) = b1*cp(l)/a1
   30 l = l-1
      if(l .le. 1) return
      fk = fk-2.
      a1 = (fk-2.)*(fk-1.)-fnnp1
      b1 = -2.*(fk*fk-fnmsq)
      c1 = (fk+1.)*(fk+2.)-fnnp1
      cp(l-1) = -(b1*cp(l)+c1*cp(l+1))/a1
      go to 30
      end
      subroutine lft (m,n,theta,cp,pb)
      double precision cp(*),pb,theta,cdt,sdt,cth,sth,chh
      cdt = dcos(theta+theta)
      sdt = dsin(theta+theta)
      nmod=mod(n,2)
      mmod=mod(m,2)
      if(nmod)1,1,2
1     if(mmod)3,3,4
c
c     n even, m even
c
3     kdo=n/2
      pb = .5*cp(1)
      if(n .eq. 0) return
      cth = cdt
      sth = sdt
      do 170 k=1,kdo
c     pb = pb+cp(k+1)*dcos(2*k*theta)
      pb = pb+cp(k+1)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  170 continue
      return
c
c     n even, m odd
c
    4 kdo = n/2
      pb = 0.
      cth = cdt
      sth = sdt
      do 180 k=1,kdo
c     pb = pb+cp(k)*dsin(2*k*theta)
      pb = pb+cp(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  180 continue
      return
2     if(mmod)13,13,14
c
c     n odd, m even
c
13    kdo = (n+1)/2
      pb = 0.
      cth = dcos(theta)
      sth = dsin(theta)
      do 190 k=1,kdo
c     pb = pb+cp(k)*dcos((2*k-1)*theta)
      pb = pb+cp(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  190 continue
      return
c
c     n odd, m odd
c
   14 kdo = (n+1)/2
      pb = 0.
      cth = dcos(theta)
      sth = dsin(theta)
      do 200 k=1,kdo
c     pb = pb+cp(k)*dsin((2*k-1)*theta)
      pb = pb+cp(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  200 continue
      return
      end
      subroutine dlft (m,n,theta,cp,pb)
c
c     computes the derivative of pmn(theta) with respect to theta
c
      dimension cp(1)
      double precision cp,pb,theta,cdt,sdt,cth,sth,chh
      cdt = dcos(theta+theta)
      sdt = dsin(theta+theta)
      nmod=mod(n,2)
      mmod=mod(abs(m),2)
      if(nmod)1,1,2
1     if(mmod)3,3,4
c
c     n even, m even
c
3     kdo=n/2
      pb = 0.d0
      if(n .eq. 0) return
      cth = cdt
      sth = sdt
      do 170 k=1,kdo
c     pb = pb+cp(k+1)*dcos(2*k*theta)
      pb = pb-2.d0*k*cp(k+1)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  170 continue
      return
c
c     n even, m odd
c
    4 kdo = n/2
      pb = 0.
      cth = cdt
      sth = sdt
      do 180 k=1,kdo
c     pb = pb+cp(k)*dsin(2*k*theta)
      pb = pb+2.d0*k*cp(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  180 continue
      return
2     if(mmod)13,13,14
c
c     n odd, m even
c
13    kdo = (n+1)/2
      pb = 0.
      cth = dcos(theta)
      sth = dsin(theta)
      do 190 k=1,kdo
c     pb = pb+cp(k)*dcos((2*k-1)*theta)
      pb = pb-(2.d0*k-1)*cp(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  190 continue
      return
c
c     n odd, m odd
c
   14 kdo = (n+1)/2
      pb = 0.
      cth = dcos(theta)
      sth = dsin(theta)
      do 200 k=1,kdo
c     pb = pb+cp(k)*dsin((2*k-1)*theta)
      pb = pb+(2.d0*k-1)*cp(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  200 continue
      return
      end



