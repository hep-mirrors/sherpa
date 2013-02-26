      double precision function df(s,ma,mb)
*
c derivative of the real part of the function f(s,ma,mb)
*
      implicit real*8(a-z)
      parameter(eps=1.d-6)
      if(s.lt.(ma-mb)**2) then
         rplus = dsqrt((ma+mb)**2-s)
         rmin = dsqrt((ma-mb)**2-s)
         df = ((mb**2-ma**2)*dlog(mb**2/ma**2)/2d0/s-
     .   ((rmin**2+2.d0*rmin**2*rplus**2/s+rplus**2)/(2.d0*rmin*rplus))*
     .        dlog((rmin+rplus)**2/(4.d0*ma*mb))-1.d0 )/s
      elseif(s.lt.(ma+mb)**2) then
         rplus = dsqrt((ma+mb)**2/s-1.d0)
         rmin = dsqrt(1.d0-(ma-mb)**2/s)
         df = ((mb**2-ma**2)*dlog(mb**2/ma**2)/2d0/s+
     .        ((rmin**2+2.d0*rmin**2*rplus**2-rplus**2)/(rmin*rplus))*
     .        datan(rmin/rplus)-1.d0 )/s
      else
         rplus = dsqrt(s-(ma+mb)**2)
         rmin = dsqrt(s-(ma-mb)**2)
         df = ((mb**2-ma**2)*dlog(mb**2/ma**2)/2d0/s-
     .   ((rmin**2-2.d0*rmin**2*rplus**2/s+rplus**2)/(2.d0*rmin*rplus))*
     .        dlog((rmin+rplus)**2/(4.d0*ma*mb))-1.d0 )/s
*        y = dsqrt(s-4.d0*ma**2)
*        z = dsqrt(s)
*        df = ((y**2-s)*dlog((z+y)/(z-y))/(2.d0*z*y)-1.d0)/s
      endif
      end
*
      double precision function g(s,ma,mb)
*
c imaginary part of the function f(s,ma,mb)
*
      implicit double precision (a-z)
      pi = 4d0*datan(1d0)
      g = 0.d0
      pm = (dabs(ma)+dabs(mb))**2
      mm = (dabs(ma)-dabs(mb))**2
      if(s.gt.pm) g = pi*dsqrt((s-pm)*(s-mm))/s
      end
*
      double precision function p(s,m)
***************************************************************
* real part of the 1-loop qed vacuumpolarisation contribution *
* from a fermion with mass m                                  *
*                                                             *
*   relation with the function f:                             *
*     p(s,m) = 1/3 - (1 + 2m**2/s)*f(s,m,m)                   *
***************************************************************
      implicit double precision (a-z)
      if(s.eq.0) then
        p = 0.d0
      else if(s.lt.0) then
        x = dsqrt(1.d0-4.d0*m**2/s)
        p = -8.d0/3.d0+x**2-
     .       x*(3.d0-x**2)*dlog(-4.d0*m**2/(s*(1.d0+x)**2))/2.d0
      else if(s.lt.4.d0*m**2) then
        x = dsqrt(4.d0*m**2/s-1.d0)
        p = -8.d0/3.d0-x**2+x*(3.d0+x**2)*datan(1.d0/x)
      else
        x = dsqrt(1.d0-4.d0*m**2/s)
        p = -8.d0/3.d0+x**2-
     .       x*(3.d0-x**2)*dlog(4.d0*m**2/(s*(1.d0+x)**2))/2.d0
      endif
      end
*
      complex*16 function l2(x)
      implicit double precision (a-z)
      parameter(z2=1.64493406684823d0)
* statement function (sp(u)=li2(x), u=-log(1-x) ):
      sp(u) = ((((((-1.d0/10886400.d0)*u**2+1.d0/211680.d0)*u**2-
     .     1.d0/3600.d0)*u**2+1.d0/36.d0)*u-1.d0/4.d0)*u+1.d0)*u
*
      pi = 4d0*datan(1d0)
      rel2 = -7.d0/2.d0-2.d0*x-(x+3.d0/2.d0)*dlog(x**2)
      if(x.gt.1) then
        u  = -dlog(1.d0+1.d0/x)
        iml2 = -pi*(2.d0*(1.d0+x)**2*u+3.d0+2.d0*x)
        rel2 = rel2+2.d0*(1.d0+x)**2*(-sp(u)-u*dlog(x))
      elseif(x.gt.0.d0) then
        u    = -dlog(1.d0+x)
        u1   =  dlog(x)
        iml2 = -pi*(2.d0*(1.d0+x)**2*(u+log(x))+3.d0+2.d0*x)
        rel2 =  rel2+2.d0*(1.d0+x)**2*(sp(u)+z2-u1**2/2.d0-u*u1)
      elseif(x.gt.-1.d0/2.d0) then
        u    = -log(1.d0+x)
        u1   =  log(-x)
        iml2 =  0.d0
        rel2 =  rel2+2.d0*(1.d0+x)**2*(sp(u)-2.d0*z2-u1**2/2.d0-u*u1)
      elseif(x.gt.-2.d0) then
        u    = -log(-x)
        iml2 =  0.d0
        rel2 =  rel2+2.d0*(1.d0+x)**2*(-sp(u)-z2-u**2/2.d0)
      else
        u    = -log(1.d0+1.d0/x)
        iml2 =  0.d0
        rel2 =  rel2+2.d0*(1.d0+x)**2*(-sp(u)-u*log(-x))
      endif
      l2=dcmplx(rel2,iml2)
      end
*
      complex*16 function l3(x)
      implicit double precision (a-z)
      pi = 4d0*datan(1d0)
      if(x.gt.1.d0/4.d0) then
        sq  =sqrt(4.d0*x-1.d0)
        f   =atan(1.d0/sq)
        iml3=0.d0
        rel3=(2.5d0-2.d0*x+(4.d0*x+2.d0)*sq*f-
     $       8.d0*x*(x+2.d0)*f**2)/3.d0
      elseif(x.gt.0.d0) then
        sq  =sqrt(1.d0-4.d0*x)
        f   =log((1.d0+sq)/(1.d0-sq))
        iml3=-pi*((2.d0*x+1.d0)*sq+2.d0*x*(x+2.d0)*f)/3.d0
        rel3=(2.5d0-2.d0*x+(2.d0*x+1.d0)*sq*f+
     .       2.d0*x*(x+2.d0)*(f**2-pi**2))/3.d0
      else
        sq  =sqrt(1.d0-4.d0*x)
        f   =log((sq+1.d0)/(sq-1.d0))
        iml3=0.d0
        rel3=(2.5d0-2.d0*x+(2.d0*x+1.d0)*sq*f+
     $       2.d0*x*(x+2.d0)*f**2)/3.d0
      endif
      l3=dcmplx(rel3,iml3)
      end
*
      double precision function f(s,ma,mb)
*
c real part of the function f(s,ma,mb)
*
      implicit real*8(a-z)
      parameter(eps=1.d-6)
      ma2=ma**2
      mb2=mb**2
      if(abs(s).lt.eps) then
         f=0.0d0
      elseif(dabs(ma).lt.eps) then
         if(s.gt.mb2+eps) then
            f=1.d0+(1.d0-mb2/s)*dlog(1./(s/mb2-1.d0))
         elseif(s.lt.mb2-eps) then
            f=1.d0+(1.d0-mb2/s)*dlog(1./(1.d0-s/mb2))
         else
            f=1.d0
         endif
      elseif(dabs(mb).lt.eps) then
         if(s.gt.ma2+eps) then
            f=1.d0+(1.d0-ma2/s)*dlog(1./(s/ma2-1.d0))
         elseif(s.lt.ma2-eps) then
            f=1.d0+(1.d0-ma2/s)*dlog(1./(1.d0-s/ma2))
         else
            f=1.d0
         endif
      else
         if(dabs(dabs(mb)-dabs(ma)).lt.eps) then
            f=2.d0
         else
            f=1.d0+((ma2-mb2)/s-(ma2+mb2)/(ma2-mb2))*dlog(mb2/ma2)/2d0
         endif
         if(s.ge.(dabs(ma)+dabs(mb))**2) then
            rplus=dsqrt(s-(dabs(ma)+dabs(mb))**2)
            rmin =dsqrt(s-(dabs(ma)-dabs(mb))**2)
            f=f-rplus*rmin*dlog((rplus+rmin)**4/(16.d0*ma2*mb2))/2d0/s
         elseif(s.lt.(dabs(ma)-dabs(mb))**2) then
            rplus=dsqrt((dabs(ma)+dabs(mb))**2-s)
            rmin =dsqrt((dabs(ma)-dabs(mb))**2-s)
            f=f+rplus*rmin*dlog((rplus+rmin)**4/(16.d0*ma2*mb2))/2d0/s
         else
            rplus=dsqrt((dabs(ma)+dabs(mb))**2-s)
            rmin =dsqrt(s-(dabs(ma)-dabs(mb))**2)
            f=f-2.d0*rplus*rmin*datan(rmin/rplus)/s
         endif
      endif
      end
c*********************************************************
*
      subroutine cfunc(s,mf,m1,m2,m3,c0,c1p,c1m,c20,c2p,c2m)
*
c  definition of the invariant functions in the 3-point integrals
c  with equal external masses mf.
c  s = momentum transfer; m1,m2,m3 are the internal masses
c
c                                   p2 ( = pf)
c                        m2   .
c
c               s .     .       m3
c
c                        m1   .
c                                  p1 ( = -pf)
c
c
c  c0 = scalar integral, c1p/m = c1+/-,  c2p/m= c2+/-
c
      implicit real*8(a-z)
      complex*16 c0, c1p,c1m,c20,c2p,c2m,cscal,b0s12,b031,b032,
     &           b132,b131,b1s12,c1,c2
      xmf=mf*mf
      c0=cscal(s,mf,m1,m2,m3)
      call bquer(s,m1,m2,b0s12,b1s12)
      call bquer(xmf,m3,m1,b031,b131)
      call bquer(xmf,m3,m2,b032,b132)
c
c  b0jk := b0(xmf,mj,mk),  b1jk:= b1(xmf,mf,mk)
c  b0s12 := b0(s,m1,m2)
c
c  the c1 functions:
c
      ms=4.d0*xmf-s
      xm1=m1**2
      xm2=m2**2
      xm3=m3**2
      c1=.25d0*dlog(xm3**2/xm1/xm2)+b0s12-.5d0*(b031+b032)
     &   +(xmf+xm3-xm1/2.d0-xm2/2.d0)*c0
      c1p=c1/ms
      c1m=(dlog(xm2/xm1)/2d0+b031-b032+(xm2-xm1)*c0)/2.d0/s
c
c  the c2 functions:
c
      c2=1.d0+b0s12+(xm1+xm2-2.d0*xm3-2.d0*xmf)*c1p+(xm1-xm2)*c1m
     &   +2.d0*xm3*c0
      c20=c2/4.d0
      c2=                     (b131+b132+2.d0*b0s12-.5d0)/4.d0
     &   +(2.d0*xm3-xm1-xm2+2.d0*xmf)/2.d0*c1p-c20
      c2p=c2/ms
      c2=                    -(b131+b132-.5d0)/4.d0
     &   -(xm1-xm2)/2.d0*c1m-c20
      c2m=c2/s
c
      return
      end
c
c********************************************************
c
      subroutine bquer(s,m1,m2,b0,b1)
c
c  b0 and b1 are the (finite) invariant functions in the
c  2-point integrals.
c  s = q**2;  m1,m2 are the internal masses
c
      implicit real*8(a-z)
      external f,g
      complex*16 b0,b1,cf
      xm1=m1**2
      xm2=m2**2
      lm=dlog(xm2/xm1)/2d0
      if ((s.eq.0d0).and.(dabs(dabs(m1)-dabs(m2)).gt.1d-8)
     $     .and.(dabs(m1*m2).gt.0d0)) then
         b0 = dcmplx(1d0+(xm1+xm2)/(xm1-xm2)*lm,0d0)
         b1 = dcmplx(-1d0/4d0-(xm1+xm2)/(xm1-xm2)/4d0-
     $        xm1/(xm1-xm2)*lm-xm1*xm2/(xm1-xm2)**2*lm,0d0)
         goto 20
      endif
      cf=dcmplx(f(s,m1,m2),g(s,m1,m2))
      if (dabs(m1).eq.dabs(m2)) goto 10
      b0=1.d0-(xm2+xm1)/(xm2-xm1)*lm+cf
      b1=-.25d0+xm1/(xm2-xm1)*lm+(xm2-xm1-s)/2.d0/s*cf
      goto 20
10    continue
      b0=cf
c achtung: Aenderung am 14.8.94: factor 1/4 dazugezaehlt!
      b1=-.5d0*cf+1d0/4d0
20    return
      end
c
c********************************************************
c
      subroutine bquer1(x,m1,m2,b0,b1,p0,p1)
c
c  b0 and b1 are the (finite) inariant functions in the
c  2-point integrals, p0 and p1 their derivatives.
c  real parts only, needed for fermion renormalization.
c  x = q**2;  m1,m2 are the internal masses
c
      implicit real*8(a-z)
      external f
      cf=f(x,m1,m2)
      xm1=m1**2
      xm2=m2**2
      lm=dlog(xm2/xm1)/2d0
      if (m1.eq.m2) goto 10
      b0=1.d0-(xm2+xm1)/(xm2-xm1)*lm+cf
      b1=-.25d0+xm1/(xm2-xm1)*lm+(xm2-xm1-x)/2.d0/x*cf
      goto 20
10    b0=cf
      b1=-b0/2.d0+1d0/4d0
20    continue
c
c   calculation of the derivatives:
c
      sm=xm1+xm2
      dm=xm2-xm1
      sm12=(m1+m2)**2
      dm12=(m1-m2)**2
      s=dsqrt(dabs(sm12-x))
      d=dsqrt(dabs(dm12-x))
      klam=(dm*dm/(x*x)-sm/x)/s/d
      anf=-1.d0/x+dm/(x*x)*lm
      if (x.lt.dm12) goto 30
      if (x.gt.sm12) goto 40
      fact=2.d0*datan(d/s)
      goto 41
30    eps=1.d0
      fact=dlog(dabs((s+d)/(s-d)))
      goto 41
40    eps=-1.d0
      fact=-dlog(dabs((s+d)/(s-d)))
41    continue
      deriv=anf-klam*fact
      p0=deriv
      b1p=.5d0-lm-2.d0*b1-b0+(xm2-xm1-x)*deriv
      p1=b1p/2.d0/x
      return
      end
c
c**************************************************************
c                                                             *
c  the scalar vertex integral with equal external masses mf   *
c                                                             *
c**************************************************************
c
      complex*16 function cscal(s,mf,m1,m2,m3)
c
c  s = momentum transfer; m1,m2,m3  are the internal masses
c
      implicit real*8 (a-y)
      complex*16 z1,z2,z11,z12,z21,z22,cl1,cl2,cl3,cspen,spence,int
      xmf=mf*mf
      if (mf.lt.1d-1) then
         mfstrich = 1d-1
         xmf=mfstrich**2
      end if
c.........xmf etc.   are fermion and boson masses squared
      xm1=m1*m1
      xm2=m2*m2
      xm3=m3*m3
c
c..t'hooft-veltman parameters
      a=1.d0
      b=xmf/s
      c=-1.d0
      d=xm1-xm2-s
      e=xm3-xm1-xmf+s
      f=xm2/s
      d=d/s
      e=e/s
c..discriminante for alpha-equation
      disc=c*c-4.d0*a*b
      if (disc .lt. 0.d0) goto 500
      al=(-c-dsqrt(disc))/2.d0/b
      nenner=c+2.d0*al*b
c..the first integral.............................................
      y0=-(d+e*al+2.d0*a+c*al)/nenner
      y01=y0-1.d0
      d1=(c+e)**2-4.d0*b*(a+d+f)
      x1=-(c+e)/2.d0/b
      if (d1.gt.0.d0) goto 10
c.......complex zeroes of logarithms
      sq1=dsqrt(-d1)
      x2=sq1/2.d0/b
      z1=dcmplx(x1,x2)
      z2=dcmplx(x1,-x2)
      z11=y0/(y0-z1)
      z12=y01/(y0-z1)
      z21=y0/(y0-z2)
      z22=y01/(y0-z2)
      cl1=spence(z11)-spence(z12)+spence(z21)-spence(z22)
      goto 15
10    continue
c........real zeroes
      sq1=dsqrt(d1)
      x2=sq1/2.d0/b
      y1=x1+x2
      y2=x1-x2
      sig1= y0/dabs(y0)
      sig2= y01/dabs(y01)
      y11=y0/(y0-y1)
      y12=y01/(y0-y1)
      y21=y0/(y0-y2)
      y22=y01/(y0-y2)
      cl1=cspen(y11,sig1)-cspen(y12,sig2)+cspen(y21,-sig1)
     &   -cspen(y22,-sig2)
15    continue
c..the second integral............................................
      y0=-(d+e*al)/nenner/(1.d0-al)
      y01=y0-1.d0
      d2=(e+d)**2-4.d0*f*(a+b+c)
      x1=-(e+d)/2.d0/(a+b+c)
      if(d2.gt.0.d0) goto 20
c.......complex zeroes of logarithms
      sq2=dsqrt(-d2)
      x2=sq2/2.d0/(a+b+c)
      z1=dcmplx(x1,x2)
      z2=dcmplx(x1,-x2)
      z11=y0/(y0-z1)
      z12=y01/(y0-z1)
      z21=y0/(y0-z2)
      z22=y01/(y0-z2)
      cl2=spence(z11)-spence(z12)+spence(z21)-spence(z22)
      goto 25
20    continue
c........real zeroes
      x2=dsqrt(d2)/2.d0/(a+b+c)
      y1=x1+x2
      y2=x1-x2
      y11=y0/(y0-y1)
      y12=y01/(y0-y1)
      y21=y0/(y0-y2)
      y22=y01/(y0-y2)
      sig1= y0/dabs(y0)
      sig2= y01/dabs(y01)
      cl2=cspen(y11,sig1)-cspen(y12,sig2)+cspen(y21,-sig1)
     &   -cspen(y22,-sig2)
25    continue
c..the third integral............................................
      y0=(d+e*al)/nenner/al
      y01=y0-1.d0
      d3=d*d-4.d0*a*f
      x1=-d/2.d0/a
      if (d3.gt.0.d0) goto 30
c........complex zeroes of logarithms
      sq3=dsqrt(-d3)
      x2=sq3/2.d0/a
      z1=dcmplx(x1,x2)
      z2=dcmplx(x1,-x2)
      z11=y0/(y0-z1)
      z12=y01/(y0-z1)
      z21=y0/(y0-z2)
      z22=y01/(y0-z2)
      cl3=spence(z11)-spence(z12)+spence(z21)-spence(z22)
      goto 35
30    continue
c........real zeroes
      x2=dsqrt(d3)/2.d0/a
      y1=x1+x2
      y2=x1-x2
 31   format(1h ,3e12.4)
      y11=y0 /(y0-y1)
      y12=y01/(y0-y1)
      y21=y0/(y0-y2)
      y22=y01/(y0-y2)
      sig1= y0/dabs(y0)
      sig2= y01/dabs(y01)
      cl3=cspen(y11,sig1)-cspen(y12,sig2)+cspen(y21,-sig1)
     &   -cspen(y22,-sig2)
35    continue
c..summation of the 3 integrals ....................................
      int=-cl1+cl2-cl3
      cscal=int/nenner/s
      goto 501
500   continue
c..error message for complex alpha................................
      write(6,21)
21    format(1h ,'  i cannot handle a complex alpha')
      stop
501   return
      end
c
c*******************************************************************
c
      complex*16 function cspen(x,sig)
      implicit real*8(a-y)
      complex*16 z,cpi,spence,zx
      pi=4d0*datan(1d0)
      pi6=pi*pi/6.d0
      cpi=dcmplx(0.d0,pi)
      if (x.lt.1.d0) goto 10
      if(x.eq.1.d0) goto 11
      lx=dlog(x)
      x1=1.d0-x
      lx1=dlog(-x1)
      z=dcmplx(x1,0.d0)
      if (sig.gt.0.d0) goto 5
      cspen=-spence(z)+pi6-lx*(lx1+cpi)
      goto 20
5     cspen=-spence(z)+pi6-lx*(lx1-cpi)
      goto 20
10    zx=dcmplx(x,0.d0)
      cspen=spence(zx)
      goto 20
11    cspen=dcmplx(pi6,0.d0)
20    return
      end
c
c**************************************************************
c
      complex*16 function cln(x,sig)
      implicit real*8(a-z)
      complex*16 cpi
      pi=4d0*datan(1d0)
      cpi=dcmplx(0.d0,pi)
      if (x.gt.0.d0) goto 10
      x1=-x
      if (sig.gt.0.d0) goto 5
      cln=dlog(x1)-cpi
      goto 20
5     cln=dlog(x1)+cpi
      goto 20
10    y=dlog(x)
      cln=dcmplx(y,0.d0)
20    return
      end
      complex*16 function spence(xx)
c  the dilogarithm for general complex argument.
c  not allowed: real(xx) gt 1 with imag(xx)=0.
      implicit real*8(a-z)
      integer n
      complex*16 xx,x,z,d,p,cdlog
      dimension a(19)
      pi=4d0*datan(1d0)
      x=xx
      xr=dreal(x)
      xi=dimag(x)
      if(xr.ne.1.) goto 111
      if(xi.eq.0.) goto 20
111   continue
c    projection into the convergence radius
      vor=1.d0
      p=dcmplx(0.d0,0.d0)
      r=dreal(x)
      if (r .le. 0.5d0) goto 1
      p=pi*pi/6.d0- cdlog(x)*cdlog(1.d0-x)
      vor=-1.d0
      x=1.d0-x
    1 continue
      b=cdabs(x)
      if (b .lt. 1.d0) goto 2
      p=p - (pi*pi/6.d0+ cdlog(-x)*cdlog(-x)/2.d0)*vor
      vor=vor*(-1.d0)
      x=1.d0/x
    2 continue
c    calculation of the spence function
      a(1)=1.d0
      a(2)=-0.5d0
      a(3)=1.d0/6.d0
      a(5)=-1.d0/30.d0
      a(7)=1.d0/42.d0
      a(9)=-1.d0/30.d0
      a(11)=5.d0/66.d0
      a(13)=-691.d0/2730.d0
      a(15)=7.d0/6.d0
      a(17)=-3617.d0/510.d0
      a(19)=43867.d0/798.d0
      do 5 n=2,9,1
      a(2*n)=0.d0
    5 continue
      z=(-1.d0)*cdlog(1.d0-x)
      d=dcmplx(a(19),0.d0)
      do 10 n=1,18,1
      d=d*z/(20.d0-n) + a(19-n)
   10 continue
      d=d*z
      spence=d*vor + p
      goto 30
   20 continue
      spence=pi*pi/6.d0
   30 continue
      return
      end
