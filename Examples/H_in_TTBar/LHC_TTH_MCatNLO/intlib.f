c---------------------------------------------------------------------
c     checked version of "integrallib.f"
c---------------------------------------------------------------------
c====================================================================
c     checked with hollik's routines ; status 16.2.93 ; 
c====================================================================
c--------------------------------------------------------------------
      subroutine fint(s,m1,m2,f)
c--------------------------------------------------------------------
c     the f-function as specified in : boehm et.al., 
c     fortschr.phys.34 ; and in : hollik, fortschr.phys. 38 
c--------------------------------------------------------------------
      implicit real*8(a-z)
      parameter (eps=1d-6)
      complex*16 f    
      pi = 4d0*datan(1d0)
      x1 = m1**2
      x2 = m2**2
      p = s-(dabs(m1)+dabs(m2))**2
      m = s-(dabs(m1)-dabs(m2))**2
      imagf = 0d0
      rest = 0d0
      if (abs(s).lt.eps) then
         realf = 0d0
         goto 999
      endif
c
      if (x1.lt.eps) then
         if (x2.lt.eps) then
            realf = 0d0
         else
            if (s.gt.x2+eps) then
               realf = 1d0+(1d0-x2/s)*dlog(1d0/(s/x2-1d0))
               imagf = pi*(1d0-x2/s)
            else if (s.lt.x2-eps) then
               realf = 1d0+(1d0-x2/s)*dlog(1d0/(1d0-s/x2))
            else
               realf = 1d0
            endif 
         endif
         goto 999
      else if (x2.lt.eps) then
         if (x1.lt.eps) then
            realf = 0d0
         else
            if (s.gt.x1+eps) then
               realf = 1d0+(1d0-x1/s)*dlog(1d0/(s/x1-1d0))
               imagf = pi*(1d0-x1/s)
            else if (s.lt.x1-eps) then
               realf = 1d0+(1d0-x1/s)*dlog(1d0/(1d0-s/x1))
            else
               realf = 1d0
            endif 
         endif
         goto 999
      endif
c
      if ((abs(s).lt.x1/1d3).or.(abs(s).lt.x2/1d3)) then
         if (dabs(dabs(m1)-dabs(m2)).lt.eps) then
            realf = s/6d0/x1
         else
            realf = s/(x1-x2)**2*((x1+x2)/2d0
     $           -x1*x2/(x1-x2)*dlog(x1/x2))
         endif
         goto 999
      endif
c
      if (dabs(dabs(m1)-dabs(m2)).lt.eps) then
         rest = 2d0
      else
         rest = 1d0-((x1-x2)/s-(x1+x2)/(x1-x2))*dlog(x1/x2)/2d0
      endif
c
      if (m.lt.0d0) then
         realf = dsqrt(p*m)*dlog((dsqrt(-p)+dsqrt(-m))**2
     $        /(4d0*dsqrt(x1*x2)))/s
      else if (p.lt.0d0) then
         realf = -2d0*dsqrt(-p*m)*datan(dsqrt(-m/p))/s
      else
         realf = -dsqrt(p*m)*dlog((dsqrt(p)+dsqrt(m))**2
     $        /(4d0*dsqrt(x1*x2)))/s
         imagf = dsqrt(p*m)/s*pi
      endif
 999  continue
      f = dcmplx(rest+realf,imagf)
      end
c---------------------------------------------------------------------
c---------------------------------------------------------------------
      subroutine afunc(m,a0)
      implicit real*8(a-z)
      parameter (eps=1d-4)
      common / dimreg / delta,muedr
      if (m**2.lt.eps) then 
         a0 = 0d0
      else
         a0 = m**2*(1d0-dlog(m**2/muedr**2))
      endif
      end
c---------------------------------------------------------------------
      subroutine bfunc(s,m1,m2,b0,b1,b20)
c---------------------------------------------------------------------
c     b0 and b1 are the b0"quer" and b1"quer"-quantities defined 
c     in : hollik, fortschr.phys. 38 (but with log-terms incl.)
c     b20 is the g_{mu}.{nu} coefficient of the 
c     tensor-2-point-integral (log-terms incl, too)
c---------------------------------------------------------------------
      implicit real*8(a-z)
      parameter (eps=1d-4)
      complex*16 b0,b1,b20,ff
      common / dimreg / delta,muedr
      pi = 4d0*datan(1d0)
      x1 = m1**2
      x2 = m2**2
      if (dabs(s).lt.eps) then
         call bnull(m1,m2,realb0)
         call beins(m1,m2,realb1)
         call bzwnull(m1,m2,realb20)
         b0 = dcmplx(realb0,0d0)
         b1 = dcmplx(realb1,0d0)
         b20 = dcmplx(realb20,0d0)
         goto 999
      endif
      if ((x1.lt.eps).and.(x2.lt.eps)) goto 47 
      call fint(s,m1,m2,ff)
 47   continue
      if (x1.lt.eps) then
         if (x2.lt.eps) then
            b0 = 2d0-dlog(s/muedr**2)+dcmplx(0d0,1d0)*pi
            b1 = -b0/2d0
         else
            b0 = 1d0+ff-dlog(x2/muedr**2)
            b1 = -(1d0+(1d0-x2/s)*ff-dlog(x2/muedr**2))/2d0
         endif
         goto 666
      endif
      if (x2.lt.eps) then
         if (x1.lt.eps) then
            b0 = 2d0-dlog(s/muedr**2)+dcmplx(0d0,1d0)*pi
            b1 = -b0/2d0
         else
            b0 = 1d0+ff-dlog(x1/muedr**2)
            b1 = -(1d0+(1d0+x1/s)*ff-dlog(x1/muedr**2))/2d0
         endif
         goto 666
      endif
      if (dabs(x1-x2).lt.1d3*eps) then
         b0 = ff
         b1 = -ff/2d0
      else
         b0 = 1d0-(x1+x2)/(x1-x2)*dlog(x1/x2)/2d0+ff
         b1 = -1d0/2d0+x1/(x1-x2)*dlog(x1/x2)/2d0-(s+x1-x2)/2d0/s*ff
      endif
      b0 = b0 - dlog(x1*x2/muedr**4)/2d0
      b1 = b1 + dlog(x2/muedr**2)/2d0 
c     
 666  continue
      call afunc(m2,a0)
      b20 = a0/6d0+x1/3d0*b0
     $     +(s+x1-x2)/6d0*b1+(x1+x2-s/3d0)/6d0
 999  continue
      end
c---------------------------------------------------------------------
      subroutine bderiv(x,m1,m2,db0,db1,db20)
c---------------------------------------------------------------------
c     real parts of d(b0(s,m1,m2))/ds , d(b1(s,m1,m2))/ds
c     and d(b20(s,m1,m2))/ds
c     (s -> x in the code)
c---------------------------------------------------------------------
      implicit real*8(a-z)
      complex*16 cf,b0,b1,b20
      parameter (eps=1d-1)
*
      pi = 4d0*datan(1d0)
      xm1 = m1**2
      xm2 = m2**2
*
      if (x.lt.eps) then
         if (dabs(xm1-xm2).lt.eps) then
            db0 = 1d0/6d0/xm2
            db1 = -1d0/12d0/xm2
         else if (xm1.lt.eps) then
            db0 = 1d0/2d0/xm2
            db1 = -1d0/6d0/xm2
         else if (xm2.lt.eps) then
            db0 = 1d0/2d0/xm1
            db1 = -1d0/6d0/xm1
         else
            db0 = (xm1**2-xm2**2-2d0*xm1*xm2*dlog(xm1/xm2))
     $           /2d0/(xm1-xm2)**3
            fss0 = (xm1**3-xm2**3+9d0*xm1*xm2*(xm1-xm2)
     $           -6d0*xm1*xm2*dlog(xm1/xm2)*(xm1+xm2))
     $           /3d0/(xm1-xm2)**5
            db1 = -db0/2d0-(xm1-xm2)/4d0*fss0
         endif
         goto 11
      endif
 13   continue
      if ((x.lt.xm1/1d3).or.(x.lt.xm2/1d3)) then
         if (dabs(xm1-xm2).lt.eps) then
            db0 = 1d0/6d0/xm1
         else
            db0 = ((xm1+xm2)/2d0-xm1*xm2/(xm1-xm2)*dlog(xm1/xm2))
     $           /(xm1-xm2)**2
         endif
         db1 = -db0/2d0
         goto 11
      endif
c
      if ((xm1.lt.eps).and.(xm2.lt.eps)) then
         db0 = -1d0/x
         db1 = -db0/2d0
         goto 11
      endif
c
      if (xm2.lt.eps) then
         if (x.gt.xm1+eps) then
            deriv = -(1d0+xm1/x*dlog(x/xm1-1d0))/x
         else if (x.lt.xm1-eps) then
            deriv = -(1d0+xm1/x*dlog(1d0-x/xm1))/x
         else
            deriv = 0d0
         endif 
         goto 10
      endif
c
      if (xm1.lt.eps) then
         if (x.gt.xm2+eps) then
            deriv = -(1d0+xm2/x*dlog(x/xm2-1d0))/x
         else if (x.lt.xm2-eps) then
            deriv = -(1d0+xm2/x*dlog(1d0-x/xm2))/x
         else
            deriv = 0d0
         endif 
         goto 10
      endif
c
      sm = xm1+xm2
      dm = xm2-xm1
      sm12 = (dabs(m1)+dabs(m2))**2
      dm12 = (dabs(m1)-dabs(m2))**2
      lm = dlog(xm2/xm1)/2d0
*
      if (x.lt.dm12) then
         s = dsqrt(sm12-x)
         d = dsqrt(dm12-x)
         fact = dlog((d+s)**2/(4d0*dabs(m1)*dabs(m2)))
         deriv = (dm*lm/x-((d**2+s**2+2d0*d**2*s**2/x)
     $        /(2d0*s*d))*fact-1d0)/x
      else if (x.lt.sm12) then
         if (dabs(xm1-xm2).lt.eps) then
            if (x.lt.eps) then
               deriv = 1d0/6d0/xm1
            else
               sx = 4d0*xm1/x
               deriv = (sx/dsqrt(sx-1d0)
     $              *datan(1d0/dsqrt(sx-1d0))-1d0)/x
            endif
         else
*
c Achtung: fuer x=sm12 geht deriv -> infty! 
*
            s = dsqrt(sm12/x-1d0)
            d = dsqrt(1d0-dm12/x)
            fact = datan(d/s)
            deriv = (dm*lm/x+((d**2-s**2+2d0*d**2*s**2)
     $           /(s*d))*fact-1d0)/x
         endif
      else
         s = dsqrt(x-sm12)
         d = dsqrt(x-dm12)
         fact = dlog((d+s)**2/(4d0*dabs(m1)*dabs(m2)))
         deriv = (dm*lm/x-((d**2+s**2-2d0*d**2*s**2/x)
     $        /(2d0*s*d))*fact-1d0)/x
      endif
*
 10   continue
c--------------------------------------------------------------------
      db0 = deriv
      call fint(x,m1,m2,cf)
      f = dreal(cf)
      db1 = -(x+xm1-xm2)/2d0/x*deriv+(xm1-xm2)/2d0/x**2*f
c     
 11   continue
c
      call bfunc(x,m1,m2,b0,b1,b20)
      db20 = xm1/3d0*db0+dreal(b1)/6d0+(x+xm1-xm2)/6d0*db1
     $     -1d0/18d0
c
 99   continue
      end
c---------------------------------------------------------------------
      subroutine cfuncnew(m,s,m1,m2,m3,c0,c1p,c1m,c20,c2p,c2m)
c---------------------------------------------------------------------
c     contains the scalar coefficient functions of the 
c     tensor three point integral
c---------------------------------------------------------------------
      implicit real*8(a-z)
      complex*16 c0,c1p,c1m,c20,c2p,c2m
     $     ,b0x32,b0x31,b0s12,b1x32,b1x31,adummy,bdummy
      x = m**2
      x1 = m1**2
      x2 = m2**2
      x3 = m3**2
c
      call bfunc(s,m1,m2,b0s12,adummy,bdummy)
      call bfunc(x,m3,m2,b0x32,b1x32,adummy)
      call bfunc(x,m3,m1,b0x31,b1x31,adummy)
      call c0int(m,s,m1,m2,m3,c0)
c
      c1m = ((b0x32+b0x31)/2d0-b0s12
     $     -(x+x3-(x1+x2)/2d0)*c0)/(4d0*x-s)
      c1p = ((b0x32-b0x31)/2d0+(x1-x2)/2d0*c0)/s
      c20 = b0s12/4d0+x3/2d0*c0-(x1-x2)/4d0*c1p
     $     -(x1+x2-2d0*(x+x3))/4d0*c1m+1d0/4d0
      c2m = (-c20+(b1x32+b1x31)/4d0+b0s12/2d0
     $     -(x+x3-(x1+x2)/2d0)*c1m)/(4d0*x-s)
      c2p = (-c20-(b1x32+b1x31)/4d0+(x1-x2)/2d0*c1p)/s
      end
c---------------------------------------------------------------------
      subroutine c0int(mf,s,m1,m2,m3,c0)
c**************************************************************
c                                                             *
c  the scalar vertex integral with equal external masses mf   *
c                                                             *
c-------------------------------------------------------------------
c     s = momentum transfer; m1,m2,m3  are the internal masses
c
      implicit real*8 (a-y)
      complex*16 z1,z2,z11,z12,z21,z22,cl1,cl2,cl3,cspen,spence,
     &     int,c0
      xmf=mf*mf
      if (mf.lt.1d-1) then
         mfstrich = 1d-1
         xmf = mfstrich**2
      endif
c     xm's : are fermion and boson masses squared
      xm1=m1*m1
      xm2=m2*m2
      xm3=m3*m3
c     the t'hooft-veltman parameters
      a=1.d0
      b=xmf/s
      c=-1.d0
      d=xm1-xm2-s
      e=xm3-xm1-xmf+s
      f=xm2/s
      d=d/s
      e=e/s
c     discriminante for alpha-equation
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
      int = -cl1+cl2-cl3
      c0 = int/nenner/s
      goto 501
500   continue
c..error message for complex alpha................................
      write(6,21)
21    format(1h ,'  i cannot handle a complex alpha')
      stop
501   return
      end
c--------------------------------------------------------------------
c     the x"null" subroutines calculate the b0-, b1-, c0-,
c     d0- and d20-integrals at zero external momenta  
c--------------------------------------------------------------------
      subroutine bnull(a,b,b0)
      implicit real*8(a-z)
      parameter (eps=1d-8)      
      common / dimreg / delta,muedr
      xa = a**2
      xb = b**2
      x = xa+xb
      if (x.lt.eps) then
c        write(6,*)'all args zero is not allowed for b0-function !'
         b0 = 0d0
         goto 2
      endif
      if (xa*xb.eq.0d0) then
         zw = 1d0-dlog(x/muedr**2)
      else
         zw = -dlog(xa*xb/muedr**4)/2d0
         if (dabs(dabs(a)-dabs(b)).gt.eps) then
            zw = zw+1d0-(xa+xb)/(xa-xb)*dlog(xa/xb)/2d0
         endif
      endif
      b0 = zw
 2    continue
      end
c-----------------------------------------------------------------     
      subroutine beins(a,b,b1)
      implicit real*8(a-z)
      parameter (eps=1d-8)
      common / dimreg / delta,muedr
      xa = a**2
      xb = b**2
      x = xa+xb
      if (x.lt.eps) then
c        write(6,*)'all args zero is not allowed for b1-function !'
         b1 = 0d0
         goto 2
      endif      
      if (xa.eq.0d0) then
         zw = -(1d0/2d0-dlog(xb/muedr**2))/2d0
      else if (xb.eq.0d0) then
         zw = -(3d0/2d0-dlog(xa/muedr**2))/2d0
      else
         zw = dlog(xb/muedr**2)/2d0
         if (dabs(dabs(a)-dabs(b)).gt.eps) then
            zw = zw -(xa+xb)/(xa-xb)/4d0-1d0/2d0
     $           +xa/(xa-xb)*dlog(xa/xb)/2d0
     $           +xa*xb/(xa-xb)**2*dlog(xa/xb)/2d0
         endif
      endif
      b1 = zw
 2    continue
      end
c-----------------------------------------------------------------     
      subroutine bzwnull(a,b,b20)
      implicit real*8(a-z)
      common / dimreg / delta,muedr
      xa = a**2
      xb = b**2
      call bnull(a,b,b0)
      call beins(a,b,b1)
      zw = xa*b0/3d0+(xa-xb)*b1/6d0+(xa+xb)/6d0
      if (xb.ne.0d0) then 
         zw = zw+xb*(1d0-dlog(xb/muedr**2))/6d0
      endif
      b20 = zw
      end
c-----------------------------------------------------------------
      subroutine cnull(a,b,c,c0)
      implicit real*8(a-z)
      parameter (eps=1d-8)
      xa = a**2
      xb = b**2
      xc = c**2
      if (dabs(dabs(a)-dabs(b)).lt.eps) then
         if (dabs(dabs(a)-dabs(c)).lt.eps) then
            zw = -1d0/2d0/xa
         else 
            zw = (-1d0+xc/(xa-xc)*dlog(xa/xc))/(xa-xc)
         endif
      else 
         call bnull(a,c,bac)
         call bnull(b,c,bbc)
         zw = (bac-bbc)/(xa-xb)
      endif
      c0 = zw
      end
c-----------------------------------------------------------------     
      subroutine dnull(a,b,c,d,d0)
      implicit real*8(a-z)
      integer i,j
      real*8 m(4)
      parameter (eps=1d-8)
      m(1) = a
      m(2) = b
      m(3) = c
      m(4) = d
      do 30 i = 1,4
         do 40 j = i+1,4
            if (dabs(dabs(m(i))-dabs(m(j))).lt.eps) then
               m(i) = m(i)+eps
               m(j) = m(j)-eps
            endif
 40      continue
 30   continue
      call cnull(a,b,c,c0abc)
      call cnull(a,b,d,c0abd)
      zw = (c0abc-c0abd)/(c**2-d**2)
      d0 = zw
      end
c-----------------------------------------------------------------
      subroutine dzwnull(a,b,c,d,d20)
      implicit real*8(a-z)
      call cnull(b,c,d,c0bcd)
      call dnull(a,b,c,d,d0abcd)
      zw = c0bcd+a**2*d0abcd
      d20 = zw
      end
c-----------------------------------------------------------------     
