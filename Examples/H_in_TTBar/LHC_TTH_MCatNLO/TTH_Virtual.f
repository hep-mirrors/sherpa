
      subroutine m2_virt_tth(pphy,res,mode)
      implicit none
      integer nlegs, j, mode
      parameter (nlegs=5)  
      real*8 pphy(0:3,nlegs)
      real*8 res(4)

      real*8 sig0(2),sig0qq(2),sig0gg(2)
      real*8 sigv(2),sigvqq(2),sigvgg(2)
      real*8 mtop,alphas,mb
      common/topas/mtop,alphas,mb
      real*8 mh
      common/smhiggs/mh
      real*8 delta,muedr
      common/dimreg/delta,muedr
c
c four-momenta and invariants
c
      real*8 pt(4),ptb(4),ph(4),b1(4),b2(4)
      real*8 sp,t1(2),t2(2),u1(2),u2(2)
c
      delta=1d0
c
c four-momenta:
c
      do j=1,3
         b1(j)=pphy(j,1)
         b2(j)=pphy(j,2)
         ph(j)=pphy(j,3)
         pt(j)=pphy(j,4)
         ptb(j)=pphy(j,5)
      enddo
      b1(4)=pphy(0,1)
      b2(4)=pphy(0,2)
      ph(4)=pphy(0,3)
      pt(4)=pphy(0,4)
      ptb(4)=pphy(0,5)
c
c calculate the invariants :
c
      sp=(b1(4)+b2(4))**2-((b1(1)+b2(1))**2+
     $     (b1(2)+b2(2))**2+(b1(3)+b2(3))**2)
      t1(1)=(b1(4)-pt(4))**2-((b1(1)-pt(1))**2+
     $     (b1(2)-pt(2))**2+(b1(3)-pt(3))**2)
      t2(1)=(b2(4)-ptb(4))**2-((b2(1)-ptb(1))**2+
     $     (b2(2)-ptb(2))**2+(b2(3)-ptb(3))**2)
      u1(1)=(b2(4)-pt(4))**2-((b2(1)-pt(1))**2+
     $     (b2(2)-pt(2))**2+(b2(3)-pt(3))**2)
      u2(1)=(b1(4)-ptb(4))**2-((b1(1)-ptb(1))**2+
     $     (b1(2)-ptb(2))**2+(b1(3)-ptb(3))**2)
*
      t1(2)=(b2(4)-pt(4))**2-((b2(1)-pt(1))**2+
     $     (b2(2)-pt(2))**2+(b2(3)-pt(3))**2)
      t2(2)=(b1(4)-ptb(4))**2-((b1(1)-ptb(1))**2+
     $     (b1(2)-ptb(2))**2+(b1(3)-ptb(3))**2)
      u1(2)=(b1(4)-pt(4))**2-((b1(1)-pt(1))**2+
     $     (b1(2)-pt(2))**2+(b1(3)-pt(3))**2)
      u2(2)=(b2(4)-ptb(4))**2-((b2(1)-ptb(1))**2+
     $     (b2(2)-ptb(2))**2+(b2(3)-ptb(3))**2)
c     
c calculate the matrix elements squared at the parton level:	
c qqg=1 -> qqbar , qqg=2 -> gg, qqg=3 -> qqbar and gg
c
      call tthvirt(mode,sp,t1,t2,u1,u2,sig0qq,sig0gg,sigvqq,sigvgg)
      if (mode.eq.1) then
         res(1)=sigvqq(1)
         res(2)=0
         res(3)=0
         res(4)=sig0qq(1)
      elseif (mode.eq.2) then
         res(1)=sigvgg(1)
         res(2)=0
         res(3)=0
         res(4)=sig0gg(1)
      endif

      end

      subroutine tthvirt(qqg,sp,t1,t2,u1,u2,sig0qq,sig0gg,sigvqq,sigvgg)
***********************************************************
c virtual O(alpha_s) corrections to ttbarH production
c October 2011, D.W.
***********************************************************
      implicit none
      integer i,j
      real*8 sig0(2),sig0qq(2),sig0gg(2)
      real*8 sigv(2),sigvqq(2),sigvgg(2)
      real*8 sp,t1(2),t2(2),u1(2),u2(2)
      real*8 mtop,alphas,mb
      common/topas/mtop,alphas,mb
      real*8 mh
      common/smhiggs/mh
      real*8 pi,pi2,w2
      common/para/pi,pi2,w2
      real*8 delta,muedr
      common/dimreg/delta,muedr
      real*8 yukt,yukb,yukf
      common/yukawa/yukt,yukb,yukf
c Coulomb singularity
c      real*8 betamax,dcoul
c      complex*16 spence
c
c o_fac: overall factor (coupling constant)
c
      real*8 o_fac
c matrixelements squared
      real*8 matb_gg,matb_qq
      real*8 matv_gg,matv_qq
c switch for qqbar and/or gluon fusion
      integer qqg
c
c initialization
      do i=1,2
         sig0qq(i)=0d0
         sig0gg(i)=0d0
         sigvqq(i)=0d0
         sigvgg(i)=0d0
      enddo
c constants
      pi=4d0*datan(1d0)
      pi2=pi*pi
      w2=dsqrt(2d0)
c
c overall factor and Yukawa couplings
c
c MSM-Higgs-fermion coupling:
c          eta
      yukf=yukt
c
      o_fac=yukf**2*(4d0*pi*alphas)**2
c     
c Born-matrix element squared:
c     
c qqbar annihilation
      do i=1,2
         sig0(i)=0d0
      enddo
      if (qqg.eq.1.or.qqg.eq.3) then
         sig0(1)=matb_qq(sp,t1(1),t2(1),u1(1),u2(1),mtop,mh)
         sig0(2)=matb_qq(sp,t1(2),t2(2),u1(2),u2(2),mtop,mh)
         do i=1,2
            sig0qq(i)=o_fac*sig0(i)
         enddo
      end if
c gluon fusion
      do i=1,2
         sig0(i)=0d0
      enddo
      if (qqg.eq.2.or.qqg.eq.3) then
         sig0(1)=matb_gg(sp,t1(1),t2(1),u1(1),u2(1),mtop,mh)
         sig0(2)=matb_gg(sp,t1(2),t2(2),u1(2),u2(2),mtop,mh)
         do i=1,2
            sig0gg(i)=o_fac*sig0(i)
         enddo
      endif
c
c Coulomb singularity:
c      betamax=dsqrt((dsqrt(sp)-mh)**2-4d0*mtop**2)/2d0/mtop
c      dcoul=-8d0/6d0*alphas/3d0/betamax
c      do i=1.2
c      sig0qq(i)=sig0qq(i)*dcoul
c      sig0gg(i)=sig0gg(i)*dcoul
c      enddo
c
c virtual O(alpha_s) corrections 
c
c renormalization constants      
c      do j=1,3
c         muedr=10d0**(j**2)
      call renorm(mtop,mh)
c qqbar annihilation:
      do i=1,2
         sigv(i)=0d0
      enddo
      if (qqg.eq.1.or.qqg.eq.3) then
         sigv(1)=matv_qq(sp,t1(1),t2(1),u1(1),u2(1),mtop,mh)
         sigv(2)=matv_qq(sp,t1(2),t2(2),u1(2),u2(2),mtop,mh)
         do i=1,2
            sigvqq(i)=o_fac*sigv(i)*alphas/4d0/pi
         enddo
      endif
c gluon fusion:
      do i=1,2
         sigv(i)=0d0
      enddo
      if (qqg.eq.2.or.qqg.eq.3) then
         sigv(1)=matv_gg(sp,t1(1),t2(1),u1(1),u2(1),mtop,mh)
         sigv(2)=matv_gg(sp,t1(2),t2(2),u1(2),u2(2),mtop,mh)
c        sigv(2)=sigv(1)
c         nlf=5d0
         do i=1,2
            sigvgg(i)=o_fac*sigv(i)*alphas/4d0/pi
c            sigvgg(i)=sigvgg(i)-2d0*sig0gg(i)*alphas/4d0/pi*
c     $        (-2d0/3d0*5d0+11d0/3d0*3d0)*dlog(muedr**2/sp)
         enddo
      endif
c      write(6,*)muedr,sigvqq(1),sigvgg(1)
c      enddo
      return
      end
