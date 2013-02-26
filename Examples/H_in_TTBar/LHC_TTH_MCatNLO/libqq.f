      double precision function matb_qq(s,t1,t2,u1,u2,mt,mh)
      implicit none
      real*8 mt,mh,s,t1,t2,u1,u2
      real*8 pt(4),ptb(4),ph(4),p3(4),p4(4)
      real*8 cs1,cs2,s3,s1,s2,p1p2,p1p3,p2p3,p1p4,p2p4

      matb_qq=0d0

      s3=4d0*mt**2+mh**2-s-t1-t2-u1-u2
      p1p2=(s3-2d0*mt**2)/2d0
      p1p3=(mt**2-t2)/2d0
      p1p4=(mt**2-u2)/2d0
      p2p3=(mt**2-u1)/2d0
      p2p4=(mt**2-t1)/2d0
      s1=s+t2+u2-mt**2
      s2=s+t1+u1-mt**2
      cs1=1d0/(s1-mt**2)
      cs2=1d0/(s2-mt**2)

      matb_qq=cs1**2*(32*mt**4*s + 16*mt**2*s**2 - 32*mt**2*s*p1p2 - 
     $     32*s*p1p2*p1p3 - 32*s*p1p2*p1p4 - 
     $     128*mt**2*p1p3*p1p4 + 128*p1p2*p1p3*p1p4 + 
     $     32*mt**2*s*p2p3 + 32*s*p1p3*p2p3 - 
     $     64*p1p3*p1p4*p2p3 + 64*p1p4**2*p2p3 + 
     $     32*mt**2*s*p2p4 + 64*p1p3**2*p2p4 + 
     $     32*s*p1p4*p2p4 - 64*p1p3*p1p4*p2p4) + 
     $     cs1*cs2*(32*mt**2*s**2 - 64*mt**2*s*p1p2 + 64*s*p1p2**2 + 
     $     32*mt**2*s*p1p3 - 32*s*p1p2*p1p3 + 
     $     32*mt**2*s*p1p4 - 32*s*p1p2*p1p4 + 
     $     32*mt**2*s*p2p3 - 32*s*p1p2*p2p3 + 
     $     64*s*p1p3*p2p3 + 128*mt**2*p1p4*p2p3 - 
     $     128*p1p2*p1p4*p2p3 - 64*p1p3*p1p4*p2p3 + 
     $     64*p1p4**2*p2p3 + 64*p1p4*p2p3**2 + 
     $     32*mt**2*s*p2p4 - 32*s*p1p2*p2p4 + 
     $     128*mt**2*p1p3*p2p4 - 128*p1p2*p1p3*p2p4 + 
     $     64*p1p3**2*p2p4 + 64*s*p1p4*p2p4 - 
     $     64*p1p3*p1p4*p2p4 - 64*p1p3*p2p3*p2p4 - 
     $     64*p1p4*p2p3*p2p4 + 64*p1p3*p2p4**2) + 
     $     cs2**2*(32*mt**4*s + 16*mt**2*s**2 - 32*mt**2*s*p1p2 + 
     $     32*mt**2*s*p1p3 + 32*mt**2*s*p1p4 - 
     $     32*s*p1p2*p2p3 + 32*s*p1p3*p2p3 + 
     $     64*p1p4*p2p3**2 - 32*s*p1p2*p2p4 + 
     $     32*s*p1p4*p2p4 - 128*mt**2*p2p3*p2p4 + 
     $     128*p1p2*p2p3*p2p4 - 64*p1p3*p2p3*p2p4 - 
     $     64*p1p4*p2p3*p2p4 + 64*p1p3*p2p4**2)
      
      matb_qq=matb_qq*2d0/9d0/s**2/4d0
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      double precision function matv_qq(s,t1,t2,u1,u2,mt,mh)
      implicit none
      integer i,j
      real*8 s,s1,s2,t1,t2,u1,u2,mt,mh
      real*8 cs1,cs2
      real*8 ncf,cg,cf
      real*8 mvs,mat(15)
      real*8 colfacb,colfac_b3(4),colfac_p(2)
      complex*16 kmvs1(15),kmvs2(15),kmvs3(15),kmvb3(4),kmvp(2),zero
      real*8 pt(4),ptb(4),ph(4),p3(4),p4(4)

      matv_qq=0d0
      zero=dcmplx(0d0,0d0)
c initialization of coefficients:
      do i=1,15
         kmvs1(i)=zero
         kmvs2(i)=zero
         kmvs3(i)=zero
      enddo
c color factors:
      ncf=3d0
      cg=ncf
      cf=(ncf**2-1d0)/ncf/2d0
c s s*
      colfacb=cg*cf/2d0
*
      s1=s+t2+u2-mt**2
      s2=s+t1+u1-mt**2
      cs1=1d0/(s1-mt**2)/s
      cs2=1d0/(s2-mt**2)/s
      call matri_qq(s,t1,t2,u1,u2,mt,mh,mat)
c      if (mt.gt.0d0) goto 999
c Form factors which multiply the SMEs:
c contributions from S1,S2,V1,V2,V3,B2:
c Higgs radiated from top: kmvs1(15)
c Higgs radiated from anti-top: kmvs2(15)
c kmvs3(15): UV+IR finite contributions (V5,B3)
c
      call qqform_s(s,t1,t2,u1,u2,mt,mh,kmvs1,kmvs2,kmvs3)
c
c interference with the Born matrix element:
      mvs=0d0
c interference of s-channel contribution with Born matrix element:
      do i=1,15
         mvs=mvs+2d0*(cs1*dreal(kmvs1(i))+
     $        cs2*dreal(kmvs2(i))+dreal(kmvs3(i)))*mat(i)
      enddo
      matv_qq=matv_qq+colfacb*mvs/9d0
 999  continue
c      if(mt.gt.0d0) goto 9999
c contribution from B3:
      colfac_b3(1)=ncf-2d0/ncf
      colfac_b3(2)=ncf-2d0/ncf
      colfac_b3(3)=-2d0/ncf
      colfac_b3(4)=-2d0/ncf
      colfac_p(1)=ncf-2d0/ncf
      colfac_p(2)=-2d0/ncf
      call qqbox_ir(s,t1,t2,u1,u2,mt,mh,kmvb3)
      do i=1,4
         matv_qq=matv_qq+2d0*dreal(kmvb3(i))*colfac_b3(i)/9d0
      enddo
c contribution from P1,P2:
      call pentagons_qq(s,t1,t2,u1,u2,mt,mh,kmvp)
      do i=1,2
         matv_qq=matv_qq+2d0*dreal(kmvp(i))*colfac_p(i)/9d0
      enddo
 9999 continue
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine qqform_s(s,t1,t2,u1,u2,mt,mh,kmv1,kmv2,kmv3)
      implicit none
      integer i,j
      real*8 s,s3,p3p4,p1p2,p1p3,p1p4,p2p3,p2p4
      real*8 t1,t2,u1,u2,mt,mh,s1,s2
      complex*16 zero,kmv1(15),kmv2(15),kmv3(15),
     $     kmvf1(6,15),kmvf2(6,15),kmvf3(2,15)
      complex*16 fvf1(6),fvf2(6),fh1(3),fh2(3),fh3(3),
     $     fs1(2),fs2(2),fs3(2),fs4(2)
      complex*16 fvb1(6),fvb2(6),kmvs(2,15)
      complex*16 sigmagg,vertexqq

      zero=dcmplx(0d0,0d0)
      do i=1,15
         kmv1(i)=zero
         kmv2(i)=zero
         do j=1,6
            kmvf1(j,i)=zero
            kmvf2(j,i)=zero
         enddo
         do j=1,2
            kmvf3(j,i)=zero
         enddo
      enddo
      s3=4d0*mt**2+mh**2-s-t1-t2-u1-u2
      s1=s+t2+u2-mt**2
      s2=s+t1+u1-mt**2
      p3p4=s/2d0
      p1p2=(s3-2d0*mt**2)/2d0
      p1p3=(mt**2-t2)/2d0
      p1p4=(mt**2-u2)/2d0
      p2p3=(mt**2-u1)/2d0
      p2p4=(mt**2-t1)/2d0

c contribution from V2(1) and V2(2):
      call formfvf(s,t1,t2,u1,u2,mt,mh,fvf1,fvf2)
      
      kmvf1(1,1)=2d0*fvf1(2)*(p1p3+p1p4-p3p4)
      kmvf1(1,2)=-fvf1(1)+2d0*mt*fvf1(2)
      kmvf1(1,3)=-fvf1(1)+2d0*mt*fvf1(2)
      kmvf1(1,6)=-2d0*fvf1(2)+fvf1(3)
      kmvf1(1,7)=-2d0*fvf1(2)+fvf1(3)
      kmvf1(1,12)=2d0*(-fvf1(1)+mt*fvf1(3)+fvf1(5)*(-p1p3-p1p4+p3p4))

      kmvf2(1,1)=2d0*fvf2(2)*(p2p3+p2p4-p3p4)
      kmvf2(1,2)=-fvf2(1)+2d0*mt*fvf2(2)
      kmvf2(1,3)=-fvf2(1)+2d0*mt*fvf2(2)
      kmvf2(1,4)=-2d0*fvf2(2)-fvf2(3)
      kmvf2(1,5)=-2d0*fvf2(2)-fvf2(3)
      kmvf2(1,13)=2d0*(fvf2(1)+mt*fvf2(3)+fvf2(5)*(p2p3+p2p4-p3p4))

c contribution from V3(1) and V3(2):
c fh3 does not contribute (only contributes to the t,u-channel)
      call formh(s,t1,t2,u1,u2,mt,mh,fh1,fh2,fh3)
      
      kmvf1(2,1)=2d0*fh1(2)*(-p1p3-p1p4+p3p4)
      kmvf1(2,2)=-fh1(1)-mt*fh1(2)-mt*fh1(3)
      kmvf1(2,3)=-fh1(1)-mt*fh1(2)-mt*fh1(3)
      kmvf1(2,12)=2d0*(-fh1(1)-mt*fh1(2)-mt*fh1(3))

      kmvf2(2,1)=2d0*fh2(3)*(-p2p3-p2p4+p3p4)
      kmvf2(2,2)=-fh2(1)-mt*fh2(2)-mt*fh2(3)
      kmvf2(2,3)=-fh2(1)-mt*fh2(2)-mt*fh2(3)
      kmvf2(2,13)=2d0*(fh2(1)+mt*fh2(2)+mt*fh2(3))

c contribution from S2(1) and S2(2):
c fs3,fs4 do not contribute (they only contribute to the t,u channel)
      call forms(s,t1,t2,u1,u2,mt,mh,fs1,fs2,fs3,fs4)

      kmvf1(3,1)=2d0*(p3p4-p1p3-p1p4)*fs1(2)
      kmvf1(3,2)=-(fs1(1)+mt*fs1(2))
      kmvf1(3,3)=-(fs1(1)+mt*fs1(2))
      kmvf1(3,12)=-2d0*(fs1(1)+mt*fs1(2))

      kmvf2(3,1)=2d0*(p3p4-p2p3-p2p4)*fs2(2)
      kmvf2(3,2)=-(fs2(1)+mt*fs2(2))
      kmvf2(3,3)=-(fs2(1)+mt*fs2(2))
      kmvf2(3,13)=2d0*(fs2(1)+mt*fs2(2))
      do i=1,13
         kmvf1(3,i)=-1d0/(s1-mt**2)*kmvf1(3,i)
         kmvf2(3,i)=-1d0/(s2-mt**2)*kmvf2(3,i)
      enddo

c contribution from S1(1) and S1(2):
c i sigma_mu_nu=i (q_mu q_nu - q^2 g_mu_nu) sigmagg 
      call ggself(s,mt,mh,sigmagg)
      kmvf1(4,2)=-(-sigmagg)
      kmvf1(4,3)=-(-sigmagg)
      kmvf1(4,12)=-2d0*(-sigmagg)

      kmvf2(4,2)=-(-sigmagg)
      kmvf2(4,3)=-(-sigmagg)
      kmvf2(4,13)=2d0*(-sigmagg)

c contribution from V1(1) and V1(2):
c -i g_s T^c Lambda_nu:
      call qqvert(s,mt,mh,vertexqq)
      kmvf1(5,2)=-vertexqq
      kmvf1(5,3)=-vertexqq
      kmvf1(5,12)=-2d0*vertexqq

      kmvf2(5,2)=-vertexqq
      kmvf2(5,3)=-vertexqq
      kmvf2(5,13)=2d0*vertexqq

c contribution from B2(1) and B2(2):
      call formb2(s,t1,t2,u1,u2,mt,mh,fvb1,fvb2)
c B2(1)
      kmvf1(6,1)=fvb1(1)
      kmvf1(6,2)=fvb1(2)
      kmvf1(6,3)=fvb1(2)
      kmvf1(6,4)=fvb1(6)
      kmvf1(6,5)=fvb1(6)
      kmvf1(6,6)=fvb1(5)
      kmvf1(6,7)=fvb1(5)
      kmvf1(6,12)=fvb1(3)
      kmvf1(6,13)=fvb1(4)
c B2(2)
      kmvf2(6,1)=fvb2(1)
      kmvf2(6,2)=fvb2(2)
      kmvf2(6,3)=fvb2(2)
      kmvf2(6,4)=fvb2(6)
      kmvf2(6,5)=fvb2(6)
      kmvf2(6,6)=fvb2(5)
      kmvf2(6,7)=fvb2(5)
      kmvf2(6,12)=fvb2(3)
      kmvf2(6,13)=fvb2(4)
      do i=1,13
         kmvf1(6,i)=kmvf1(6,i)*(s1-mt**2)
         kmvf2(6,i)=kmvf2(6,i)*(s2-mt**2)
      enddo
c UV+IR finite s-channel contributions:
      call qqs_finite(s,t1,t2,u1,u2,mt,mh,kmvs)
c contribution from B3
      do i=1,15
         kmvf3(1,i)=kmvs(1,i)/s
      enddo
c contribution from V5(1-2):
      do i=1,15
         kmvf3(2,i)=2d0*kmvs(2,i)/(2d0*mt**2+2d0*p1p2)/s
      enddo
c sum all contributions:
      do i=1,15
         do j=1,6
c         do j=1,5
            kmv1(i)=kmv1(i)+kmvf1(j,i)
            kmv2(i)=kmv2(i)+kmvf2(j,i)
         enddo
         do j=1,2
            kmv3(i)=kmv3(i)+kmvf3(j,i)
         enddo
      enddo
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine qqvert(s,mt,mh,vertex)
      implicit none
      real*8 ncf,cg,cf,nlf,mt,mh,s
      complex*16 zero,vertex
      complex*16 b0s,b00,c00
      real*8 pi,pi2,w2
      common/para/pi,pi2,w2
      real*8 delta,muedr
      common/dimreg/delta,muedr

      zero=dcmplx(0d0,0d0)
      vertex=zero
c color factors:
      ncf=3d0
      cg=ncf
      cf=(ncf**2-1d0)/ncf/2d0

c b00 and C00 are IR divergent:
c b00 is also UV divergent:
      b00=dcmplx(-dlog(mt**2/muedr**2),0d0)
      c00=1d0/s*dcmplx(dlog(mt**2/s)**2/2d0-2d0/3d0*pi**2,0d0)
      b0s=dcmplx(-dlog(s/muedr**2)+2d0,0d0)
      vertex=(cf-cg/2d0)*(b0s+4d0*(b00-b0s)-2d0*s*c00-2d0)+
     $     cg/2d0*(3d0*b0s+4d0*(b00-b0s))
c add counterterm:
      vertex=vertex-cf*b00
c check of leftover mu-dep.:
c      vertex=vertex+cg*dlog(s/muedr**2)
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine qqs_finite(s,t1,t2,u1,u2,mt,mh,kmv)
c form factors for the s-channel UV+IR finite contribution
c box diagram B3, vertex corr. V5(1-2)
      implicit none
      integer i,j,k
      real*8 ncf,cg,cf,mf,mt,mh
      real*8 mi,yuki
      real*8 s,t1,t2,u1,u2,s3,p3p4,p1p2,p1p3,p1p4,p2p3,p2p4,p1pq,p2pq
      real*8 s1,s2
      complex*16 zero,kmv(2,15),kmvi(2,15)
      complex*16 c0,d0,c11,c12,c20,c21,c22,c23
      complex*16 b0,b1,b20
      complex*16 d1(0:3),d2(0:6),d3(0:3,0:3),c1(0:2,4),c2(0:3,4)
      complex*16 d11,d12,d13,d21,d22,d23,d20,d212,d213,d223
      complex*16 d301,d302,d303,d31,d32,d33,d312,d313,d321,d323,
     $     d331,d332,d3123

      real*8 yukt,yukb,yukf
      common/yukawa/yukt,yukb,yukf
      real*8 mtop,alphas,mb
      common/topas/mtop,alphas,mb

      zero=dcmplx(0d0,0d0)
      do i=1,15
         do j=1,2
            kmv(j,i)=zero
            kmvi(j,i)=zero
         enddo
      enddo
      s3=4d0*mt**2+mh**2-s-t1-t2-u1-u2
      s1=s+t2+u2-mt**2
      s2=s+t1+u1-mt**2
      p3p4=s/2d0
      p1p2=(s3-2d0*mt**2)/2d0
      p1p3=(mt**2-t2)/2d0
      p1p4=(mt**2-u2)/2d0
      p2p3=(mt**2-u1)/2d0
      p2p4=(mt**2-t1)/2d0
      p1pq=p1p3+p1p4
      p2pq=p2p3+p2p4
c color factors:
      ncf=3d0
      cg=ncf
      cf=(ncf**2-1d0)/ncf/2d0
c top and bottom loop in V5
      mf=mt
      do k=1,2
         if(k.eq.1)then
            mi=mtop
            yuki=yukt
         endif
         if(k.eq.2)then
            mi=mb
            yuki=yukb
         endif
c V5(1),V5(2)
         call ccfunc2(2d0*mf**2+2d0*p1p2,mh**2,
     $        2d0*mf**2+2d0*p1p2-p1p4-p2p4-p1p3-p2p3,mi,mi,mi,
     $        c0,c11,c12,c21,c22,c20,c23)
         call bfunc(s,mi,mi,b0,b1,b20)
         kmvi(2,1)=(-4*b0+16*c20-8*c0*mf**2-16*c11*mf**2-16*c12*mf**2- 
     $        8*c0*p1p2-16*c11*p1p2-16*c12*p1p2+4*c0*p1p3+ 
     $        8*c12*p1p3+4*c0*p1p4+8*c12*p1p4+4*c0*p2p3+ 
     $        8*c12*p2p3+4*c0*p2p4+8*c12*p2p4)
         kmvi(2,4)=(-4*c0 - 16*c12 - 16*c22 - 16*c23)
         kmvi(2,5)=(-4*c0 - 16*c12 - 16*c22 - 16*c23)
         kmvi(2,6)=(-4*c0 - 16*c12 - 16*c22 - 16*c23)
         kmvi(2,7)=(-4*c0 - 16*c12 - 16*c22 - 16*c23)
c 1/2 comes from color factor
         do i=1,15
            kmv(2,i)=kmv(2,i)-kmvi(2,i)*mi/2d0*(yuki/yukf)
         enddo
      enddo
c B3
      call ddfunc(s,s+mt**2-2d0*p1pq,mt**2,s-p1pq,p2pq,p2pq-p1p2,
     $      0d0,0d0,mt,mt,d1,d2,d3,c1,c2)
      c0=c1(0,1)
      c11=c1(1,1)
      c12=c1(2,1)
      d0=d1(0)
      d11=d1(1)
      d12=d1(2)
      d13=d1(3)
      d20=d2(0)
      d21=d2(1)
      d22=d2(2)
      d23=d2(3)
      d212=d2(4)
      d213=d2(5)
      d223=d2(6)
      kmv(1,1)=(-12*d20*mt-2*d22*mt**3+4*d223*mt**3-2*d23*mt**3-
     $     2*d22*mt*p1p2+4*d223*mt*p1p2-2*d23*mt*p1p2-
     $     2*d12*mt*p1p3+2*d13*mt*p1p3+2*d212*mt*p1p3-
     $     2*d213*mt*p1p3+4*d22*mt*p1p3-4*d223*mt*p1p3-
     $     2*d12*mt*p1p4+2*d13*mt*p1p4+2*d212*mt*p1p4-
     $     2*d213*mt*p1p4+4*d22*mt*p1p4-4*d223*mt*p1p4+
     $     4*d12*mt*p2p3-4*d13*mt*p2p3+2*d212*mt*p2p3-
     $     2*d213*mt*p2p3+2*d22*mt*p2p3-4*d223*mt*p2p3+
     $     2*d23*mt*p2p3+4*d12*mt*p2p4-4*d13*mt*p2p4+
     $     2*d212*mt*p2p4-2*d213*mt*p2p4+2*d22*mt*p2p4-
     $     4*d223*mt*p2p4+2*d23*mt*p2p4-2*d12*mt*p3p4+
     $     2*d13*mt*p3p4-4*d212*mt*p3p4+4*d213*mt*p3p4-
     $     4*d22*mt*p3p4+4*d223*mt*p3p4)
      kmv(1,2)=(5*c0-2*d20+2*d13*mt**2-2*d212*mt**2+2*d213*mt**2-
     $     2*d22*mt**2+2*d223*mt**2-2*d12*p1p2-2*d212*p1p2+
     $     2*d213*p1p2-2*d22*p1p2+2*d223*p1p2-2*d11*p1p3-
     $     6*d12*p1p3+2*d21*p1p3+6*d212*p1p3+4*d22*p1p3-
     $     2*d11*p1p4-6*d12*p1p4+2*d21*p1p4+6*d212*p1p4+
     $     4*d22*p1p4+4*d0*p2p3+6*d11*p2p3+6*d12*p2p3+
     $     4*d13*p2p3+2*d21*p2p3+4*d212*p2p3-2*d213*p2p3+
     $     2*d22*p2p3-2*d223*p2p3+4*d0*p2p4+6*d11*p2p4+
     $     6*d12*p2p4+4*d13*p2p4+2*d21*p2p4+4*d212*p2p4-
     $     2*d213*p2p4+2*d22*p2p4-2*d223*p2p4+6*d11*p3p4+
     $     6*d12*p3p4-4*d21*p3p4-8*d212*p3p4-4*d22*p3p4)
      kmv(1,3)=(5*c0-2*d20+2*d13*mt**2-2*d212*mt**2+2*d213*mt**2-
     $     2*d22*mt**2+2*d223*mt**2-2*d12*p1p2-2*d212*p1p2+
     $     2*d213*p1p2-2*d22*p1p2+2*d223*p1p2-2*d11*p1p3-
     $     6*d12*p1p3+2*d21*p1p3+6*d212*p1p3+4*d22*p1p3-
     $     2*d11*p1p4-6*d12*p1p4+2*d21*p1p4+6*d212*p1p4+
     $     4*d22*p1p4+4*d0*p2p3+6*d11*p2p3+6*d12*p2p3+
     $     4*d13*p2p3+2*d21*p2p3+4*d212*p2p3-2*d213*p2p3+
     $     2*d22*p2p3-2*d223*p2p3+4*d0*p2p4+6*d11*p2p4+
     $     6*d12*p2p4+4*d13*p2p4+2*d21*p2p4+4*d212*p2p4-
     $     2*d213*p2p4+2*d22*p2p4-2*d223*p2p4+6*d11*p3p4+
     $     6*d12*p3p4-4*d21*p3p4-8*d212*p3p4-4*d22*p3p4)
      kmv(1,4)=(-4*d12*mt-8*d213*mt-8*d223*mt)
      kmv(1,5)=(-4*d12*mt-8*d213*mt-8*d223*mt)
      kmv(1,6)=(8*d12*mt-4*d13*mt+8*d212*mt+8*d22*mt)
      kmv(1,7)=(8*d12*mt-4*d13*mt+8*d212*mt+8*d22*mt)
      kmv(1,12)=(2*c0-6*c11-8*d20+4*d13*mt**2+4*d223*mt**2+
     $     4*d12*p1p2+4*d22*p1p2+8*d223*p1p2+8*d212*p1p3+
     $     16*d22*p1p3+8*d212*p1p4+16*d22*p1p4+8*d0*p2p3+
     $     4*d11*p2p3+4*d12*p2p3+8*d13*p2p3-4*d212*p2p3-
     $     4*d22*p2p3-8*d223*p2p3+8*d0*p2p4+4*d11*p2p4+
     $     4*d12*p2p4+8*d13*p2p4-4*d212*p2p4-4*d22*p2p4-
     $     8*d223*p2p4+8*d11*p3p4-16*d212*p3p4-16*d22*p3p4)
      kmv(1,13)=(-2*c0+6*c12+8*d20-4*d12*mt**2-4*d223*mt**2-
     $     4*d13*p1p2-8*d223*p1p2-4*d23*p1p2-4*d0*p1p3+
     $     4*d11*p1p3-4*d13*p1p3-4*d213*p1p3-8*d223*p1p3-
     $     4*d0*p1p4+4*d11*p1p4-4*d13*p1p4-4*d213*p1p4-
     $     8*d223*p1p4+4*d13*p2p3+8*d213*p2p3+8*d223*p2p3+
     $     4*d23*p2p3+4*d13*p2p4+8*d213*p2p4+8*d223*p2p4+
     $     4*d23*p2p4+4*d0*p3p4+4*d13*p3p4+8*d213*p3p4+
     $     8*d223*p3p4)
      do i=1,15
         kmv(1,i)=-cg/2d0*kmv(1,i)
      enddo
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine matri_qq(s,t1,t2,u1,u2,mt,mh,mat)
c----------------------------------------------------------------------
c mat(i) = (summe ueber quark spins)*trace(
c            born-matrixelem.*(p1s-mt)*i-entwicklungskoeff.*(p2+mt))/4
c  i: 1-15 SMEs
c----------------------------------------------------------------------
      implicit none
      integer i
      real*8 s,s1,s2,s3,t1,t2,u1,u2,mt,mh
      real*8 cs1,cs2
      real*8 p1p2,p1p3,p1p4,p2p3,p2p4
      real*8 mat(15)
c initialization
      do i=1,15
         mat(i)=0d0
      end do
      s3=4d0*mt**2+mh**2-s-t1-t2-u1-u2
      p1p2=(s3-2d0*mt**2)/2d0
      p1p3=(mt**2-t2)/2d0
      p1p4=(mt**2-u2)/2d0
      p2p3=(mt**2-u1)/2d0
      p2p4=(mt**2-t1)/2d0
      s1=s+t2+u2-mt**2
      s2=s+t1+u1-mt**2
      cs1=1d0/(s1-mt**2)/s
      cs2=1d0/(s2-mt**2)/s

      mat(1)=(cs1*(16*mt**3*s-16*mt*s*p1p2+16*mt*s*p1p3+
     $     16*mt*s*p1p4-64*mt*p1p3*p1p4+16*mt*s*p2p3+
     $     32*mt*p1p4*p2p3+16*mt*s*p2p4+32*mt*p1p3*p2p4)+ 
     $     cs2*(16*mt**3*s-16*mt*s*p1p2+16*mt*s*p1p3+16*mt*s*p1p4+
     $     16*mt*s*p2p3+32*mt*p1p4*p2p3+16*mt*s*p2p4+
     $     32*mt*p1p3*p2p4-64*mt*p2p3*p2p4))
      mat(2)=(cs1*(-8*mt**2*s**2+16*s*p1p2*p1p4-
     $     32*p1p4**2*p2p3-16*mt**2*s*p2p4-
     $     32*s*p1p4*p2p4+32*p1p3*p1p4*p2p4)+
     $     cs2*(-8*mt**2*s**2-16*mt**2*s*p1p4+16*s*p1p2*p2p4-
     $     32*s*p1p4*p2p4+32*p1p4*p2p3*p2p4-
     $     32*p1p3*p2p4**2))
      mat(3)=(cs1*(-8*mt**2*s**2+16*s*p1p2*p1p3-
     $     16*mt**2*s*p2p3-32*s*p1p3*p2p3+
     $     32*p1p3*p1p4*p2p3-32*p1p3**2*p2p4)+
     $     cs2*(-8*mt**2*s**2-16*mt**2*s*p1p3+16*s*p1p2*p2p3-
     $     32*s*p1p3*p2p3-32*p1p4*p2p3**2+
     $     32*p1p3*p2p3*p2p4))
      mat(4)=(cs1*(-4*mt**3*s**2-4*mt*s**2*p1p2+
     $     16*mt*s*p1p2*p1p4+8*mt*s*p1p4*p2p3-
     $     32*mt*p1p4**2*p2p3-16*mt*s*p1p2*p2p4+
     $     8*mt*s*p1p3*p2p4-32*mt*p1p3*p1p4*p2p4+
     $     16*mt*s*p2p3*p2p4+32*mt*p1p4*p2p3*p2p4+
     $     32*mt*p1p3*p2p4**2)+
     $     cs2*(-4*mt**3*s**2-4*mt*s**2*p1p2-16*mt**3*s*p1p4+
     $     8*mt*s*p1p4*p2p3+16*mt**3*s*p2p4+
     $     8*mt*s*p1p3*p2p4+16*mt*s*p2p3*p2p4+
     $     64*mt*p1p4*p2p3*p2p4-64*mt*p2p3*p2p4**2))
      mat(5)=(cs1*(-4*mt**3*s**2-4*mt*s**2*p1p2+
     $     16*mt*s*p1p2*p1p3-16*mt*s*p1p2*p2p3+
     $     8*mt*s*p1p4*p2p3-32*mt*p1p3*p1p4*p2p3+
     $     32*mt*p1p4*p2p3**2+8*mt*s*p1p3*p2p4-
     $     32*mt*p1p3**2*p2p4+16*mt*s*p2p3*p2p4+
     $     32*mt*p1p3*p2p3*p2p4)+
     $     cs2*(-4*mt**3*s**2-4*mt*s**2*p1p2-16*mt**3*s*p1p3+
     $     16*mt**3*s*p2p3+8*mt*s*p1p4*p2p3+
     $     8*mt*s*p1p3*p2p4+16*mt*s*p2p3*p2p4+
     $     64*mt*p1p3*p2p3*p2p4-64*mt*p2p3**2*p2p4))
      mat(6)=(cs1*(-4*mt**3*s**2-4*mt*s**2*p1p2+16*mt**3*s*p1p4+
     $     16*mt*s*p1p3*p1p4-64*mt*p1p3*p1p4**2+
     $     8*mt*s*p1p4*p2p3-16*mt**3*s*p2p4+
     $     8*mt*s*p1p3*p2p4+64*mt*p1p3*p1p4*p2p4)+
     $     cs2*(-4*mt**3*s**2-4*mt*s**2*p1p2-16*mt*s*p1p2*p1p4+
     $     16*mt*s*p1p3*p1p4+8*mt*s*p1p4*p2p3+
     $     32*mt*p1p4**2*p2p3+16*mt*s*p1p2*p2p4+
     $     8*mt*s*p1p3*p2p4+32*mt*p1p3*p1p4*p2p4-
     $     32*mt*p1p4*p2p3*p2p4-32*mt*p1p3*p2p4**2))
      mat(7)=(cs1*(-4*mt**3*s**2-4*mt*s**2*p1p2+16*mt**3*s*p1p3+
     $     16*mt*s*p1p3*p1p4-64*mt*p1p3**2*p1p4-
     $     16*mt**3*s*p2p3+8*mt*s*p1p4*p2p3+
     $     64*mt*p1p3*p1p4*p2p3+8*mt*s*p1p3*p2p4)+
     $     cs2*(-4*mt**3*s**2-4*mt*s**2*p1p2-16*mt*s*p1p2*p1p3+
     $     16*mt*s*p1p3*p1p4+16*mt*s*p1p2*p2p3+
     $     8*mt*s*p1p4*p2p3+32*mt*p1p3*p1p4*p2p3-
     $     32*mt*p1p4*p2p3**2+8*mt*s*p1p3*p2p4+
     $     32*mt*p1p3**2*p2p4-32*mt*p1p3*p2p3*p2p4))
      do i=8,11
         mat(i)=0d0
      enddo
      mat(12)=(cs1*(-16*mt**4*s+16*mt**2*s*p1p2+8*s*p1p2*p1p3+
     $     8*s*p1p2*p1p4+64*mt**2*p1p3*p1p4-
     $     64*p1p2*p1p3*p1p4-8*mt**2*s*p2p3+
     $     16*p1p3*p1p4*p2p3-16*p1p4**2*p2p3-
     $     8*mt**2*s*p2p4-16*p1p3**2*p2p4+
     $     16*p1p3*p1p4*p2p4)+
     $     cs2*(16*mt**2*s*p1p2-16*s*p1p2**2+8*s*p1p2*p1p3+
     $     8*s*p1p2*p1p4-8*mt**2*s*p2p3-
     $     32*mt**2*p1p4*p2p3+32*p1p2*p1p4*p2p3+
     $     16*p1p3*p1p4*p2p3-16*p1p4**2*p2p3-
     $     8*mt**2*s*p2p4-32*mt**2*p1p3*p2p4+
     $     32*p1p2*p1p3*p2p4-16*p1p3**2*p2p4+
     $     16*p1p3*p1p4*p2p4))
      mat(13)=(cs1*(-16*mt**2*s*p1p2+16*s*p1p2**2+
     $     8*mt**2*s*p1p3+8*mt**2*s*p1p4-8*s*p1p2*p2p3+
     $     32*mt**2*p1p4*p2p3-32*p1p2*p1p4*p2p3+
     $     16*p1p4*p2p3**2-8*s*p1p2*p2p4+
     $     32*mt**2*p1p3*p2p4-32*p1p2*p1p3*p2p4-
     $     16*p1p3*p2p3*p2p4-16*p1p4*p2p3*p2p4+
     $     16*p1p3*p2p4**2)+
     $     cs2*(16*mt**4*s-16*mt**2*s*p1p2+8*mt**2*s*p1p3+
     $     8*mt**2*s*p1p4-8*s*p1p2*p2p3+
     $     16*p1p4*p2p3**2-8*s*p1p2*p2p4-
     $     64*mt**2*p2p3*p2p4+64*p1p2*p2p3*p2p4-
     $     16*p1p3*p2p3*p2p4-16*p1p4*p2p3*p2p4+
     $     16*p1p3*p2p4**2))
      do i=14,15
         mat(i)=0d0
      enddo
      do i=1,15
         mat(i)=mat(i)/4d0
      end do
      return
      end
*      
      double precision function matbox(s,t1,t2,u1,u2,mt,mh)
      implicit none
      integer i
      real*8 pi,s,t1,t2,u1,u2,mt,mh
      real*8 ncf,cg,cf
      real*8 colfac,mat(25),mbox(25)

      matbox=0d0
      pi=4d0*datan(1d0)
c color factors:
      ncf=3d0
      cg=ncf
      cf=(ncf**2-1d0)/ncf/2d0
c color factor:
      colfac=cf*cg/2d0*(cf-cg/2d0)

c interference of Box with Born matrix element:
      call matribox(s,t1,t2,u1,u2,mt,mh,mat)
c form factors :
      call box31(s,t1,t2,u1,u2,mt,mh,mbox)

      do i=1,25
         matbox=matbox+2d0*mbox(i)*mat(i)
      end do
      matbox=colfac/9d0*matbox
      return
      end

      subroutine box31(s,t1,t2,u1,u2,mt,mh,mbox)
      implicit none
      integer i
      real*8 s,s3,t1,t2,u1,u2,mt,mh
      real*8 p1p2,p1p3,p1p4,p2p3,p2p4,s1,s2
      real*8 d1(0:3),d2(0:6),d3(0:3,0:3),ch1(0:2),mbox(25)
      real*8 d0,d11,d12,d13,d20,d21,d22,d23,d212,d213,d223
      real*8 c12,c0,det
      real*8 lambda,lambda2
      complex*16 b0s,b0h,b0s3,b1

      COMMON/IR/LAMBDA,LAMBDA2

c initialization
      do i=1,25
         mbox(i)=0d0
      end do
      s3=4d0*mt**2+mh**2-s-t1-t2-u1-u2
      p1p2=(s3-2d0*mt**2)/2d0
      p1p3=(mt**2-t2)/2d0
      p1p4=(mt**2-u2)/2d0
      p2p3=(mt**2-u1)/2d0
      p2p4=(mt**2-t1)/2d0
      s1=s+t2+u2-mt**2
      s2=s+t1+u1-mt**2
c calculation of the coefficients
      call dmn(mt**2,mh**2,s,s1,s2,mt**2,0d0,
     $     mt,mt,mt,d1,d2,d3,ch1)
      det=((s-s3)**2+mh**4-2d0*mh**2*(s+s3))/4d0
      c0 = ch1(0)
      d0 = d1(0)
      call bquer(s,mt,mt,b0s,b1)
      call bquer(mh**2,mt,mt,b0h,b1)
      call bquer(s3,mt,mt,b0s3,b1)
      c12 = 1d0/det/2d0*dreal(-s*b0s+
     $     (s-s3+mh**2)/2d0*b0h+
     $     (s+s3-mh**2)/2d0*b0s3+s*(s3-s+mh**2)/2d0*c0)
c      write(6,*)'d0=',d0
c      write(6,*)'c01=',c0
c      write(6,*)'b0s=',b0s
c      write(6,*)'b0h=',b0h
c      write(6,*)'b0s3=',b0s3
      d11 = d1(1)
      d12 = d1(2)
      d13 = d1(3)
      d20 = d2(0)
      d21 = d2(1)
      d22 = d2(2)
      d23 = d2(3)
      d212 = d2(4)
      d213 = d2(5)
      d223 = d2(6)
      mbox(1)=(4*p1p2*d0 + 4*d11*p1p2 + 4*d12*p1p2 + 4*d13*p1p2)
      mbox(2)=(2*c0*mt - 8*d20*mt + d11*(-4*mt*p1p2 + 4*mt*p1p3) + d13
     $     *(4*mt**3 + 4*mt*p2p3) + d12*(4*mt**3 - 4*mt*p1p3 + 4*mt*p2p3
     $     ))
      mbox(3)=(-4*d11*mt + 8*d12*mt + 8*d22*mt + 8*d223*mt)
      mbox(4)=(-4*d12*mt - 4*d13*mt - 8*d212*mt)
      mbox(5)=(-4*c0 + 4*c12 + 8*p1p2*d0 + 16*d13*p1p2 + 
     $     d23*(8*p1p2 + 8*mt**2) + d22*(8*p1p2 - 8*p2p3 + 8*mt**2) + 
     $     d223*(16*p1p2 - 8*p2p3 + 16*mt**2) + 
     $     d11*8*p1p2 + d12*(8*p1p3 - 4*s + 16*p1p2 - 8*p2p3))
      mbox(6)=(-4*c0 + 4*c12 + 8*d20 - 8*d13*mt**2 + 
     $     d213*(-8*p1p2 - 8*mt**2) + 
     $     d212*(-8*p1p2 + 8*p2p3 - 8*mt**2) + 
     $     d12*(-8*mt**2 - 4*s + 8*p1p3))
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine matribox(s,t1,t2,u1,u2,mt,mh,mat)
c----------------------------------------------------------------------
c mat(i,j) = trace({mu}*p4*{nu}*p3)*trace(
c            born-matrixelem.*(p1s-mt)*i-entwicklungskoeff.*(p2+mt))/4
c i: 1-25 SMEs
c----------------------------------------------------------------------
      implicit none
      integer i,j
      real*8 s,s1,s2,s3,t1,t2,u1,u2,mt,mh,cs1,cs2
      real*8 p1p2,p1p3,p1p4,p2p3,p2p4
      real*8 mat(25)
c initialization
      do i=1,25
         mat(i)=0d0
      end do
      s3=4d0*mt**2+mh**2-s-t1-t2-u1-u2
      p1p2=(s3-2d0*mt**2)/2d0
      p1p3=(mt**2-t2)/2d0
      p1p4=(mt**2-u2)/2d0
      p2p3=(mt**2-u1)/2d0
      p2p4=(mt**2-t1)/2d0
      s1=s+t2+u2-mt**2
      s2=s+t1+u1-mt**2
      cs1=1d0/(s1-mt**2)
      cs2=1d0/(s2-mt**2)
*
c G[{ro},p3]:
      mat(1)=cs1*(-16*mt**2*s**2 + 16*s*p1p2*p1p3 + 
     $        16*s*p1p2*p1p4 - 16*mt**2*s*p2p3 - 
     $        32*s*p1p3*p2p3 + 32*p1p3*p1p4*p2p3 - 
     $        32*p1p4**2*p2p3 - 16*mt**2*s*p2p4 - 
     $        32*p1p3**2*p2p4 - 32*s*p1p4*p2p4 + 
     $        32*p1p3*p1p4*p2p4) + 
     $        cs2*(-16*mt**2*s**2 - 16*mt**2*s*p1p3 -16*mt**2*s*p1p4+ 
     $        16*s*p1p2*p2p3 - 32*s*p1p3*p2p3 - 
     $        32*p1p4*p2p3**2 + 16*s*p1p2*p2p4 - 
     $        32*s*p1p4*p2p4 + 32*p1p3*p2p3*p2p4 + 
     $        32*p1p4*p2p3*p2p4 - 32*p1p3*p2p4**2)
c G[{ro}]:
      mat(2)=cs1*(16*mt**3*s - 16*mt*s*p1p2 + 16*mt*s*p1p3 + 
     $     16*mt*s*p1p4 - 64*mt*p1p3*p1p4 + 16*mt*s*p2p3 + 
     $     32*mt*p1p4*p2p3 + 16*mt*s*p2p4 + 
     $     32*mt*p1p3*p2p4) + 
     $     cs2*(16*mt**3*s - 16*mt*s*p1p2 + 16*mt*s*p1p3 + 
     $     16*mt*s*p1p4 + 16*mt*s*p2p3 + 32*mt*p1p4*p2p3 + 
     $     16*mt*s*p2p4 + 32*mt*p1p3*p2p4 - 64*mt*p2p3*p2p4)
c G[p3] p1.{ro}:
      mat(3)=cs1*(-8*mt**3*s**2 - 8*mt*s**2*p1p2 + 
     $     16*mt**3*s*p1p3 + 16*mt**3*s*p1p4 + 
     $     32*mt*s*p1p3*p1p4 - 64*mt*p1p3**2*p1p4 - 
     $     64*mt*p1p3*p1p4**2 - 16*mt**3*s*p2p3 + 
     $     16*mt*s*p1p4*p2p3 + 64*mt*p1p3*p1p4*p2p3 - 
     $     16*mt**3*s*p2p4 + 16*mt*s*p1p3*p2p4 + 
     $     64*mt*p1p3*p1p4*p2p4) + 
     $     cs2*(-8*mt**3*s**2 - 8*mt*s**2*p1p2 - 16*mt*s*p1p2*p1p3 - 
     $     16*mt*s*p1p2*p1p4 + 32*mt*s*p1p3*p1p4 + 
     $     16*mt*s*p1p2*p2p3 + 16*mt*s*p1p4*p2p3 + 
     $     32*mt*p1p3*p1p4*p2p3 + 32*mt*p1p4**2*p2p3 - 
     $     32*mt*p1p4*p2p3**2 + 16*mt*s*p1p2*p2p4 + 
     $     16*mt*s*p1p3*p2p4 + 32*mt*p1p3**2*p2p4 + 
     $     32*mt*p1p3*p1p4*p2p4 - 
     $     32*mt*p1p3*p2p3*p2p4 - 
     $     32*mt*p1p4*p2p3*p2p4 - 32*mt*p1p3*p2p4**2)
c G[p3] p2.{ro}:
      mat(4)=cs1*(-8*mt**3*s**2 - 8*mt*s**2*p1p2 + 
     $     16*mt*s*p1p2*p1p3 + 16*mt*s*p1p2*p1p4 - 
     $     16*mt*s*p1p2*p2p3 + 16*mt*s*p1p4*p2p3 - 
     $     32*mt*p1p3*p1p4*p2p3 - 32*mt*p1p4**2*p2p3 + 
     $     32*mt*p1p4*p2p3**2 - 16*mt*s*p1p2*p2p4 + 
     $     16*mt*s*p1p3*p2p4 - 32*mt*p1p3**2*p2p4 - 
     $     32*mt*p1p3*p1p4*p2p4 + 32*mt*s*p2p3*p2p4 + 
     $     32*mt*p1p3*p2p3*p2p4 + 
     $     32*mt*p1p4*p2p3*p2p4 + 32*mt*p1p3*p2p4**2) + 
     $     cs2*(-8*mt**3*s**2 - 8*mt*s**2*p1p2 - 16*mt**3*s*p1p3 - 
     $     16*mt**3*s*p1p4 + 16*mt**3*s*p2p3 + 
     $     16*mt*s*p1p4*p2p3 + 16*mt**3*s*p2p4 + 
     $     16*mt*s*p1p3*p2p4 + 32*mt*s*p2p3*p2p4 + 
     $     64*mt*p1p3*p2p3*p2p4 + 
     $     64*mt*p1p4*p2p3*p2p4 - 64*mt*p2p3**2*p2p4 - 
     $     64*mt*p2p3*p2p4**2)
c p1.{ro}:
      mat(5)=cs1*(-16*mt**4*s + 16*mt**2*s*p1p2 + 8*s*p1p2*p1p3 + 
     $     8*s*p1p2*p1p4 + 64*mt**2*p1p3*p1p4 - 
     $     64*p1p2*p1p3*p1p4 - 8*mt**2*s*p2p3 + 
     $     16*p1p3*p1p4*p2p3 - 16*p1p4**2*p2p3 - 
     $     8*mt**2*s*p2p4 - 16*p1p3**2*p2p4 + 
     $     16*p1p3*p1p4*p2p4) + 
     $     cs2*(16*mt**2*s*p1p2 - 16*s*p1p2**2 + 8*s*p1p2*p1p3 + 
     $     8*s*p1p2*p1p4 - 8*mt**2*s*p2p3 - 
     $     32*mt**2*p1p4*p2p3 + 32*p1p2*p1p4*p2p3 + 
     $     16*p1p3*p1p4*p2p3 - 16*p1p4**2*p2p3 - 
     $     8*mt**2*s*p2p4 - 32*mt**2*p1p3*p2p4 + 
     $     32*p1p2*p1p3*p2p4 - 16*p1p3**2*p2p4 + 
     $     16*p1p3*p1p4*p2p4)
c p2.{ro}:
      mat(6)=cs1*(-16*mt**2*s*p1p2 + 16*s*p1p2**2 + 
     $     8*mt**2*s*p1p3 + 8*mt**2*s*p1p4 - 8*s*p1p2*p2p3 + 
     $     32*mt**2*p1p4*p2p3 - 32*p1p2*p1p4*p2p3 + 
     $     16*p1p4*p2p3**2 - 8*s*p1p2*p2p4 + 
     $     32*mt**2*p1p3*p2p4 - 32*p1p2*p1p3*p2p4 - 
     $     16*p1p3*p2p3*p2p4 - 16*p1p4*p2p3*p2p4 + 
     $     16*p1p3*p2p4**2) + 
     $     cs2*(16*mt**4*s - 16*mt**2*s*p1p2 + 8*mt**2*s*p1p3 + 
     $     8*mt**2*s*p1p4 - 8*s*p1p2*p2p3 + 
     $     16*p1p4*p2p3**2 - 8*s*p1p2*p2p4 - 
     $     64*mt**2*p2p3*p2p4 + 64*p1p2*p2p3*p2p4 - 
     $     16*p1p3*p2p3*p2p4 - 16*p1p4*p2p3*p2p4 + 
     $     16*p1p3*p2p4**2)

      do i=1,25
         mat(i)=mat(i)/4d0
      end do
      return
      end
************************************************************************
        function C0_(q01,q12,q20,m0,m1,m2,ext)                            
************************************************************************
*	scalar 3-point function
*-----------------------------------------------------------------------
*	General result from A.Denner, Fortschr. Phys. 41 (1993) 307
* d0=q^2-m0^2, d1=(q+p1)^2-m1^2, d2=(q+p2)^2-m2^2
* q01=p1^2, q20=p2^2, q12=(p1-p2)^2
*-----------------------------------------------------------------------
************************************************************************
        implicit real*8 (a-z)                                         
	complex*16 c0_,cspens,etass,ieps,alpha,alp(0:2),sc,xs
	complex*16 y0(0:2),y(0:2,-1:1),x(0:2,-1:1)
	real*8 thp(0:2),thy(0:2)
	integer i,j,ext

	COMMON/IR/LAMBDA,LAMBDA2

	kappa2(a,b,c) = a**2+b**2+c**2-2d0*a*b-2d0*a*c-2d0*b*c
        pi   = 4d0*datan(1d0)                                               
	eps  = 1d-17
	ieps = dcmplx(0d0,eps)
	m02 = m0**2
	m12 = m1**2
	m22 = m2**2
	p01 = q01
	p12 = q12
	p20 = q20
C*** Check for IR divergence
c    Permutate until m0=0
10	continue
	if ((m02*m12*m22.eq.0d0).and.(m02.ne.0d0)) then
	  m   = m02
	  m02 = m12
	  m12 = m22
	  m22 = m
	  p   = p01
	  p01 = p12
	  p12 = p20
	  p20 = p
	else
	  goto 20
	endif
	goto 10  
20	continue
	if ((p01.eq.m12).and.(p20.eq.m22).and.(m02.eq.0d0)) goto 500

C****** Regular C0 function
c    Permutate until p01=0
11	continue
	if ((p01*p12*p20.eq.0d0).and.(p01.ne.0d0)) then
	  m   = m02
	  m02 = m12
	  m12 = m22
	  m22 = m
	  p   = p01
	  p01 = p12
	  p12 = p20
	  p20 = p
	else
	  goto 21
	endif
	goto 11  
21	continue
	if (p01.eq.0d0) goto 600

C*** Regular C0 function with p01,p12,p20 =/= 0
	alpha  = sqrt( abs(kappa2(p01,p12,p20)) )
	alp(0) = sqrt( kappa2(p12,m12,m22)*(1d0+ieps*sign(1d0,p12)) )
	alp(1) = sqrt( kappa2(p20,m22,m02)*(1d0+ieps*sign(1d0,p20)) )
	alp(2) = sqrt( kappa2(p01,m02,m12)*(1d0+ieps*sign(1d0,p01)) )

	do 99 i=0,2
	  if (alp(i).eq.dcmplx(0d0,0d0)) alp(i) = ieps*abs(alpha)
99	continue

	y0(0)  = ( p12*(p12-p01-p20+2d0*m02-m12-m22)
     &	  - (p20-p01)*(m12-m22)+alpha*(p12-m12+m22) )/2d0/alpha/p12
	y0(1)  = ( p20*(p20-p12-p01+2d0*m12-m22-m02)
     &	  - (p01-p12)*(m22-m02)+alpha*(p20-m22+m02) )/2d0/alpha/p20
	y0(2)  = ( p01*(p01-p20-p12+2d0*m22-m02-m12)
     &	  - (p12-p20)*(m02-m12)+alpha*(p01-m02+m12) )/2d0/alpha/p01

	do 100 j=-1,1,2
	  x(0,j) = (p12-m12+m22+j*alp(0))/2d0/p12
	  x(1,j) = (p20-m22+m02+j*alp(1))/2d0/p20
	  x(2,j) = (p01-m02+m12+j*alp(2))/2d0/p01
	do 100 i=0,2
	  y(i,j) = y0(i)-x(i,j)
100	continue

	do 200 i=0,2
	  thp(i) = 0d0
	  thy(i) = 0d0
	  if (dimag(y(i,+1)*y(i,-1)).le.0d0) thy(i) = 1d0
200	continue
	if (p12.le.0d0) thp(0) = 1d0
	if (p20.le.0d0) thp(1) = 1d0
	if (p01.le.0d0) thp(2) = 1d0

	c0_ = 0d0
	do 300 i=0,2
	do 400 j=-1,1,2
	  c0_ = c0_ + cspens((y0(i)-1d0)/y(i,j)) - cspens(y0(i)/y(i,j))
     &	         + etass(1d0-x(i,j),1d0/y(i,j))*log((y0(i)-1d0)/y(i,j))
     &	            - etass(   -x(i,j),1d0/y(i,j))*log(y0(i)/y(i,j))
400	continue
	  c0_ = c0_ - log((y0(i)-1d0)/y0(i))*(
     &		        etass(-x(i,+1),-x(i,-1))-etass(y(i,+1),y(i,-1))
     &		      - dcmplx(0d0,2d0*pi)*thp(i)*thy(i) )
300	continue
	c0_ = c0_/alpha

	return
600	continue

	if (m02*m12.eq.0d0) goto 700
C*** Regular C0 function with p01=0 and m0,m1=/=0
	alpha   = 1d0+(p12-p20)/(m02-m12-ieps*(m02+m12))
	alp(0)  = sqrt( kappa2(p20,m02,m22)+ieps*p20*(p20-m02-m22) )
	alp(1)  = sqrt( kappa2(p12,m12,m22)+ieps*p12*(p12-m12-m22) )
	x(0,+1) = (p20-m02-m22+alp(0))/2d0/m02
	x(0,-1) = (p20-m02-m22-alp(0))/2d0/m02
	x(1,+1) = (p12-m12-m22+alp(1))/2d0/m12
	x(1,-1) = (p12-m12-m22-alp(1))/2d0/m12
	c0_ = 0d0
	do 610 i=-1,1,2
	do 620 j=0,1
	c0_ = c0_ + (1d0-2d0*j)*( 
     &		cspens(1d0+x(j,i)/alpha) - cspens(1d0+x(j,i)) 
     &		+ etass(1d0/alpha,-x(j,i))*log(1d0+x(j,i)/alpha) )
620	continue
610	continue
	c0_ = c0_ + log(alpha)*log(m02/m12)
	c0_ = c0_/(p20-p12) 

	return
700	continue

C*** Regular C0 function with p01=0 and m0=0
	if (m12.eq.0d0) then
	  m12 = m02
	  m02 = 0d0
	  p   = p12
	  p12 = p20
	  p20 = p
	endif
	alpha   = 1d0+(p12-p20)/(-m12-ieps*m12)
	alp(1)  = sqrt( kappa2(p12,m12,m22)+ieps*p12*(p12-m12-m22) )
	x(0,-1) = m22/(p20-m22+ieps*m12)
	x(1,+1) = (p12-m12-m22+alp(1))/2d0/m12
	x(1,-1) = (p12-m12-m22-alp(1))/2d0/m12
	c0_ = 0d0
	do 710 i=-1,1,2
	c0_ = c0_ - cspens(1d0+x(1,i)/alpha) + cspens(1d0+x(1,i)) 
     &		  - etass(1d0/alpha,-x(1,i))*log(1d0+x(1,i)/alpha) 
710	continue
	c0_ = c0_ + cspens(1d0+x(0,-1)/alpha) - cspens(1d0+x(0,-1)) 
     &		  + etass(1d0/alpha,-x(0,-1))*log(1d0+x(0,-1)/alpha) 
	c0_ = c0_ + log(alpha)*log((m22-ieps*m12-p20)/m12)
     &		  - log(alpha)**2/2d0
	c0_ = c0_/(p20-p12) 

	return
500	continue

C*** IR-divergent C0 function
        write(6,*)'warning: C0 is IR divergent'
        stop
	sc  = p12+abs(p12)*ieps
	mm1 = sqrt(m12)
	mm4 = sqrt(m22)
	xs = -4d0*mm1*mm4/(sc-(mm1-mm4)**2) /
     &	     ( sqrt(1d0-4d0*mm1*mm4/( sc-(mm1-mm4)**2))+1d0 )**2 
	c0_ = xs/mm1/mm4/(1d0-xs**2)*(
     &	    log(xs)*( -log(xs)/2d0+2d0*log(1d0-xs**2)
     &		    -log(lambda2/mm1/mm4) )
     &	  - pi**2/6d0+cspens(xs**2)+log(mm1/mm4)**2/2d0
     &	  + cspens(1d0-xs*mm1/mm4) + cspens(1d0-xs*mm4/mm1) )
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        FUNCTION CSPENS(Z)                                              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SPENCE-FUNKTION KOMPLEX, FREI NACH HOLLIK                     C
C---------------------------------------------------------------------C
C       20.07.83    LAST CHANGED 10.05.89        ANSGAR DENNER        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        COMPLEX*16 CSPENS,W,SUM,Z,U                                     
        REAL*8 RZ,AZ,A1                                                
        REAL*8 B(9)/                                                   
     1   0.1666666666666666666666666667D0,                             
     2  -0.0333333333333333333333333333D0,                             
     3   0.0238095238095238095238095238D0,                             
     4  -0.0333333333333333333333333333D0,                             
     5   0.0757575757575757575757575758D0,                             
     6  -0.2531135531135531135531135531D0,                             
     7   1.1666666666666666666666666667D0,                             
     8  -7.09215686274509804D0         ,                               
     9  54.97117794486215539D0         /                               
C     BEACHTE:                 B(N)=B2N                                
C     B(1)=1./6.                                                       
C     B(2)=-1./30.                                                     
C     B(3)=1./42.                                                      
C     B(4)=-1./30.                                                     
C     B(5)=5./66.                                                      
C     B(6)=-691./2730.                                                 
C     B(7)=7./6.                                                       
C     B(8)=-3617./510.                                                 
C     B(9)=43867./798.                                                 
C     B(10)=-174611./330.                                              
C     B(11)=854513./138.                                               
C     PI=3.1415926535897932384                                         
C     PI*PI/6.=1.6449..., PI*PI/3=3.28986...                           
C                                                                      
      Z =Z*DCMPLX(1D0)                                                 
      RZ=DREAL(Z)                                                      
      AZ=CDABS(Z)                                                      
      A1=CDABS(1D0-Z)                                                  
C     IF((SNGL(RZ) .EQ. 0.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN     
C ---> CHANGED  10.5.89                                                
      IF(AZ .LT. 1D-20) THEN                                           
        CSPENS=-CDLOG(1D0-Z)                                            
        RETURN                                                         
      END IF                                                           
c      IF((SNGL(RZ) .EQ. 1.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN     
c ---> changed 5.7.94
       IF( (ABS(RZ-1D0).LT.1D-18) .AND. (ABS(DIMAG(Z)).LT.1D-18) ) THEN     
        CSPENS=1.64493406684822643D0                                    
        RETURN                                                         
      END IF                                                           
      IF(RZ.GT.5D-1) GOTO 20                                           
      IF(AZ.GT.1D0) GOTO 10                                            
      W=-CDLOG(1D0-Z)                                                  
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 2                                     
      DO 1 K=1,9                                                       
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 2                            
      SUM=SUM+U*B(K)                                                   
 1    CONTINUE                                                         
 2    CSPENS=SUM                                                        
      RETURN                                                           
10    W=-CDLOG(1D0-1D0/Z)                                              
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 12                                    
                                                                       
      DO 11 K=1,9                                                      
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(B(K)*U/SUM).LT.1D-20) GOTO 12                           
      SUM=SUM+U*B(K)                                                   
11    CONTINUE                                                         
12    CSPENS=-SUM-1.64493406684822643D0-.5D0*CDLOG(-Z)**2               
      RETURN                                                           
20    IF(A1.GT.1D0) GOTO 30                                            
      W=-CDLOG(Z)                                                      
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 22                                    
      DO 21 K=1,9                                                      
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 22                           
      SUM=SUM+U*B(K)                                                   
21    CONTINUE                                                         
22    CSPENS=-SUM+1.64493406684822643D0-CDLOG(Z)*CDLOG(1D0-Z)           
      RETURN                                                           
30    W=CDLOG(1D0-1D0/Z)                                               
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 32                                    
      DO 31 K=1,9                                                      
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 32                           
      SUM=SUM+U*B(K)                                                   
31    CONTINUE                                                         
32    CSPENS=SUM+3.28986813369645287D0                                  
     *               +.5D0*CDLOG(Z-1D0)**2-CDLOG(Z)*CDLOG(1D0-Z)       
50    CONTINUE                                                         
        END                                                            
***********************************************************************
        FUNCTION ETASS(C1,C2)                                            
***********************************************************************
*       COMPLEX ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       8.06.90    ANSGAR DENNER                                       
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETASS,C1,C2                                           
        REAL*8     PI,IM1,IM2,IM12                                     
                                                                       
        PI     = 4D0*DATAN(1D0)                                        
        IM1    = DIMAG(C1)                                             
        IM2    = DIMAG(C2)                                             
        IM12   = DIMAG(C1*C2)                                          

	if (((IM1.eq.0d0).and.(DREAL(C1).lt.0d0)).or.
     &	    ((IM2.eq.0d0).and.(DREAL(C2).lt.0d0)).or.
     &	    ((IM12.eq.0d0).and.(DREAL(C1*C2).lt.0d0))) then
	  write(*,*) 'etass function on cut !!!'
	  write(*,*) 'C1    = ',C1
	  write(*,*) 'C2    = ',C2
	  write(*,*) 'C1*C2 = ',C1*C2
	  stop
	endif
                                                                       
        IF(IM1.LT.0D0.AND.IM2.LT.0D0.AND.IM12.GT.0D0) THEN             
            ETASS = DCMPLX(0D0,2D0*PI)                                   
        ELSE IF (IM1.GT.0D0.AND.IM2.GT.0D0.AND.IM12.LT.0D0) THEN       
            ETASS = DCMPLX(0D0,-2D0*PI)                                  
        ELSE                                                           
            ETASS = DCMPLX(0D0)                                          
        END IF                                                         
        END                                                            

***********************************************************************
        FUNCTION ETAS(Y,R,RS)                                            
***********************************************************************
*       MODIFIED ETA-FUNKTION                                           
*---------------------------------------------------------------------*
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETA,ETAS,Y,R,RS
        REAL*8     PI,IMY,IMRS
                                                                       
        PI     = 4D0*DATAN(1D0)                                        

	IF( DIMAG(R).NE.0D0 ) THEN
	    ETAS = ETA(Y,R)
	ELSE	    
	    IF( DREAL(R).GT.0D0 ) THEN
		ETAS = DCMPLX(0D0,0D0)
	    ELSE
	 	IMY  = DIMAG(Y)
		IMRS = DIMAG(RS)
		ETAS = 2D0*DCMPLX(0D0,PI)*(
     *			(1D0+SIGN(1D0,-IMY))*(1D0+SIGN(1D0,-IMRS))-
     *			(1D0+SIGN(1D0, IMY))*(1D0+SIGN(1D0, IMRS))
     *					  )/4D0
	    ENDIF
	ENDIF
        END                                                            

***********************************************************************
        FUNCTION SQE(A,B,C)                                            
***********************************************************************
*       SOLUTION OF QUADRATIC EQUATION				      *
*---------------------------------------------------------------------*
***********************************************************************
        IMPLICIT REAL*8 (A-Z)                                        
        COMPLEX*16 A,B,C,SQE,X1,X2

	X1=(-B+SQRT(B**2-4D0*A*C))/2D0/A
	X2=(-B-SQRT(B**2-4D0*A*C))/2D0/A

	IF (ABS(X1).GT.ABS(X2)) THEN
	   SQE=X1
	ELSE
	   SQE=X2
	ENDIF

        END                                                            

************************************************************************
        FUNCTION D0_(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4,ext) 
************************************************************************
*  SCALAR 4-POINT FUNCTION                                             *
*  P1,P2,P3,P4 = SQUARED EXTERNAL MOMENTA			       *
*  P12 = (p1+p2)**2,  P23 = (p2+p3)**2				       *
*----------------------------------------------------------------------*
*	General result from                                            *
*        A.Denner, U.Nierste and R.Scharf, Nucl. Phys. B367 (1991) 637 *
*	IR-divergent case from                                         *
*        W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349         *
*----------------------------------------------------------------------*
************************************************************************
        IMPLICIT REAL*8 (A-Z) 
	REAL*8 M(4),P(4,4),K(4,4)                                     
	COMPLEX*16 A1,A2,A3,A4,SWAP
	COMPLEX*16 SS(4), XX(2), X(2,4),RS(4,4)
	COMPLEX*16 S0(4),XX0(2),X0(2,4), R(4,4),G(2)
        COMPLEX*16 D0_,D0_ext,CSPENS,ETA,SQE,ETAS
	COMPLEX*16 AA,BB,CC,DD,IEPS,H,HH,L1,L2,L3,L4
	COMPLEX*16 SC,TC,XS,X2,X3,q2c,q3c,Y
	INTEGER I,J,ext

	COMMON/IR/LAMBDA,LAMBDA2


	if (ext.ne.0) then
c	  D0_ = D0_ext(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4,ext)
          write(6,*)'call to D0_ext not allowed'
          stop
	  return
	endif

        PI = 4D0*DATAN(1D0)                                               

        MM1=M1    
        MM2=M2    
        MM3=M3    
        MM4=M4    
        M12=M1*M1 
        M22=M2*M2 
        M32=M3*M3 
        M42=M4*M4 
        Q1=P1 
        Q2=P2   
        Q3=P3   
	Q4=P4
        Q12=P12   
        Q23=P23
C	IS AT LEAST ONE MASS ZERO ???
	IF (MM1*MM2*MM3*MM4.NE.0D0) GOTO 130

C--->	****** IR-divergent CASE ******
	IF ( ((Q1.EQ.M12).AND.(Q2.EQ.M32).AND.(M22.EQ.0D0)).OR.
     *	     ((Q2.EQ.M22).AND.(Q3.EQ.M42).AND.(M32.EQ.0D0)).OR.
     *	     ((Q3.EQ.M32).AND.(Q4.EQ.M12).AND.(M42.EQ.0D0)).OR.
     *	     ((Q4.EQ.M42).AND.(Q1.EQ.M22).AND.(M12.EQ.0D0)) ) goto 50

C--->	****** REGULAR CASE with at least one mass zero ******

C	PERMUTATE UNTIL MM3=0D0
	GOTO 20
10	CONTINUE
	MM0=MM1
	MM1=MM2
	MM2=MM3
	MM3=MM4
	MM4=MM0
	M02=M12
	M12=M22
	M22=M32
	M32=M42
	M42=M02
	Q00=Q12
	Q12=Q23
	Q23=Q00
	Q0=Q1
	Q1=Q2
	Q2=Q3
	Q3=Q4
	Q4=Q0
20	IF (MM3.NE.0D0) GOTO 10
C	ONLY MM3 IS ZERO
	IF (MM1*MM2*MM4.NE.0D0) GOTO 30
C	ONLY MM3 AND MM4 ARE ZERO ==> 3->2, 4->3...   
	IF ((MM1*MM2.NE.0D0).AND.(MM4.EQ.0D0)) GOTO 10
C	ONLY MM2 AND MM3 ARE ZERO
	IF ((MM1*MM4.NE.0D0).AND.(MM2.EQ.0D0)) GOTO 40
C	ONLY MM1 AND MM3 ARE ZERO ==> m1 <-> m2
	IF ((MM2*MM4.NE.0D0).AND.(MM1.EQ.0D0)) then
	  mm0 = mm1
	  mm1 = mm2
	  mm2 = mm0
	  q0  = q2
	  q2  = q12
	  q12 = q0
	  q0  = q4
	  q4  = q23
	  q23 = q0
	  goto 40
	endif
C 	check whether all masses are zero
	if ((mm1.eq.0d0).and.(mm2.eq.0d0).and.(mm4.eq.0d0)) then
	  WRITE(*,*) 'D0 case with all mi=0 not implemented !'
	  stop
	endif
C 	permutate until mm1 is non-zero
	if (mm1.eq.0d0) goto 10
	goto 70

C	****** NO MASS EQUAL TO ZERO ******
130	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)

	IF( ABS((MM1**2+MM3**2-Q12)/MM1/MM3).LT.2D0 ) THEN
C	R13 WOULD BE NOT REAL. -> PERMUTATION! -> R(2,4) IS NOT REAL.
	   M(1)=MM2
	   M(2)=MM3
	   M(3)=MM4
	   M(4)=MM1
	   P(1,2)=Q2
	   P(1,3)=Q23
	   P(1,4)=Q1
	   P(2,3)=Q3
	   P(2,4)=Q12
	   P(3,4)=Q4
	ELSE
C	R(1,3) IS REAL.
	   M(1)=MM1
	   M(2)=MM2
	   M(3)=MM3
	   M(4)=MM4
	   P(1,2)=Q1
	   P(1,3)=Q12
	   P(1,4)=Q4
	   P(2,3)=Q2
	   P(2,4)=Q23
	   P(3,4)=Q3
	ENDIF

	DO 11 J=2,4
	DO 11 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
11	CONTINUE

	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	S0(1)=R(1,2)
	S0(2)=R(2,3)
	S0(3)=R(3,4)
	S0(4)=R(1,4)
	AA=K(3,4)/R(2,4)+R(1,3)*K(1,2)-K(1,4)*R(1,3)/R(2,4)-K(2,3)
	BB=(R(2,4)-1D0/R(2,4))*(R(1,3)-1D0/R(1,3))
     *		+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)/R(1,3)+R(2,4)*K(3,4)-K(1,4)*R(2,4)/R(1,3)-K(2,3)
	DD=K(2,3)-R(1,3)*K(1,2)-R(2,4)*K(3,4)+R(1,3)*R(2,4)*K(1,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	XX0(1)=SQE(AA,BB,CC)
	XX0(2)=CC/AA/XX0(1)

c	XX(1)=XX0(1)-IEPS*DD/(XX0(1)-XX0(2))
c	XX(2)=XX0(2)+IEPS*DD/(XX0(1)-XX0(2))

c	IF (ABS(DREAL(XX0(1)-XX(2))).LT.ABS(DREAL(XX0(1)-XX(1)))) THEN
	IF (ABS(XX0(1)-XX(2)).LT.ABS(XX0(1)-XX(1))) THEN
	  SWAP  =XX0(1)
	  XX0(1)=XX0(2)
	  XX0(2)=SWAP
	ENDIF

	DO 12 I=1,2
	G(I)  =SIGN( 1D0,DREAL(AA*(XX(I)-XX(3-I))) )
	 X(I,1)= XX(I)/R(2,4)
	X0(I,1)=XX0(I)/R(2,4)
	 X(I,2)= XX(I)/R(2,4)*R(1,3)
	X0(I,2)=XX0(I)/R(2,4)*R(1,3)
	 X(I,3)= XX(I)*R(1,3)
	X0(I,3)=XX0(I)*R(1,3)
	 X(I,4)= XX(I)
	X0(I,4)=XX0(I)
12	CONTINUE

	D0_ = DCMPLX(0D0,0D0)
	DO 13 I=1,2
	DO 13 J=1,4
	A1 = 1D0+X0(I,J)*S0(J) + ABS(1D0+X0(I,J)*S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)*SS(J)))
	A2 = 1D0+X0(I,J)/S0(J) + ABS(1D0+X0(I,J)/S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)/SS(J)))
	D0_ = D0_ + (-1D0)**(I+J)*(
     *		CSPENS(A1)+ETA(-X(I,J),SS(J))*LOG(A1)
     *	       +CSPENS(A2)+ETA(-X(I,J),1D0/SS(J))*LOG(A2)     )
13	CONTINUE

	IF( DIMAG(R(1,3)).EQ.0D0 ) THEN
	DO 14 I=1,2
	   A1 = (K(1,3)-2D0*R(1,3))/XX0(I)
     *		      -R(1,3)*K(1,4)+K(3,4)
     	   A2 = ((K(2,4)-2D0*R(2,4))*R(1,3)*XX0(I)
     *		      -R(2,4)*K(3,4)+K(2,3))/DD
	   A3 = (K(1,3)-2D0*R(1,3))*R(2,4)/XX0(I)
     *		      -R(1,3)*K(1,2)+K(2,3)
	   A4 = ((K(2,4)-2D0*R(2,4))*XX0(I)
     *		      -R(2,4)*K(1,4)+K(1,2))/DD
	   L1 = LOG( A1-ABS(A1)*IEPS )
     	   L2 = LOG( A2+ABS(A2)*IEPS*G(I)*SIGN(1D0,DREAL(R(1,3))
     *				        	  *DIMAG(RS(2,4))) ) 
	   L3 = LOG( A3-ABS(A3)*IEPS )
	   L4 = LOG( A4+ABS(A4)*IEPS*G(I)*SIGN(1D0,DIMAG(RS(2,4))) ) 

	   D0_ = D0_ + (3D0-2D0*I)*(
     *		 ETAS(-XX(I),R(1,3),RS(1,3))
     *		   *( LOG(R(1,3)*XX(I)) + L1 + L2 )
     *		+ETAS(-XX(I),1D0/R(2,4),1D0/RS(2,4))
     *		   *( LOG(XX(I)/R(2,4)) + L3 + L4 )
     *		-( ETAS(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))
     *		  +ETA(RS(1,3),1D0/RS(2,4)) )
     *		   *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 + L2 )
     *	  	+ETA(RS(1,3),1D0/RS(2,4))
     *		   *ETAS(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))   )
14	CONTINUE
	ELSE
	DO 15 I=1,2
	   L1 = LOG( R(2,4)/XX0(I)+XX0(I)/R(2,4)+K(1,2)
     *		     -XX0(I)/R(2,4)*EPS*BB*G(I) )
	   L2 = LOG( R(1,3)*XX0(I)+1D0/XX0(I)/R(1,3)+K(3,4)
     *		     -XX0(I)*R(1,3)*EPS*BB*G(I) )
	   L3 = LOG( R(1,3)/R(2,4)*XX0(I)+R(2,4)/XX0(I)/R(1,3)+K(2,3)
     *		     -XX0(I)*R(1,3)/R(2,4)*EPS*BB*G(I) )
	   D0_ = D0_ + (3D0-2D0*I)*(
     *		+ETA(-XX(I),1D0/R(2,4))
     *		   *( LOG(XX(I)/R(2,4)) + L1 )
     *		+ETA(-XX(I),R(1,3))
     *		   *( LOG(R(1,3)*XX(I)) + L2 )
     *		-( ETA(-XX(I),R(1,3)/R(2,4))
     *		  +ETA(R(1,3),1D0/R(2,4)) )
     *		   *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 )
     *	  	+ETA(R(1,3),1D0/R(2,4))
     *		   *ETA(-XX(I),-R(1,3)/R(2,4))
     *		   *(1D0-G(I)*SIGN(1D0,DREAL(BB)))	    )
15	CONTINUE
	ENDIF

	D0_ = D0_/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN

C	****** ONLY MM3 IS ZERO ******
30	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)
	M(1)=MM1
	M(2)=MM2
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 1 J=2,4
	DO 1 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
1	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(3,4)/R(2,4)-K(2,3)
	BB=K(1,3)*(1D0/R(2,4)-R(2,4))+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(1,3)*K(1,4)*R(2,4)+R(2,4)*K(3,4)-K(2,3)
	DD=K(2,3)-R(2,4)*K(3,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 2 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
2	CONTINUE
	D0_ = DCMPLX(0D0,0D0)
	DO 3 I=1,2
	D0_ = D0_ + (2D0*I-3D0)*(
     *		CSPENS(1D0+SS(4)*X(I,4))
     *	       -CSPENS(1D0+SS(1)*X(I,1))
     *	       +CSPENS(1D0+X(I,4)/SS(4))
     *	       -CSPENS(1D0+X(I,1)/SS(1))
     *	       +ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       -ETA(-X(I,1),SS(1))*LOG(1D0+SS(1)*X(I,1))
     *	       +ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -ETA(-X(I,1),1D0/SS(1))*LOG(1D0+X(I,1)/SS(1))
     *	       -CSPENS(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +CSPENS(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA(-X(I,4),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +ETA(-X(I,1),(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))   )
	IF (DIMAG(R(2,4)).NE.0D0) THEN
	   H=ETA(-1D0/XX(I),R(2,4))
	ELSE
	   H=DCMPLX(0D0,0D0)
	   IF (DREAL(R(2,4)).LT.0D0) THEN
	      HH=-1D0/XX(I)
	      IM1=DIMAG(HH)
	      IM2=DIMAG(RS(2,4))
	      IF ((IM1.GT.0D0).AND.(IM2.GT.0D0)) THEN
	         H=-DCMPLX(0D0,2D0*PI)
	      ENDIF
	      IF ((IM1.LT.0D0).AND.(IM2.LT.0D0)) THEN
	         H=+DCMPLX(0D0,2D0*PI)
	      ENDIF
	   ENDIF
	ENDIF
	D0_ = D0_ + (2D0*I-3D0)*
     *	          H*( LOG( (K(1,2)-R(2,4)*K(1,4)
     *			  +XX(I)*(1D0/R(2,4)-R(2,4)))/DD )
     *		     +LOG(K(1,3)-IEPS) )
3	CONTINUE
	D0_ = D0_/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN

C	****** ONLY MM2 AND MM3 ARE ZERO ******
40	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)

	M(1)=MM1
	M(2)=10D0
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 4 J=2,4
	DO 4 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.2) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.2) K(I,J)=K(I,J)-M(J)/M(I)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
4	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(2,4)*K(3,4)-K(2,3)
	BB=K(1,3)*K(2,4)+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(2,3)
	DD=K(2,3)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 5 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
5	CONTINUE
	D0_ = DCMPLX(0D0,0D0)
	DO 6 I=1,2
	D0_ = D0_ + (2D0*I-3D0)*(
     *		CSPENS(1D0+SS(4)*X(I,4))
     *	       +CSPENS(1D0+X(I,4)/SS(4))
     *	       +ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       +ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -CSPENS(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -CSPENS(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       -ETA(-XX(I),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA(-XX(I),(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS)) 
     *	       +LOG(-XX(I))*( LOG(K(1,2)-IEPS)
     *			     +LOG(K(1,3)-IEPS)-LOG(K(2,3)-IEPS) ) )
6	CONTINUE
	D0_ = D0_/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))

	return

C	****** ONLY MM1 IS NON-ZERO ******
70	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)

	M(1)=MM1
	M(2)=10D0
	M(3)=10D0
	M(4)=10D0
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	k(1,2) = (m(1)**2-p(1,2))/m(1)/m(2)
	k(1,3) = (m(1)**2-p(1,3))/m(1)/m(3)
	k(1,4) = (m(1)**2-p(1,4))/m(1)/m(4)
	k(2,3) = -p(2,3)/m(2)/m(3)
	k(2,4) = -p(2,4)/m(2)/m(4)
	k(3,4) = -p(3,4)/m(3)/m(4)
	DO 74 J=2,4
	DO 74 I=1,J-1
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
74	CONTINUE
	AA=K(2,4)*K(3,4)
	BB=K(1,3)*K(2,4)+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(2,3)
	DD=K(2,3)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	D0_ = DCMPLX(0D0,0D0)
	DO 76 I=1,2
	D0_ = D0_ + (2D0*I-3D0)*(
     *	        CSPENS(1D0+(K(1,4)-IEPS)*xx(i))
     *	       +eta(-xx(i),K(1,4)-IEPS)*log(1D0+(K(1,4)-IEPS)*xx(i))
     *	       -CSPENS(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -CSPENS(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       -ETA(-XX(I),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA(-XX(I),(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS)) 
     *	       +LOG(-XX(I))*( LOG(K(1,2)-IEPS)
     *			     +LOG(K(1,3)-IEPS)-LOG(K(2,3)-IEPS) ) )
76	CONTINUE
	D0_ = D0_/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	return

C	****** general IR-divergent D0-function ******
50	CONTINUE
        write(6,*)'warning: D0 is IR divergent'
        stop
C	PERMUTATE UNTIL MM1 IS THE PHOTON
	GOTO 52
51	CONTINUE
	MM0=MM1
	MM1=MM2
	MM2=MM3
	MM3=MM4
	MM4=MM0
	M02=M12
	M12=M22
	M22=M32
	M32=M42
	M42=M02
	Q00=Q12
	Q12=Q23
	Q23=Q00
	Q0=Q1
	Q1=Q2
	Q2=Q3
	Q3=Q4
	Q4=Q0
52	IF ((mm1.ne.0d0).or.(q1.ne.m22).or.(q4.ne.m42)) GOTO 51

	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)
	sc  = q23+abs(q23)*ieps
	tc  = q12+abs(q12)*ieps
	q2c = q2+abs(q2)*ieps
	q3c = q3+abs(q3)*ieps
	xs = -4d0*mm2*mm4/(sc-(mm2-mm4)**2) /
     &	     ( sqrt(1d0-4d0*mm2*mm4/( sc-(mm2-mm4)**2))+1d0 )**2 

	if (mm3.eq.0d0) goto 60

C	*** general case ***
	if (q2.ne.(mm2-mm3)**2) then
	  x2 = -4d0*mm2*mm3/(q2c-(mm2-mm3)**2) /
     &	     ( sqrt(1d0-4d0*mm2*mm3/( q2c-(mm2-mm3)**2))+1d0 )**2 
	else
	  x2 = 1d0
	endif
	if (q3.ne.(mm4-mm3)**2) then
	  x3 = -4d0*mm4*mm3/(q3c-(mm4-mm3)**2) /
     &	     ( sqrt(1d0-4d0*mm4*mm3/( q3c-(mm4-mm3)**2))+1d0 )**2 
	else
	  x3 = 1d0
	endif

	d0_ = xs/mm2/mm4/(q12-m32)/(1d0-xs**2)*(
     &	 2d0*cdlog(xs)*(cdlog(1d0-xs**2)-cdlog(lambda*mm3/(m32-tc)))
     &	+pi**2/2d0+cspens(xs**2)+cdlog(x2)**2+cdlog(x3)**2
     &	-cspens(xs*x2*x3)-(cdlog(xs)+cdlog(x2)+cdlog(x3))
     $       *cdlog(1d0-xs*x2*x3)
     &	-cspens(xs*x2/x3)-(cdlog(xs)+cdlog(x2)-cdlog(x3))
     $       *cdlog(1d0-xs*x2/x3)
     &	-cspens(xs/x2*x3)-(cdlog(xs)-cdlog(x2)+cdlog(x3))
     $       *cdlog(1d0-xs/x2*x3)
     &	-cspens(xs/x2/x3)-(cdlog(xs)-cdlog(x2)-cdlog(x3))
     $       *cdlog(1d0-xs/x2/x3) )
	return

60	continue
C	*** special case: mass mm3 opposite to photon is 0 ***
	if ((q2.eq.m22).or.(q3.eq.m42)) goto 61
	Y = mm2/mm4*(q3c-m42)/(q2c-m22)
	d0_ = xs/mm2/mm4/q12/(1d0-xs**2)*(
     &	 log(xs)*( -log(xs)/2d0+2d0*log(1d0-xs**2)-log(lambda2/mm2/mm4)
     &		   -log((q2-m22)/tc)-log((q3-m42)/tc) )
     &	+pi**2/6d0+cspens(xs**2)+log(y)**2/2d0
     &	-cspens(xs*y)-(log(xs)+log(y))*log(1d0-xs*y)
     &	-cspens(xs/y)-(log(xs)-log(y))*log(1d0-xs/y) )
	return

61	continue
C	*** special case: doubly IR-divergent D0 ***
	if ((q2.eq.m22).and.(q3.eq.m42)) then
	d0_ = -xs/mm2/mm4/q12/(1d0-xs**2)*2d0*log(xs)*log(-lambda2/tc)
	else
	  write(*,*) 'Special case of IR-divergent D0 not implemented!'
          stop
	endif

	END
***********************************************************************
      complex*16 function fs(q2,rm,rn)
***********************************************************************
*       skalares einschleifenintegral, doppeltlang, 'regulaerer anteil'
*       f(q2,rm,rn)=b0(q2,rm,rn)-b0(0d0,rm,rn)  'subtrahiertes f'
*       q2=quadrat des die schleife durchlaufenden impulses
*       rm,rn: massen der teilchen auf beiden armen
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       19.10.83
************************************************************************
      real*8 m,n,pi,a,s,t,b,c,d,q2,rm,rn,u,v,w,m2,n2
      data pi/3.1415926535897932384626438d0/
      m=rm
      n=rn
      m2 = m**2
      n2 = n**2
      if (dabs(m) .eq. dabs(n) ) goto 30
      if (n2 .eq. 0d0) goto 310
      if (m2 .eq. 0d0) goto 300
c---- >  allgemeiner fall
      if (q2 .ne. 0d0) goto 520
      b=0d0
      a=0d0
      goto 560
 520  u=m*m+n*n
      v=m*m-n*n
      w=m*n
      if (dabs(q2/v).le.1d-4) then
         b = u/v/v/2d0+2d0*w*w/v/v/v*dlog(n2/m2)/2d0
         b = b*q2
         a = 0d0
         goto 570
      end if
      s=dabs(m)+dabs(n)
      t=dabs(m)-dabs(n)
      c=dsqrt(dabs(s*s-q2))
      d=dsqrt(dabs(t*t-q2))
      b=1d0+(v/q2-u/v)*dlog(n2/m2)/2d0
      if (2d0*w .le. dabs(q2-u)) goto 550
      b=b-2d0*c*d/q2*datan(d/c)
      a=0d0
      goto 560
 550  a=c*d/q2
      b=b-dsign(1d0,q2-u)*a*dlog((c+d)*(c+d)/(4d0*w))
      a=pi*a
      if (q2 .ge. u) goto 560
      a=0d0
 560  continue
 570  fs=dcmplx(b,a)
      return
c---- >  gleiche massen
 30   if (q2 .ne. 0d0) goto 40
      b=0d0
      a=0d0
      goto 560
 40   u=4d0*m*m
      v=dsqrt(dabs(1d0-u/q2))
      if ((q2 .ge. 0d0) .and. (q2 .lt. u)) goto 50
      b=2d0-v*dlog((v+1d0)*(v+1d0)/u*dabs(q2))
      a=pi*v
      if (q2 .ge. u) goto 560
      a=0d0
      goto 560
 50   b=2d0-2d0*v*datan(1d0/v)
      a=0d0
      goto 560
c---- >  eine masse null
 300  m=n
 310  if (q2 .ne. m*m) goto 320
      a=0d0
      b=1d0
      goto 560
 320  b=1d0
      if(q2 .eq. 0d0) b=0d0
      a=b*(1d0-m*m/(q2+(1d0-b)))
      b=b-a*dlog(dabs(1d0-q2/m/m))
      a=pi*a
      if (q2 .gt. m*m) goto 560
      a=0d0
      goto 560
      end
***********************************************************************
      complex*16 function b0(q2,rm,rn)
***********************************************************************
*       skalares einschleifenintegral b0 minus
*       divergenter anteil (only (delta-log mue**2) subtracted)
*       rm,rn: massen der teilchen auf beiden armen
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       04.02.87 sa
************************************************************************
      implicit real*8(a,c,g,q,m,s,r)
      complex*16 fs
*      
      b0 = dcmplx(0d0)
      r12 = rm*rm
      r22 = rn*rn
      if (dabs(rm).eq.dabs(rn)) goto 110
      if ((r22.eq.0d0).or.(r12.eq.0d0)) goto 120
      b0 = dcmplx(-dlog(r12)/2d0-dlog(r22)/2d0+1d0-(r12+r22)/(r12-r22)
     1     *dlog(r12/r22)/2d0,0d0)
      goto 200
c---  >   beide massen gleich
 110  if (r12.eq.0d0) goto 300
      b0 = dcmplx(-dlog(r12))
      goto 200
c---  >   eine masse null
 120  b0 = 1d0-dcmplx(dlog(r12+r22))
 200  b0 = b0+fs(q2,rm,rn)
      return
c---  >   beide massen gleich null
 300  b0 = 2d0-cdlog(dcmplx(q2,1d-20))
      return
      end
***********************************************************************
      complex*16 function b1(q2,rm,rn)
***********************************************************************
*     skalares einschleifenintegral b1 minus
*     divergenter anteil (only -1/2(delta-log(mue**2)) subtracted)
*     rm,rn: massen der teilchen auf beiden armen
*     d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       04.02.87 sa
************************************************************************
      implicit real*8(a,c,g,q,m,s,r)
      complex*16 fs,b0,xb0,xb1
*
      r12= rm*rm
      r22= rn*rn
      b1 = dcmplx(0d0)
      if (dabs(rm).eq.dabs(rn)) goto 200
      if (r12.eq.0d0) goto 120
      if (r22.eq.0d0) goto 130
      if (q2.eq.0d0) goto 140
      b1 = (r22-r12)*b0(q2,rm,rn)+r12*(1d0-dlog(r12))
     1     -r22*(1d0-dlog(r22))
      b1 = b1/2d0/q2
      goto 200
 120  b1 = r22*fs(q2,rm,rn)/2d0/q2
      goto 200
 130  b1 = -r12*fs(q2,rm,rn)/2d0/q2
 200  b1 = b1-b0(q2,rm,rn)/2d0
      goto 300
c Achtung: original code has been changed here!
 140  call bquer(q2,rm,rn,xb0,xb1)
      b1 = xb1+dlog(r22)/2d0
 300  continue
      return
      end
***********************************************************************
      complex*16 function b20(q2,rm,rn)
***********************************************************************
*       skalares einschleifenintegral b0 minus
*       divergenter anteil
*       rm,rn: massen der teilchen auf beiden armen
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       04.02.87 sa
************************************************************************
      implicit real*8(r,m,q)
      complex*16 b0,b1,b21
*      
      r12= rm*rm
      r22= rn*rn
      b20 = (q2+r12-r22)*b1(q2,rm,rn)+2d0*r12*b0(q2,rm,rn)
     1     +r12+2d0*r22-q2/3d0
      if (r22.ne.0d0) b20 = b20-r22*dlog(r22)
      b20 = b20/6d0
      return
      entry b21(q2,rm,rn)
      r12= rm*rm
      if ((dabs(rm).eq.dabs(rn)).and.(q2.eq.0d0)) goto 100
      r22= rn*rn
      b21 = 2d0*(r12-r22+q2)*b1(q2,rm,rn)+r12*b0(q2,rm,rn)
     1     +(r12-r22-q2/3d0)/2d0
      if (r22.ne.0d0) b21 = b21+r22*dlog(r22)
      b21 = -b21/3d0/q2
      return
 100  b21 = -dlog(r12)/3d0
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function c0(cq1,cq2,cq3,mm1,mm2,mm3)
***********************************************************************
*       skalare dreipunktfunktion, doppeltlang,
*       q1,q2,q3=quadrate der impulse
*       m1,m2,m3: massenquadrate der teilchen auf inneren linien
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       20.01.87 sa
************************************************************************
      implicit real*8(m,q,c)
      implicit complex*16(z)
      complex*16 y(3,3),spence,cywur
      real*8 cj(3)
*
      j = 3
      cj(1) = -1d0
      cj(2) = 1d0
      cj(3) = -1d0
      m1 = mm1*mm1
      m2 = mm2*mm2
      m3 = mm3*mm3
      q1 = cq1
      q2 = cq2
      q3 = cq3
      if (q2 .eq. 0d0) then
         q2 = cq1
         q1 = cq2
         m3 = mm1*mm1
         m1 = mm3*mm3
      end if
      if (q1 .eq. 0d0) then
         j = 2
      end if
      c1 = 1d0+(q1-q3)/q2
      if (c1*c1-4d0*q1/q2 .lt. 0d0) then
         write(6,*) ' neg. determinante !'
         stop
      end if
      cwurz = dsqrt(c1*c1-4d0*q1/q2)
      calphs = (c1-cwurz)/2d0
      c0 = 0d0
      ca = q1
      cb = q2
      cc = q3-q2-q1
      cd = m2-m1-q1
      ce = q1-q3+m3-m2
      cf = m1
      cnum = -(cd+ce*calphs)
      cdenom = cc+2d0*calphs*cb
      y(1,3) = dcmplx((cnum-2d0*ca-cc*calphs)/cdenom,0d0)
      cy0 = -cc-ce
      cywur = cdsqrt(dcmplx(cy0*cy0-4d0*cb*(ca+cd+cf),1d-20))
      y(1,1) = (cy0+cywur)/2d0/cb
      y(1,2) = (cy0-cywur)/2d0/cb
      cy0 = -ce-cd
      cywur = cdsqrt(dcmplx(cy0*cy0-4d0*cf*(ca+cb+cc),1d-20))
      y(2,3) = dcmplx(cnum/cdenom/(1d0-calphs),0d0)
      y(2,1) = (cy0+cywur)/2d0/(ca+cb+cc)
      y(2,2) = (cy0-cywur)/2d0/(ca+cb+cc)
      cywur = cdsqrt(dcmplx(cd*cd-4d0*ca*cf,1d-20))
      if (j .eq. 3) then
         y(3,3) = dcmplx(-cnum/calphs/cdenom,0d0)
         y(3,1) = (-cd+cywur)/2d0/ca
         y(3,2) = (-cd-cywur)/2d0/ca
      end if
      do 100 i=1,j
         c0 = c0+cj(i)*dreal(spence(y(i,3)/(y(i,3)-y(i,1)))
     1        -spence((y(i,3)-1d0)/(y(i,3)-y(i,1)))
     2        +spence(y(i,3)/(y(i,3)-y(i,2)))
     3        -spence((y(i,3)-1d0)/(y(i,3)-y(i,2))))
 100  continue
      c0 = c0/cdenom
      return
      end
***********************************************************************
      subroutine cmue(q1,q2,q3,m1,m2,m3,c1)
***********************************************************************
*     nicht-skalare dreipunktfunktion, doppeltlang,
*     q1,q2,q3=quadrate der impulse
*     m1,m2,m3: massen der teilchen auf inneren linien
*     nicht ir divergent
*     d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       20.01.87 sa
************************************************************************
      implicit real*8(a,c,g,q,m,s,d)
      implicit complex*16(z)
      complex*16 b0,spence,cscal,c0_
      real*8 c1(0:2)
*
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
c---> two gluons on shell, 3 equal internal fermions
      if ((q1.eq.0d0).and.(q2.eq.0d0).and.(q3.gt.0d0)) goto 230
      goto 240
 230  if (4d0*m12.ge.q3) then
         awur = dsqrt(4d0*m12/q3-1d0)
         z1 = dcmplx(1d0,awur)/2d0
         z2 = dcmplx(1d0,-awur)/2d0
      else
         awur = dsqrt(1d0-4d0*m12/q3)
         al1 = 1d0+awur
         al2 = 1d0-awur
         z1 = dcmplx(al1,al2*1d-10)/2d0
         z2 = dcmplx(al2,-al2*1d-10)/2d0
      end if
      c1(0) = -spence(1d0/z1)-spence(1d0/z2)
      c1(0) = c1(0)/q3
c     wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
c     write(6,*)'Check this point in ''boxlib.f'' at line 408'
c     write(6,*)'                                       ---'
c     wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
c cscal should be used for unequal internal masses !
c      c1(0) = dreal(cscal(q3,0d0,m1,m3,m2))
      c1(0)=dreal(C0_(q1,q2,q3,m1,m2,m3,1))
      goto 100
 240  qq2 = q2
      qq1 = q1
      qq3 = q3
      mm1 = m1
      mm2 = m2
      mm3 = m3
      if ((dabs(q1).gt.0d0).and.(q3.eq.0d0)) then
         qq3 = q1
         qq1 = q3
         mm2 = m3
         mm3 = m2
      end if
c 300  c1(0) = c0(qq1,qq2,qq3,mm1,mm2,mm3)
 300  c1(0)=dreal(C0_(qq1,qq2,qq3,mm1,mm2,mm3,1))
      goto 100
 100  matr = (q3-q2-q1)/2d0
      det = q1*q2-matr*matr
      cb3 = dreal(b0(q3,m1,m3))
      cb2 = dreal(b0(q2,m2,m3))
      cb1 = dreal(b0(q1,m1,m2))
      cmp1 = (cb3-cb2-(m12-m22+q1)*c1(0))/2d0
      cmp2 = (cb1-cb3-(m22-m32+q3-q1)*c1(0))/2d0
      c1(1) = (q2*cmp1-matr*cmp2)/det
      c1(2) = (q1*cmp2-matr*cmp1)/det
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cmn(q1,q2,q3,m1,m2,m3,c1,c2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     dreipunktfkt. mit 2 integrationsimpulsen in zaehler
c     realteil
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     04.02.87 sa
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(c,q,m)
      complex*16 b0,b1
      real*8 det,c1(0:2),c2(0:3)
*      
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      matr = (q3-q2-q1)/2d0
      det = q1*q2-matr*matr
      cf1 = m12-m22+q1
      cf2 = m22-m32+q3-q1
      cb0 = dreal(b0(q2,m2,m3))
      call cmue(q1,q2,q3,m1,m2,m3,c1)
      c2(0) = m12/2d0*c1(0)+(cb0+cf1*c1(1)+cf2*c1(2)+1d0)/4d0
      cb1 = dreal(b1(q3,m1,m3))
      cb2 = dreal(b1(q1,m1,m2))
      cmp1 = cb1+cb0-cf1*c1(1)-2d0*c2(0)
      cmp2 = cb2-cb1-cf2*c1(1)
      c2(1) = q2*cmp1-matr*cmp2
      c2(1) = c2(1)/2d0/det
      c2(3) =-matr*cmp1+q1*cmp2
      c2(3) = c2(3)/2d0/det
      cb3 = b1(q2,m2,m3)
      cmp1 = cb1-cb3-cf1*c1(2)
      cmp2 = -cb1-cf2*c1(2)-2d0*c2(0)
      c2(2) = -matr*cmp1+q1*cmp2
      c2(2) = c2(2)/2d0/det
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chutmn(q1,q2,q3,m1,m2,m3,ch1,ch2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     dreipunktfkt. mit 2 integrationsimpulsen in zaehler
c     realteil
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     06.05.87 sa
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(c,q,m)
      real*8 c1(0:2),c2(0:3),ch1(0:2),ch2(0:3)
*
      call cmn(q1,q2,q3,m1,m2,m3,c1,c2)
      ch1(0) = c1(0)
      ch1(1) = c1(1)-c1(2)
      ch1(2) = c1(2)
      ch2(0) = c2(0)
      ch2(1) = c2(1)-2d0*c2(3)+c2(2)
      ch2(2) = c2(2)
      ch2(3) = c2(3)-c2(2)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 function d0gen(q1,q2,q3,q4,q5,q6,m1,m2,m3,m4)
***********************************************************************
*       skalare vierpunktfunktion, doppeltlang,
*       za-zf: parameters a-f from 't hooft veltman
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       07.08.90 sa
************************************************************************
      implicit real*8(m,q,p)
      implicit complex*16(c,z)
      complex*16 cm2(4),cmm2(4),cl(4,4),cq(4,4),ca(4),ca2(4)
      real*8 p(4,4)
*
      cm2(1) = dcmplx(m1*m1,-1d-20)
      cm2(2) = dcmplx(m2*m2,-1d-20)
      cm2(3) = dcmplx(m3*m3,-1d-20)
      cm2(4) = dcmplx(m4*m4,-1d-20)
      p(1,2) = q1
      p(1,3) = q4
      p(1,4) = q6
      p(2,3) = q2
      p(2,4) = q1+q2+q3+q6-q4-q5
      p(3,4) = q3
      do 110 i1=1,4
         do 100 i2=i1+1,4
            cl(i1,i2) = p(i1,i2)-cm2(i1)-cm2(i2)
 100     continue
 110  continue
      ca(2) = dcmplx(1d0/m1/m1)
      ca2(2) = ca(2)*ca(2)
      ca(1) = cl(1,2)-cdsqrt(cl(1,2)*cl(1,2)-4d0*cm2(1)*cm2(2))
      ca(1) = -ca(1)/2d0/cm2(1)*ca(2)
      ca2(1) = ca(1)*ca(1)
      ca(3) = (-cm2(2)*ca2(2)+cm2(1)*ca2(1))
     1     /(cl(2,3)*ca(2)-cl(1,3)*ca(1))
      ca2(3) = ca(3)*ca(3)
      ca(4) = (-cm2(2)*ca2(2)+cm2(1)*ca2(1))
     1     /(cl(2,4)*ca(2)-cl(1,4)*ca(1))
      ca2(4) = ca(4)*ca(4)
      do 210 i1=1,4
         cmm2(i1) = -cm2(i1)*ca2(i1)
         do 200 i2=i1+1,4
            cq(i1,i2) = cl(i1,i2)*ca(i1)*ca(i2)+cm2(i1)*ca2(i1)
     1           +cm2(i2)*ca2(i2)
 200     continue
 210  continue
      caa = -cq(3,4)
      cb = -cq(2,3)
      cg = -cq(1,2)
      cc = -cq(2,4)+cq(2,3)+cq(3,4)
      ch = -cq(1,4)-cq(2,3)+cq(1,3)+cq(2,4)
      cj = -cq(1,3)+cq(1,2)+cq(2,3)
      cd = cmm2(3)-cmm2(4)+cq(3,4)
      ce = cmm2(2)-cmm2(3)+cq(2,4)-cq(3,4)
      ck = cmm2(1)-cmm2(2)+cq(1,4)-cq(2,4)
      cf = cmm2(4)
      d0gen = c0gen(caa,cb,cc,cd,ce,cf)
     1     -c0gen(caa,cb,cc,cd,(ce+ck),cf)
      d0gen = -ca(1)*ca(2)*ca(3)*ca(4)*d0gen/ck
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 function c0gen(za,zb,zc,zd,ze,zf)
***********************************************************************
*       skalare dreipunktfunktion, doppeltlang,
*       za-zf: parameters a-f from 't hooft veltman
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       07.08.90 sa
************************************************************************
      implicit real*8(m,q)
      implicit complex*16 (c,z)
      complex*16 y(3,3),c1(3),c2(3),yd1,yd2
      complex*16 spence
      real*8 cj(3)
      j = 3
      cj(1) = -1d0
      cj(2) = 1d0
      cj(3) = -1d0
*      
      cwurz = cdsqrt(zc*zc-4d0*za*zb)
      calphs = (-zc-cwurz)/2d0/zb
      c0 = dcmplx(0d0)
      ca = za
      cb = zb
      cc = zc
      cd = zd
      ce = ze
      cf = zf
      cnum = -(cd+ce*calphs)
      cdenom = cc+2d0*calphs*cb
      y(1,3) = (cnum-2d0*ca-cc*calphs)/cdenom
      cy0 = -cc-ce
      cywur = cdsqrt(cy0*cy0-4d0*cb*(ca+cd+cf))
      y(1,1) = (cy0+cywur)/2d0/cb
      y(1,2) = (cy0-cywur)/2d0/cb
      c1(1) = cb
      c2(1) = ca+cd+cf
      cy0 = -ce-cd
      cywur = cdsqrt(cy0*cy0-4d0*cf*(ca+cb+cc))
      y(2,3) = cnum/cdenom/(1d0-calphs)
      y(2,1) = (cy0+cywur)/2d0/(ca+cb+cc)
      y(2,2) = (cy0-cywur)/2d0/(ca+cb+cc)
      c1(2) = ca+cb+cc
      c2(1) = cf
      cywur = cdsqrt(cd*cd-4d0*ca*cf)
      y(3,3) = -cnum/calphs/cdenom
      y(3,1) = (-cd+cywur)/2d0/ca
      y(3,2) = (-cd-cywur)/2d0/ca
      c1(3) = ca
      c2(3) = cf
      do 100 i=1,3
         yd1 = y(i,3)-y(i,1)
         yd2 = y(i,3)-y(i,2)
         c0 = c0+cj(i)*(spence(y(i,3)/yd1)
     1        -spence((y(i,3)-1d0)/yd1)
     2        +spence(y(i,3)/yd2)
     3        -spence((y(i,3)-1d0)/yd2)
     5        )
 100  continue
      c0gen = c0/cdenom
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function ggttb(s,t,mex2,mib2,mif2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       mass singular scalar box integral gg-->tt
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       21.06.90 sa
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(a-b,d-y)
      implicit complex*16(c,z)
      complex*16 spence
*
      pi = 4d0*datan(1d0)
      cs = dcmplx(s,s*1d-10)
      cmh2 = dcmplx(mib2,-mib2*1d-10)
      tmt = t-mex2
      cmht = cmh2-t
      cmhmt = cmh2-mex2
      mht = mib2-t
      mhmt = mib2-mex2
      c1 = -tmt/cmht
      c2 = -s*cmh2/cmhmt/cmhmt
      zlogi = dcmplx(dlog(mif2/s),pi)
      zlogh = cdlog(cmhmt)
      zlogs = dcmplx(dlog(mib2/s),pi)
      log1 = dlog(mht)
      log2 = dlog(mib2)
      log3 = dlog((s*mib2+mhmt*mhmt)/mhmt/mhmt)
      cd0 = 4d0*spence(c1)+spence(c2)
     1     -zlogi*(2d0*zlogh-2d0*log1+zlogi/2d0)
     2     +log3*(2d0*log2-2d0*zlogh-zlogs)
      ggttb = dreal(cd0)/s/mht
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function ggttf(s,t,mex2,mib2,mif2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       scalar box integral gg-->tt
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       11.07.90 sa
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(a-b,d-y)
      implicit complex*16(c,z)
      complex*16 spence
*
      eps1 = 1d-25
      eps = mib2*eps1
      cs = dcmplx(s,eps1*s)
      ct = dcmplx(t,-eps1*t)
      cmh2 = dcmplx(mib2,-eps)
      ctm = ct-mex2
      tm = dreal(ctm)
      tmth = t+mex2-mib2
      tmth2 = tmth*tmth
      st = s*t+tm*tm
      wur1 = dsqrt(tmth2-4d0*mex2/s*st)
      w12 = wur1*s
      y1 = (tmth+wur1)/2d0*s/st
      y2 = (tmth-wur1)/2d0*s/st
      cy1 = dcmplx(y1,eps1)
      cy2 = dcmplx(y2,-eps1)
      y11 = (y1-1d0)/y1
      y21 = (y2-1d0)/y2
      wur3 = dsqrt(tmth2-4d0*mex2*t)
      cy3 = dcmplx(tmth+wur3,eps)/2d0/t
      cy4 = dcmplx(tmth-wur3,-eps)/2d0/t
      cy31 = (cy3-1d0)/cy3
      cy41 = (cy4-1d0)/cy4
      cy13 = y1-cy3
      cy14 = y1-cy4
      cy23 = y2-cy3
      cy24 = y2-cy4
      mth = 2d0*mex2-mib2
      if (4d0*mex2.le.mib2) then
         wur5 = mib2*dsqrt(1d0-4d0*mex2/mib2)
         cy5 = dcmplx(mth+wur5,eps)/2d0/mex2
         cy6 = dcmplx(mth-wur5,-eps)/2d0/mex2
      else
         wur5 = mib2*dsqrt(4d0*mex2/mib2-1d0)
         cy5 = dcmplx(mth,wur5)/2d0/mex2
         cy6 = dcmplx(mth,-wur5)/2d0/mex2
      end if
      cy51 = (cy5-1d0)/cy5
      cy61 = (cy6-1d0)/cy6
      cy15 = y1-cy5
      cy16 = y1-cy6
      cy25 = y2-cy5
      cy26 = y2-cy6
      cy7 = s/(s+ctm)
      cy17 = y1-cy7
      cy27 = y2-cy7
      cd0 =
     1     -spence((cy1-1d0)/cy1)
     2     -spence(y1/cy17)+spence((y1-1d0)/cy17)
     3     -spence(y1/cy13)+spence((y1-1d0)/cy13)
     4     -spence(y1/cy14)+spence((y1-1d0)/cy14)
     5     +spence(y1/cy15)-spence((y1-1d0)/cy15)
     6     +spence(y1/cy16)-spence((y1-1d0)/cy16)
     7     +dlog((y1-1d0)/y1)*( dlog(tm*mex2/t/(s+tm))+cdlog(cy1)
     8     -cdlog(cy17)-cdlog(cy13*cy14)+cdlog(cy15*cy16) )
      cd0 = cd0-(
     1     -spence((cy2-1d0)/cy2)
     2     -spence(y2/cy27)+spence((y2-1d0)/cy27)
     3     -spence(y2/cy23)+spence((y2-1d0)/cy23)
     4     -spence(y2/cy24)+spence((y2-1d0)/cy24)
     5     +spence(y2/cy25)-spence((y2-1d0)/cy25)
     6     +spence(y2/cy26)-spence((y2-1d0)/cy26)
     7     +dlog((y2-1d0)/y2)*( dlog(tm*mex2/t/(s+tm))+cdlog(cy2)
     8     -cdlog(cy27)-cdlog(cy23*cy24)+cdlog(cy25*cy26) )
     9     )
*     
      beta = dsqrt(1d0-4d0*mex2/s)
      al1 = (1d0+beta)/2d0
      al2 = (1d0-beta)/2d0
      sat = al1*s+tm
      rz1 = (-beta*tm+mib2+wur1)/2d0/sat
      rz2 = (-beta*tm+mib2-wur1)/2d0/sat
      z1 = dcmplx(rz1,eps1)
      z2 = dcmplx(rz2,-eps1)
      z3 = cmplx(al1,eps1)
      z4 = cmplx(al2,-eps1)
      z13 = rz1-z3
      z14 = rz1-z4
      z23 = rz2-z3
      z24 = rz2-z4
      z5 = -cy5+1d0
      z6 = -cy6+1d0
      z15 = rz1-al2*z5
      z16 = rz1-al2*z6
      z25 = rz2-al2*z5
      z26 = rz2-al2*z6
      cd0 = cd0+(
     1     -spence((rz1-al2)/z15)+spence(rz1/z15)
     2     -spence((rz1-al2)/z16)+spence(rz1/z16)
     3     -spence(z1/(z1-al2))
     7     -cdlog(z1/(z1-al2))*( dlog(-mex2/tm/al2)-cdlog(al2-z1)
     8     +cdlog(z15*z16) )
     7     )
      cd0 = cd0-(
     1     -spence((rz2-al2)/z25)+spence(rz2/z25)
     2     -spence((rz2-al2)/z26)+spence(rz2/z26)
     3     -spence(z2/(z2-al2))
     7     -cdlog(z2/(z2-al2))*( dlog(-mex2/tm/al2)-cdlog(al2-z2)
     8     +cdlog(z25*z26) )
     7     )
      z15 = rz1+al1*z5
      z16 = rz1+al1*z6
      z25 = rz2+al1*z5
      z26 = rz2+al1*z6
      z17 = rz1+al1*(1d0-cy7)
      z27 = rz2+al1*(1d0-cy7)
      zst = dcmplx(s+tm,eps)
      cd0 = cd0+(
     1     +spence((rz1+al1)/z15)-spence(rz1/z15)
     2     +spence((rz1+al1)/z16)-spence(rz1/z16)
     3     -spence((rz1+al1)/z17)+spence(rz1/z17)
     7     +cdlog(z1/(z1+al1))*( cdlog(mex2/zst/al1)
     8     -cdlog(-z17)+cdlog(z15*z16) )
     7     )
      cd0 = cd0-(
     1     +spence((rz2+al1)/z25)-spence(rz2/z25)
     2     +spence((rz2+al1)/z26)-spence(rz2/z26)
     3     -spence((rz2+al1)/z27)+spence(rz2/z27)
     7     +cdlog(z2/(z2+al1))*( cdlog(mex2/zst/al1)
     8     -cdlog(-z27)+cdlog(z25*z26) )
     7     )
      z15 = rz1+z3
      z25 = rz2+z3
      z16 = z15-1d0
      z26 = z25-1d0
*     
      cd0 = cd0+(
     1     +spence(z15/z16)
     2     +spence(z16/(z16+al1))-spence(z15/(z16+al1))
     3     +spence(z16/z1)-spence(z15/z1)
     7     -cdlog(z15/z16)*( cdlog(-z16)-cdlog(-z1)
     8     -cdlog(-z16-al1) )
     7     )
      cd0 = cd0-(
     1     +spence(z25/z26)
     2     +spence(z26/(z26+al1))-spence(z25/(z26+al1))
     3     +spence(z26/z2)-spence(z25/z2)
     7     -cdlog(z25/z26)*( cdlog(-z26)-cdlog(-z2)
     8     -cdlog(-z26-al1) )
     7     )
*     
 999  ggttf = dreal(cd0)/w12
      return
      end
************************************************************************
      function eta(c1,c2)
************************************************************************
*     complex eta-function
*----------------------------------------------------------------------*
*     8.06.90    ansgar denner       last changed   11.07.94
************************************************************************
      implicit     none
      complex*16 eta,c1,c2
      real*8     im1,im2,im12,re1,re2
      real*8     pi

      pi = 4d0*datan(1d0)
      re1    = dreal(c1)
      re2    = dreal(c2)
      im1    = dimag(c1)
      im2    = dimag(c2)
      im12   = dimag(c1*c2)
 
      if(im1.lt.0d0.and.im2.lt.0d0.and.im12.gt.0d0) then
          eta = dcmplx(0d0,2d0*pi)
      else if (im1.gt.0d0.and.im2.gt.0d0.and.im12.lt.0d0) then
          eta = dcmplx(0d0,-2d0*pi)
      else
          eta = dcmplx(0d0,0d0)
          if(.not.(im2.eq.0d0.and.re2.gt.0d0.or.
     &             im1.eq.0d0.and.re1.gt.0d0).and.
     &       (im1.eq.0.and.re1.lt.0d0 .or.
     &        im2.eq.0.and.re2.lt.0d0 .or.
     &        im12.eq.0.and.dreal(c1*c2).lt.0d0)) then
             write(80,*) ' eta not defined '
             write(80,*) ' eta:  c1  = ',c1
             write(80,*) ' eta:  c2  = ',c2
             write(80,*) ' eta:  c12 = ',c1*c2
          end if
      end if
      end
