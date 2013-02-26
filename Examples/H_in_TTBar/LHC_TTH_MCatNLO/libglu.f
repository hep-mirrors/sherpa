      double precision function matb_gg(s,t1,t2,u1,u2,mt,mh)
      implicit none
      integer i,j
      real*8 s,s1,s2,t1,t2,u1,u2,mt,mh
      real*8 cs1,cs2,ct1,ct2,ct3,cu1,cu2,cu3
      real*8 ncf,cg,cf
      real*8 colfac(3,3),born(3,3),mat(3,68),matu(3,68),matus(3,68)
      real*8 pt(4),ptb(4),ph(4),p3(4),p4(4)

      matb_gg=0d0
c color factors:
      ncf=3d0
      cg=ncf
      cf=(ncf**2-1d0)/ncf/2d0
c s s*
      colfac(1,1)=cg**2*cf
c s t*
      colfac(1,2)=cg**2*cf/2d0
      colfac(2,1)=colfac(1,2)
c s u*
      colfac(1,3)=-cg**2*cf/2d0
      colfac(3,1)=colfac(1,3)
c t t*
      colfac(2,2)=cg*cf**2
c t u*
      colfac(2,3)=cg*cf**2-cg**2*cf/2d0
      colfac(3,2)=colfac(2,3)
c u u*      
      colfac(3,3)=cg*cf**2

      s1=s+t2+u2-mt**2
      s2=s+t1+u1-mt**2
      cs1=1d0/(s1-mt**2)/s
      cs2=1d0/(s2-mt**2)/s
      ct1=1d0/(s1-mt**2)/(t2-mt**2)
      ct2=1d0/(s2-mt**2)/(t1-mt**2)
      ct3=1d0/(t1-mt**2)/(t2-mt**2)
      cu1=1d0/(s1-mt**2)/(u2-mt**2)
      cu2=1d0/(s2-mt**2)/(u1-mt**2)
      cu3=1d0/(u1-mt**2)/(u2-mt**2)
      call matri(s,t1,t2,u1,u2,mt,mh,mat)
      call matri(s,u1,u2,t1,t2,mt,mh,matus)
      do i=1,68
         matu(1,i)=-matus(1,i)
         matu(2,i)=matus(3,i)
         matu(3,i)=matus(2,i)
      end do

      do i=1,3
c interference of s-channel Born with Born matrix element:

         born(1,i)=cs1*((s+u2-t2)*mat(i,1)-2d0*mat(i,2)+2d0*mat(i,18)+
     $        2d0*mat(i,21)+2d0*mat(i,24)+2d0*mat(i,15)+4d0*mat(i,66)-
     $        4d0*mat(i,62)-4d0*mat(i,60))+
     $        cs2*((s+u1-t1)*mat(i,1)-2d0*mat(i,2)+2d0*mat(i,18)+
     $        2d0*mat(i,21)+2d0*mat(i,24)+2d0*mat(i,15)+4d0*mat(i,61)-
     $        4d0*mat(i,63)-4d0*mat(i,60))

c interference of t-channel Born with Born matrix element:

         born(2,i)=ct1*((mt**2-t2)*mat(i,5)-mat(i,6)+2d0*mat(i,24)+
     $        2d0*mat(i,17)-2d0*mat(i,26)-4d0*mat(i,62)+4d0*mat(i,68))+
     $        ct2*((mt**2-t1)*mat(i,5)-mat(i,6)+2d0*mat(i,15)+
     $        2d0*mat(i,25)-2d0*mat(i,16)-4d0*mat(i,63)+4d0*mat(i,64))+
     $        ct3*(2d0*mat(i,25)+2d0*mat(i,17)-mat(i,6)-
     $        4d0*mat(i,65))

c interference of u-channel Born with Born matrix element:

         born(3,i)=cu1*((mt**2-u2)*matu(i,5)-matu(i,6)+2d0*matu(i,24)+
     $        2d0*matu(i,17)-2d0*matu(i,26)-4d0*matu(i,62)+
     $        4d0*matu(i,68))+
     $        cu2*((mt**2-u1)*matu(i,5)-matu(i,6)+
     $        2d0*matu(i,15)+2d0*matu(i,25)-2d0*matu(i,16)-
     $        4d0*matu(i,63)+4d0*matu(i,64))+
     $        cu3*(2d0*matu(i,25)+2d0*matu(i,17)-matu(i,6)-
     $        4d0*matu(i,65))
      end do
      do i=1,3
         do j=1,3
            matb_gg=matb_gg+colfac(i,j)*born(i,j)/64d0
         end do
      end do
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      double precision function matv_gg(s,t1,t2,u1,u2,mt,mh)
      implicit none
      integer i,j
      real*8 s,s1,s2,t1,t2,u1,u2,mt,mh
      real*8 cs1,cs2,ct1,ct2,ct3,cu1,cu2,cu3
      real*8 ncf,cg,cf
      real*8 colfac(3,3),mv(3,3),mat(3,68),matu(3,68),matus(3,68)
      real*8 colfac_b(3,3),mvt(3,3),mvu(3,3)
      real*8 matv4(3)
      complex*16 kmvs1(68),kmvs2(68),kmvs3(68),kmvs4(68),kmvs5(68)
      complex*16 kmvt1(68),kmvt2(68),kmvt3(68),kmvt4(68),kmvt5(68)
      complex*16 kmvu1(68),kmvu2(68),kmvu3(68),kmvu4(68),kmvu5(68)
      complex*16 kmvb1(3,68),kmvb2(3,68),kmvub1(3,68),kmvub2(3,68),zero
      complex*16 kmvtp2(68),kmvup2(68),kmvtp3(68),kmvup3(68)
      complex*16 kmvtp4(68),kmvup4(68),kmvtp1(68),kmvup1(68)
      complex*16 kmvtp5(68),kmvup5(68),kmvtp6(68),kmvup6(68)
      complex*16 mvtp1(8),mvup1(8),mvtp2(8),mvup2(8),mvtp3(8),mvup3(8)
      complex*16 mvtp4(8),mvup4(8),mvtp5(8),mvup5(8),mvtp6(8),mvup6(8)
      complex*16 mvtp(8),mvup(8)
      real*8 pt(4),ptb(4),ph(4),p3(4),p4(4)
      real*8 dets3

      dets3=1d0

      matv_gg=0d0
      zero=dcmplx(0d0,0d0)
c initialization of coefficients:
      do i=1,68
         kmvs1(i)=zero
         kmvs2(i)=zero
         kmvs3(i)=zero
         kmvs4(i)=zero
         kmvs5(i)=zero
         kmvt1(i)=zero
         kmvt2(i)=zero
         kmvt3(i)=zero
         kmvt4(i)=zero
         kmvt5(i)=zero
         kmvu1(i)=zero
         kmvu2(i)=zero
         kmvu3(i)=zero
         kmvu4(i)=zero
         kmvu5(i)=zero
         kmvtp1(i)=zero
         kmvtp2(i)=zero
         kmvtp3(i)=zero
         kmvtp4(i)=zero
         kmvtp5(i)=zero
         kmvtp6(i)=zero
         kmvup1(i)=zero
         kmvup2(i)=zero
         kmvup3(i)=zero
         kmvup4(i)=zero
         kmvup5(i)=zero
         kmvup6(i)=zero
         do j=1,3
            kmvb1(j,i)=zero
            kmvb2(j,i)=zero
            kmvub1(j,i)=zero
            kmvub2(j,i)=zero
         enddo
      enddo
      do i=1,8
         mvtp1(i)=zero
         mvtp2(i)=zero
         mvtp3(i)=zero
         mvtp4(i)=zero
         mvtp5(i)=zero
         mvtp6(i)=zero
         mvup1(i)=zero
         mvup2(i)=zero
         mvup3(i)=zero
         mvup4(i)=zero
         mvup5(i)=zero
         mvup6(i)=zero
         mvtp(i)=zero
         mvup(i)=zero
      enddo
c color factors:
      ncf=3d0
      cg=ncf
      cf=(ncf**2-1d0)/ncf/2d0
c s s*
      colfac(1,1)=cg**2*cf
c s t*
      colfac(1,2)=cg**2*cf/2d0
      colfac(2,1)=colfac(1,2)
c s u*
      colfac(1,3)=-cg**2*cf/2d0
      colfac(3,1)=colfac(1,3)
c t t*
      colfac(2,2)=cg*cf**2
c t u*
      colfac(2,3)=cg*cf**2-cg**2*cf/2d0
      colfac(3,2)=colfac(2,3)
c u u*      
      colfac(3,3)=cg*cf**2
c
      s1=s+t2+u2-mt**2
      s2=s+t1+u1-mt**2
      cs1=1d0/(s1-mt**2)/s
      cs2=1d0/(s2-mt**2)/s
      ct1=1d0/(s1-mt**2)/(t2-mt**2)
      ct2=1d0/(s2-mt**2)/(t1-mt**2)
      ct3=1d0/(t1-mt**2)/(t2-mt**2)
      cu1=1d0/(s1-mt**2)/(u2-mt**2)
      cu2=1d0/(s2-mt**2)/(u1-mt**2)
      cu3=1d0/(u1-mt**2)/(u2-mt**2)
      call matri(s,t1,t2,u1,u2,mt,mh,mat)
      call matri(s,u1,u2,t1,t2,mt,mh,matus)
      do i=1,68
         matu(1,i)=-matus(1,i)
         matu(2,i)=matus(3,i)
         matu(3,i)=matus(2,i)
      end do
c      if (mt.gt.0d0) goto 999
c Form factors which multiply the SMEs:
c s-channel: 
c Higgs radiated from top: kmvs1(68)
c Higgs radiated from anti-top: kmvs2(68)
c kmvs3(68): UV+IR finite contributions (B1,V4,B3, parts of B4,V6)
c kmvs4(68): UV+IR finite contributions (parts of B4,V6)
c kmvs5(68): UV+IR finite contributions (parts of B4,V6)
c
      call gluform_s(s,t1,t2,u1,u2,mt,mh,kmvs1,kmvs2,kmvs3,kmvs4,kmvs5)
c      do i=1,68
c         kmvs1(i)=kmvs1(i)*dets3
c         kmvs2(i)=kmvs2(i)*dets3
c         kmvs3(i)=kmvs3(i)*dets3
c         kmvs4(i)=kmvs4(i)*dets3
c         kmvs5(i)=kmvs5(i)*dets3
c      enddo
c
c t-channel: 
c Higgs radiated from top: kmvt1(68)
c Higgs radiated from anti-top: kmvt2(68)
c Higgs radiated from virtual top: kmvt3(68)
c kmvt4(68): UV+IR finite contributions (V10,B5,B6)
c kmvt5(68): UV+IR finite contributions (parts of B9)
c
      call gluform_t(s,t1,t2,u1,u2,mt,mh,kmvt1,kmvt2,kmvt3,kmvt4,kmvt5)
c      do i=1,68
c         kmvt1(i)=kmvt1(i)*dets3
c         kmvt2(i)=kmvt2(i)*dets3
c         kmvt3(i)=kmvt3(i)*dets3
c         kmvt4(i)=kmvt4(i)*dets3
c         kmvt5(i)=kmvt5(i)*dets3
c      enddo
c
c Pentagons (P2,P3,P4):
c
c      call pentp2_t(s,t1,t2,u1,u2,mt,mh,kmvtp2)
c      call pentp3_t(s,t1,t2,u1,u2,mt,mh,kmvtp3)
c      call pentp4_t(s,t1,t2,u1,u2,mt,mh,kmvtp4)
c
c u-channel: 
c Higgs radiated from top: kmvu1(68)
c Higgs radiated from anti-top: kmvu2(68)
c Higgs radiated from virtual top: kmvu3(68)
c kmvu4(68): UV+IR finite contributions (V10,B5,B6)
c kmvu5(68): UV+IR finite contributions (parts of B9)
c
      call gluform_t(s,u1,u2,t1,t2,mt,mh,kmvu1,kmvu2,kmvu3,kmvu4,kmvu5)
c      do i=1,68
c         kmvu1(i)=kmvu1(i)*dets3
c         kmvu2(i)=kmvu2(i)*dets3
c         kmvu3(i)=kmvu3(i)*dets3
c         kmvu4(i)=kmvu4(i)*dets3
c         kmvu5(i)=kmvu5(i)*dets3
c      enddo
c
c Pentagons (P2,P3,P4):
c
c      call pentp2_t(s,u1,u2,t1,t2,mt,mh,kmvup2)
c      call pentp3_t(s,u1,u2,t1,t2,mt,mh,kmvup3)
c      call pentp4_t(s,u1,u2,t1,t2,mt,mh,kmvup4)
c
c interference with the Born matrix element:
      do i=1,3
         do j=1,3
            mv(i,j)=0d0
         enddo
      enddo   
      do i=1,3
c interference of s-channel contribution with Born matrix element:
         do j=1,68
            mv(1,i)=mv(1,i)+2d0*(cs1*dreal(kmvs1(j))+
     $           cs2*dreal(kmvs2(j))+dreal(kmvs3(j)))*mat(i,j)
         enddo
c interference of t-channel contribution with Born matrix element:
         do j=1,68
            mv(2,i)=mv(2,i)+2d0*(ct1*dreal(kmvt1(j))+
     $           ct2*dreal(kmvt2(j))+ct3*dreal(kmvt3(j))+
     $           dreal(kmvt4(j))+(cf-cg/2d0)*
     $           (dreal(kmvt5(j))+dreal(kmvtp2(j))+
     $           dreal(kmvtp3(j))+dreal(kmvtp4(j))))*mat(i,j)
         enddo
c interference of u-channel contribution with Born matrix element:
         do j=1,68
            mv(3,i)=mv(3,i)+2d0*(cu1*dreal(kmvu1(j))+
     $           cu2*dreal(kmvu2(j))+cu3*dreal(kmvu3(j))+
     $           dreal(kmvu4(j))+(cf-cg/2d0)*(dreal(kmvu5(j))+
     $           dreal(kmvup2(j))+dreal(kmvup3(j))+
     $           dreal(kmvup4(j))))*matu(i,j)
         enddo
      end do
c
c new calculation of P2,P3,P4:
c      write(6,*)'old',mv(2,1),mv(2,2),mv(2,3)
      call newp2_t(s,t1,t2,u1,u2,mt,mh,mvtp2)
      call newp3_t(s,t1,t2,u1,u2,mt,mh,mvtp3)
      call newp4_t(s,t1,t2,u1,u2,mt,mh,mvtp4)
      do i=1,8
         mvtp(i)=mvtp2(i)+mvtp3(i)+mvtp4(i)
      enddo
      mv(2,1)=mv(2,1)+2d0*(cf-cg/2d0)*
     $     (cs1*dreal(mvtp(1))+cs2*dreal(mvtp(2)))
      mv(2,2)=mv(2,2)+2d0*(cf-cg/2d0)*
     $     (ct1*dreal(mvtp(3))+ct2*dreal(mvtp(4))+
     $     ct3*dreal(mvtp(5)))
      mv(2,3)=mv(2,3)+2d0*(cf-cg/2d0)*
     $     (cu1*dreal(mvtp(6))+cu2*dreal(mvtp(7))+cu3*dreal(mvtp(8)))
c      write(6,*)'new',mv(2,1),mv(2,2),mv(2,3)
c      write(6,*)'old',mv(3,1),mv(3,2),mv(3,3)
      call newp2_t(s,u1,u2,t1,t2,mt,mh,mvup2)
      call newp3_t(s,u1,u2,t1,t2,mt,mh,mvup3)
      call newp4_t(s,u1,u2,t1,t2,mt,mh,mvup4)
      do i=1,8
         mvup(i)=mvup2(i)+mvup3(i)+mvup4(i)
      enddo
      mv(3,1)=mv(3,1)-2d0*(cf-cg/2d0)*
     $     (cs1*dreal(mvup(1))+cs2*dreal(mvup(2)))
      mv(3,2)=mv(3,2)+2d0*(cf-cg/2d0)*
     $     (ct1*dreal(mvup(6))+ct2*dreal(mvup(7))+
     $     ct3*dreal(mvup(8)))
      mv(3,3)=mv(3,3)+2d0*(cf-cg/2d0)*
     $     (cu1*dreal(mvup(3))+cu2*dreal(mvup(4))+cu3*dreal(mvup(5)))
c      write(6,*)'new',mv(3,1),mv(3,2),mv(3,3)
      do i=1,3
         do j=1,3
            matv_gg=matv_gg+colfac(i,j)*mv(i,j)/64d0
         end do
      end do
 999  continue
c      if(mt.gt.0d0) goto 9999
c contribution from B7,B8,B10,P1,P5,P6 (and parts of B4,V6,B9,P2):
c t-channel
      call kmvbox_ir(s,t1,t2,u1,u2,mt,mh,kmvb1,kmvb2)
c      do i=1,68
c         do j=1,3
c            kmvb1(j,i)=kmvb1(j,i)*dets3
c            kmvb2(j,i)=kmvb2(j,i)*dets3
c         enddo
c      enddo
c      call pentp1_t(s,t1,t2,u1,u2,mt,mh,kmvtp1)
c      call pentp5_t(s,t1,t2,u1,u2,mt,mh,kmvtp5)
c      call pentp6_t(s,t1,t2,u1,u2,mt,mh,kmvtp6)
      do i=1,3
         do j=1,3
            mvt(i,j)=0d0
         enddo
      enddo
      do i=1,3
         do j=1,68
            mvt(1,i)=mvt(1,i)+2d0*(dreal(kmvb1(1,j))/(t2-mt**2)+
     $           dreal(kmvb2(1,j))/(t1-mt**2))*mat(i,j)
            mvt(2,i)=mvt(2,i)+2d0*(dreal(kmvb1(2,j))/(s1-mt**2)+
     $           dreal(kmvb2(2,j))/(s2-mt**2))*mat(i,j)
c add contribution from B4,V6:
            mvt(2,i)=mvt(2,i)+2d0*dreal(kmvs4(j))*mat(i,j)
c add contribution from P1:
            mvt(2,i)=mvt(2,i)+2d0*dreal(kmvtp1(j))*mat(i,j)
            mvt(3,i)=mvt(3,i)+2d0*(dreal(kmvb1(3,j))/(s1-mt**2)+
     $           dreal(kmvb2(3,j))/(s2-mt**2))*mat(i,j)
c add contribution from B9,P2,P3,P4:
            mvt(3,i)=mvt(3,i)+2d0*(dreal(kmvt5(j))+
     $           dreal(kmvtp2(j))+dreal(kmvtp3(j))+
     $           dreal(kmvtp4(j)))*cg/2d0*mat(i,j)
c add contribution from P5,P6:
            mvt(3,i)=mvt(3,i)+2d0*(dreal(kmvtp5(j))+
     $           dreal(kmvtp6(j)))*mat(i,j)
        enddo
      enddo
c new calculation of P1,P5,P6:
c      write(6,*)'old',mvt(3,1),mvt(3,2),mvt(3,3)
      call newp1_t(s,t1,t2,u1,u2,mt,mh,mvtp1)
      call newp5_t(s,t1,t2,u1,u2,mt,mh,mvtp5)
      call newp6_t(s,t1,t2,u1,u2,mt,mh,mvtp6)
      mvt(2,1)=mvt(2,1)+2d0*(cs1*dreal(mvtp1(1))+cs2*dreal(mvtp1(2)))
      mvt(2,2)=mvt(2,2)+2d0*(ct1*dreal(mvtp1(3))+ct2*dreal(mvtp1(4))+
     $     ct3*dreal(mvtp1(5)))
      mvt(2,3)=mvt(2,3)+2d0*(cu1*dreal(mvtp1(6))+cu2*dreal(mvtp1(7))+
     $     cu3*dreal(mvtp1(8)))
c contribution from P5,P6 and rest of P2,P3,P4:
      do i=1,8
         mvtp(i)=cg/2d0*mvtp(i)+mvtp5(i)+mvtp6(i)
      enddo
      mvt(3,1)=mvt(3,1)+2d0*(cs1*dreal(mvtp(1))+cs2*dreal(mvtp(2)))
      mvt(3,2)=mvt(3,2)+2d0*(ct1*dreal(mvtp(3))+ct2*dreal(mvtp(4))+
     $     ct3*dreal(mvtp(5)))
      mvt(3,3)=mvt(3,3)+2d0*(cu1*dreal(mvtp(6))+cu2*dreal(mvtp(7))+
     $     cu3*dreal(mvtp(8)))
c      write(6,*)'new',mvt(3,1),mvt(3,2),mvt(3,3)
c B7: 
c t s*      
      colfac_b(1,1)=cg**2*cf/2d0
c t t*
      colfac_b(1,2)=cg*cf**2
c t u*
      colfac_b(1,3)=cg*cf**2-cg**2*cf/2d0
c B8: 
c t s*      
      colfac_b(2,1)=cg**2*cf/2d0
c t t*
      colfac_b(2,2)=cg*cf**2+cf*cg/6d0
c t u*
      colfac_b(2,3)=cg*cf**2-cg**2*cf/2d0+cf*cg/6d0
c B10: 
c t s*      
      colfac_b(3,1)=0d0
c t t*
      colfac_b(3,2)=cf*cg/6d0
c t u*
      colfac_b(3,3)=cf*cg/6d0
      do j=1,3
         do i=1,3
            matv_gg=matv_gg+colfac_b(i,j)*mvt(i,j)/64d0
         enddo
      enddo
c u-channel
      call kmvbox_ir(s,u1,u2,t1,t2,mt,mh,kmvub1,kmvub2)
c      do i=1,68
c         do j=1,3
c            kmvub1(j,i)=kmvub1(j,i)*dets3
c            kmvub2(j,i)=kmvub2(j,i)*dets3
c         enddo
c      enddo
c      call pentp1_t(s,u1,u2,t1,t2,mt,mh,kmvup1)
c      call pentp5_t(s,u1,u2,t1,t2,mt,mh,kmvup5)
c      call pentp6_t(s,u1,u2,t1,t2,mt,mh,kmvup6)
      do i=1,3
         do j=1,3
            mvu(i,j)=0d0
         enddo
      enddo   
      do i=1,3
         do j=1,68
            mvu(1,i)=mvu(1,i)+2d0*(dreal(kmvub1(1,j))/(u2-mt**2)+
     $           dreal(kmvub2(1,j))/(u1-mt**2))*matu(i,j)
            mvu(2,i)=mvu(2,i)+2d0*(dreal(kmvub1(2,j))/(s1-mt**2)+
     $           dreal(kmvub2(2,j))/(s2-mt**2))*matu(i,j)
c add contribution from B4,V6 (with s-channel SMEs !):
            mvu(2,i)=mvu(2,i)+2d0*dreal(kmvs5(j))*mat(i,j)
c add contribution from P1:
            mvu(2,i)=mvu(2,i)+2d0*dreal(kmvup1(j))*matu(i,j)
            mvu(3,i)=mvu(3,i)+2d0*(dreal(kmvub1(3,j))/(s1-mt**2)+
     $           dreal(kmvub2(3,j))/(s2-mt**2))*matu(i,j)
c add contribution from B9,P2,P3,P4:
            mvu(3,i)=mvu(3,i)+2d0*(dreal(kmvu5(j))+
     $           dreal(kmvup2(j))+dreal(kmvup3(j))+
     $           dreal(kmvup4(j)))*cg/2d0*matu(i,j)
c add contribution from P5,P6:
            mvu(3,i)=mvu(3,i)+2d0*(dreal(kmvup5(j))+
     $           dreal(kmvup6(j)))*matu(i,j)
         enddo
      enddo
c
c new calculation of P1,P5,P6:
c      write(6,*)'old',mvu(3,1),mvu(3,2),mvu(3,3)
      call newp1_t(s,u1,u2,t1,t2,mt,mh,mvup1)
      call newp5_t(s,u1,u2,t1,t2,mt,mh,mvup5)
      call newp6_t(s,u1,u2,t1,t2,mt,mh,mvup6)
c contribution from P1:
      mvu(2,1)=mvu(2,1)-2d0*(cs1*dreal(mvup1(1))+cs2*dreal(mvup1(2)))
      mvu(2,2)=mvu(2,2)+2d0*(ct1*dreal(mvup1(6))+ct2*dreal(mvup1(7))+
     $     ct3*dreal(mvup1(8)))
      mvu(2,3)=mvu(2,3)+2d0*(cu1*dreal(mvup1(3))+cu2*dreal(mvup1(4))+
     $     cu3*dreal(mvup1(5)))
c contribution from P5,P6 and rest of P2,P3,P4:
      do i=1,8
         mvup(i)=cg/2d0*mvup(i)+mvup5(i)+mvup6(i)
      enddo
      mvu(3,1)=mvu(3,1)-2d0*(cs1*dreal(mvup(1))+cs2*dreal(mvup(2)))
      mvu(3,2)=mvu(3,2)+2d0*(ct1*dreal(mvup(6))+ct2*dreal(mvup(7))+
     $     ct3*dreal(mvup(8)))
      mvu(3,3)=mvu(3,3)+2d0*(cu1*dreal(mvup(3))+cu2*dreal(mvup(4))+
     $     cu3*dreal(mvup(5)))
c      write(6,*)'new',mvu(3,1),mvu(3,2),mvu(3,3)
c B7: 
c u s*      
      colfac_b(1,1)=-cg**2*cf/2d0
c u t*
      colfac_b(1,2)=cg*cf**2-cg**2*cf/2d0
c u u*
      colfac_b(1,3)=cg*cf**2
c B8: 
c u s*      
      colfac_b(2,1)=-cg**2*cf/2d0
c u t*
      colfac_b(2,2)=cg*cf**2-cg**2*cf/2d0+cf*cg/6d0
c u u*
      colfac_b(2,3)=cg*cf**2+cf*cg/6d0
c B10: 
c u s*      
      colfac_b(3,1)=0d0
c u t*
      colfac_b(3,2)=cf*cg/6d0
c u u*
      colfac_b(3,3)=cf*cg/6d0
      do j=1,3
         do i=1,3
            matv_gg=matv_gg+colfac_b(i,j)*mvu(i,j)/64d0
         enddo
      enddo
 9999 continue
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine gluform_s(s,t1,t2,u1,u2,mt,mh,kmv1,kmv2,kmv3,
     $     kmv4,kmv5)
      implicit none
      integer i,j
      real*8 s,s3,p3p4,p1p2,p1p3,p1p4,p2p3,p2p4
      real*8 t1,t2,u1,u2,mt,mh,s1,s2
      complex*16 zero,kmv1(68),kmv2(68),kmvf1(6,68),kmvf2(6,68)
      complex*16 kmv3(68),kmvf3(15,68),kmv4(68),kmv5(68)
      complex*16 fvf1(6),fvf2(6),fh1(3),fh2(3),fh3(3),
     $     fs1(2),fs2(2),fs3(2),fs4(2)
      complex*16 fvb1(6),fvb2(6),kmvs(16,68)
      complex*16 sigmagg,vertexgg

      zero=dcmplx(0d0,0d0)
      do i=1,68
         kmv1(i)=zero
         kmv2(i)=zero
         kmv3(i)=zero
         kmv4(i)=zero
         kmv5(i)=zero
         do j=1,6
            kmvf1(j,i)=zero
            kmvf2(j,i)=zero
         enddo
         do j=1,15
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
      
      kmvf1(1,1)=2d0*fvf1(1)*(p1p3-p1p4+p3p4)-4d0*mt*p3p4*fvf1(2)
     $     +2d0*mt*fvf1(3)*(p1p4-p1p3)+2d0*fvf1(5)*(p1p3**2-p1p4**2
     $     -p1p3*p3p4+p1p4*p3p4)
      kmvf1(1,2)=2d0*(2d0*mt*fvf1(2)-fvf1(1))
      kmvf1(1,3)=2d0*fvf1(2)*(-p3p4+2d0*p1p3)+fvf1(3)*(p1p4-p1p3)
      kmvf1(1,4)=2d0*fvf1(2)*(p3p4-2d0*p1p4)+fvf1(3)*(p1p4-p1p3)

      kmvf2(1,1)=2d0*fvf2(1)*(-p2p3+p2p4+p3p4)-4d0*mt*p3p4*fvf2(2)
     $     +2d0*mt*fvf2(3)*(p2p4-p2p3)+2d0*fvf2(5)*(-p2p3**2+p2p4**2
     $     +p2p3*p3p4-p2p4*p3p4)
      kmvf2(1,2)=2d0*(2d0*mt*fvf2(2)-fvf2(1))
      kmvf2(1,3)=2d0*fvf2(2)*(-p3p4+2d0*p2p3)+fvf2(3)*(-p2p4+p2p3)
      kmvf2(1,4)=2d0*fvf2(2)*(p3p4-2d0*p2p4)+fvf2(3)*(-p2p4+p2p3)

c contribution from V3(1) and V3(2):
c fh3 does not contribute (only contributes to the t,u-channel)
      call formh(s,t1,t2,u1,u2,mt,mh,fh1,fh2,fh3)
      
      kmvf1(2,1)=2d0*(p3p4+p1p3-p1p4)*(fh1(1)+mt*(fh1(2)+fh1(3)))
      kmvf1(2,2)=-2d0*(fh1(1)+mt*(fh1(2)+fh1(3)))
      kmvf1(2,3)=2d0*fh1(2)*(p3p4-p1p3-p1p4)
      kmvf1(2,4)=2d0*fh1(2)*(-p3p4+p1p3+p1p4)

      kmvf2(2,1)=2d0*(p3p4+p2p4-p2p3)*(fh2(1)+mt*(fh2(2)+fh2(3)))
      kmvf2(2,2)=-2d0*(fh2(1)+mt*(fh2(2)+fh2(3)))
      kmvf2(2,3)=2d0*fh2(3)*(p3p4-p2p4-p2p3)
      kmvf2(2,4)=2d0*fh2(3)*(-p3p4+p2p4+p2p3)

c contribution from S2(1) and S2(2):
c fs3,fs4 do not contribute (they only contribute to the t,u channel)
      call forms(s,t1,t2,u1,u2,mt,mh,fs1,fs2,fs3,fs4)

      kmvf1(3,1)=2d0*(p3p4+p1p3-p1p4)*(fs1(1)+mt*fs1(2))
      kmvf1(3,2)=-2d0*(fs1(1)+mt*fs1(2))
      kmvf1(3,3)=2d0*fs1(2)*(p3p4-p1p3-p1p4)
      kmvf1(3,4)=2d0*fs1(2)*(-p3p4+p1p3+p1p4)

      kmvf2(3,1)=2d0*(p3p4+p2p4-p2p3)*(fs2(1)+mt*fs2(2))
      kmvf2(3,2)=-2d0*(fs2(1)+mt*fs2(2))
      kmvf2(3,3)=2d0*fs2(2)*(p3p4-p2p4-p2p3)
      kmvf2(3,4)=2d0*fs2(2)*(-p3p4+p2p4+p2p3)
      do i=1,4
         kmvf1(3,i)=-1d0/(s+t2+u2-2d0*mt**2)*kmvf1(3,i)
         kmvf2(3,i)=-1d0/(s+t1+u1-2d0*mt**2)*kmvf2(3,i)
      enddo

c contribution from S1(1) and S1(2):
c i sigma_mu_nu=i (q_mu q_nu - q^2 g_mu_nu) sigmagg 
      call ggself(s,mt,mh,sigmagg)
      kmvf1(4,1)=2d0*(p3p4+p1p3-p1p4)*(-sigmagg)
      kmvf1(4,2)=-2d0*(-sigmagg)

      kmvf2(4,1)=2d0*(p3p4+p2p4-p2p3)*(-sigmagg)
      kmvf2(4,2)=-2d0*(-sigmagg)

c contribution from V1(1) and V1(2):
c -g_s f_abc Lambda_mu_nu_ro:
      call ggvert(s,mt,mh,vertexgg)
      kmvf1(5,1)=2d0*(p3p4+p1p3-p1p4)*(vertexgg)
      kmvf1(5,2)=-2d0*(vertexgg)

      kmvf2(5,1)=2d0*(p3p4+p2p4-p2p3)*(vertexgg)
      kmvf2(5,2)=-2d0*(vertexgg)

c contribution from B2(1) and B2(2):
      call formb2(s,t1,t2,u1,u2,mt,mh,fvb1,fvb2)
c B2(1)
      kmvf1(6,1)=-2d0*p3p4*fvb1(2)+(p1p4-p1p3)*fvb1(3)
     $     +(p2p4-p2p3)*fvb1(4)
      kmvf1(6,2)=2d0*fvb1(2)
      kmvf1(6,3)=fvb1(1)+(p1p4-p1p3)*fvb1(5)
     $     +(p2p4-p2p3)*fvb1(6)
      kmvf1(6,4)=-fvb1(1)+(p1p4-p1p3)*fvb1(5)
     $     +(p2p4-p2p3)*fvb1(6)
c B2(2)
      kmvf2(6,1)=-2d0*p3p4*fvb2(2)+(p1p4-p1p3)*fvb2(3)
     $     +(p2p4-p2p3)*fvb2(4)
      kmvf2(6,2)=2d0*fvb2(2)
      kmvf2(6,3)=fvb2(1)+(p1p4-p1p3)*fvb2(5)
     $     +(p2p4-p2p3)*fvb2(6)
      kmvf2(6,4)=-fvb2(1)+(p1p4-p1p3)*fvb2(5)
     $     +(p2p4-p2p3)*fvb2(6)
      do i=1,4
         kmvf1(6,i)=kmvf1(6,i)*(s1-mt**2)
         kmvf2(6,i)=kmvf2(6,i)*(s2-mt**2)
      enddo
c UV+IR finite s-channel contributions:
      call kmvs_finite(s,t1,t2,u1,u2,mt,mh,kmvs)
c contribution from B1(1-6)
      do i=1,68
         kmvf3(1,i)=-2d0*kmvs(1,i)/(2d0*mt**2+2d0*p1p2)
         kmvf3(2,i)=-2d0*kmvs(2,i)/(2d0*mt**2+2d0*p1p2)
         kmvf3(3,i)=-2d0*kmvs(3,i)/(2d0*mt**2+2d0*p1p2)
      enddo
c contribution from V4(1-4)
      do i=1,68
         kmvf3(4,i)=2d0*kmvs(4,i)/(2d0*mt**2+2d0*p1p2)/
     $        (2d0*mt**2+2d0*p1p2-2d0*p1p4-2d0*p2p4)
         kmvf3(5,i)=2d0*kmvs(5,i)/(2d0*mt**2+2d0*p1p2)/
     $        (2d0*mt**2+2d0*p1p2-2d0*p1p3-2d0*p2p3)
      enddo
c contribution from B3
      do i=1,68
         kmvf3(6,i)=kmvs(6,i)/s
      enddo
c contribution from B4,V6 (part with s-channel color structure)
c B4:7,8,9
c V6(1): 10,11,12;  V6(2): 13,14,15  
      do i=1,68
         kmvf3(7,i)=kmvs(7,i)
         kmvf3(10,i)=kmvs(10,i)/(s1-mt**2)
         kmvf3(13,i)=kmvs(13,i)/(s2-mt**2)
      enddo
c parts with t(u) channel color structure
      do i=1,68
         kmv4(i)=kmvs(8,i)+kmvs(11,i)/(s1-mt**2)+
     $        kmvs(14,i)/(s2-mt**2)
         kmv5(i)=-kmvs(9,i)-kmvs(12,i)/(s1-mt**2)-
     $        kmvs(15,i)/(s2-mt**2)
      enddo
c contribution from V5(1-2):
      do i=1,68
         kmvf3(14,i)=2d0*kmvs(16,i)/(2d0*mt**2+2d0*p1p2)/s
      enddo
c sum all contributions:
      do i=1,68
         do j=1,6
            kmv1(i)=kmv1(i)+kmvf1(j,i)
            kmv2(i)=kmv2(i)+kmvf2(j,i)
         enddo
         do j=1,14
            kmv3(i)=kmv3(i)+kmvf3(j,i)
         enddo
      enddo
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine formb2(s,t1,t2,u1,u2,mt,mh,fv1,fv2)
c form factors for the s-channel UV finite, but IR div. contribution
c box diagram B2(1) and B2(2)
      implicit none
      integer i
      real*8 ncf,cg,cf,mt,mh
      real*8 s,t1,t2,u1,u2,s3,p3p4,p1p2,p1p3,p1p4,p2p3,p2p4,p1pq,p2pq
      real*8 d1r(0:3),d2r(0:6),d3r(0:3,0:3),ch1(0:2),s1,s2
      complex*16 zero,fv1(6),fv2(6)
      complex*16 c0,d0,D0_ir,c11,c12,c20,c21,c22,c23
      complex*16 d1(0:3),d2(0:6),d3(0:3,0:3),c1(0:2,4),c2(0:3,4)
      complex*16 d11,d12,d13,d21,d22,d23,d20,d212,d213,d223
      real*8 lambda,lambda2

      COMMON/IR/LAMBDA,LAMBDA2

      lambda=mt
      lambda2=lambda**2
      zero=dcmplx(0d0,0d0)
      do i=1,6
         fv1(i)=zero
         fv2(i)=zero
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

      d1(0)=D0_ir(mt,mt**2,mh**2,s,mt**2,s1,
     $     2d0*mt**2+2d0*p1p2,mt,mt,mt)
      call ddfunc_ir(mt,mt**2,mt**2,s1,-p1p2,mt**2-p1pq,
     $     -p1p2+p2pq,0d0,mt,mt,mt,d1,d2,d3,c1,c2)
c      call dmn(mt**2,2d0*(mt**2+p1p2),mh**2,mt**2,
c     $     3d0*mt**2+mh**2-2d0*p1pq+2d0*p1p2,
c     $     s1,0d0,mt,mt,mt,d1r,d2r,d3r,ch1)
c      write(6,*)'1',dreal(d1(0)),d1r(0)
c      write(6,*)'1',dreal(c1(0,1)),ch1(0)
c      write(6,*)'1',dreal(c1(1,1)),-ch1(0)-ch1(1)-ch1(2)
c      write(6,*)'1',dreal(d2(0)),d2r(0)
c      write(6,*)'1',dreal(d2(1)),d2r(1)
c      write(6,*)'1',dreal(d2(2)),d2r(2)
c      write(6,*)'1',dreal(d2(3)),d2r(3)
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
c B2(1)      
      fv1(1)=2d0*mt*(c0-4d0*d20+2d0*(mt**2+p2pq)*d11+
     $     (-2d0*p1p2+2d0*p1pq)*d12+
     $     (2d0*mt**2-2d0*p1pq+2d0*p2pq)*d13)
      fv1(2)=4d0*p1p2*(d0+d11+d12+d13)
      fv1(3)=4d0*(-c0+c11+2d0*p1p2*(d0+2d0*d11+d12)+
     $     2d0*(mt**2+p1p2)*d21+2d0*(mt**2+p1p2-p2pq)*d23+
     $     2d0*(2d0*mt**2+2d0*p1p2-p2pq)*d213+
     $     (4d0*p1p2+2d0*p1pq-2d0*p2pq-s)*d13)
      fv1(4)=4d0*(-c0+c11+2d0*d20-2d0*mt**2*d11-
     $     2d0*(mt**2+p1p2)*d212+2d0*(-mt**2-p1p2+p2pq)*d223+
     $     (-2d0*mt**2+2d0*p1pq-s)*d13)
      fv1(5)=4d0*mt*(-d12+2d0*(d13+d213+d23))
      fv1(6)=-4d0*mt*(d11+d13+2d0*d223)
c B2(2)
      d1(0)=D0_ir(mt,mt**2,mh**2,s,mt**2,s+mt**2-2d0*p2pq,
     $     2d0*mt**2+2d0*p1p2,mt,mt,mt)
      call ddfunc_ir(mt,mt**2,mt**2,mt**2+s-2d0*p2pq,-p1p2,p1pq-p1p2,
     $     mt**2-p2pq,0d0,mt,mt,mt,d1,d2,d3,c1,c2)
      d0=d1(0)
      c0=c1(0,1)
      c11=c1(1,1)
      c12=c1(2,1)
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
      fv2(1)=2d0*mt*(c0-4d0*d20+2d0*(mt**2+p1pq)*d12+
     $     (-2d0*p1p2+2d0*p2pq)*d11+
     $     (2d0*mt**2+2d0*p1pq-2d0*p2pq)*d13)
      fv2(2)=4d0*p1p2*(d0+d11+d12+d13)
      fv2(3)=4d0*(2d0*c0+c11+c12-2d0*d20+2d0*mt**2*d12+
     $     2d0*(mt**2+p1p2)*d212+2d0*(mt**2+p1p2-p1pq)*d213+
     $     (2d0*mt**2-2d0*p2pq+s)*d13)
      fv2(4)=4d0*(2d0*c0+c11+c12-2d0*p1p2*(d0+d11+2d0*d12)-
     $     2d0*(mt**2+p1p2)*d22+2d0*(-2d0*mt**2-2d0*p1p2+p1pq)*d223+
     $     2d0*(-mt**2-p1p2+p1pq)*d23+
     $     (-4d0*p1p2+2d0*p1pq-2d0*p2pq+s)*d13)
      fv2(5)=-4d0*mt*(d12+d13+2d0*d213)
      fv2(6)=4d0*mt*(-d11+2d0*(d13+d223+d23))
      do i=1,6
         fv1(i)=(cf-cg/2d0)*fv1(i)
         fv2(i)=(cf-cg/2d0)*fv2(i)
      enddo
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine formfvf(s,t1,t2,u1,u2,mt,mh,fv1,fv2)
c form factors for the s-channel UV divergent contribution
c vertex corrections
      implicit none
      integer i
      real*8 ncf,cg,cf,mt,mh
      real*8 s,t1,t2,u1,u2,s3,p3p4,p1p2,p1p3,p1p4,p2p3,p2p4,p1pq,p2pq
      real*8 c1123(0:2),c2123(0:3)
      complex*16 zero,fv1(6),fvV21(6),fvV31(6)
      complex*16 fv2(6),fvV22(6),fvV32(6)
      complex*16 c11,c12,c21,c22,c20,c23,b0s,b1s,b20
      complex*16 c0,C0_
      real*8 delta,muedr
      common/dimreg/delta,muedr
      real*8 dZ3,dZ2,dZm,dZ1F,dZg
      common/renconst/dZ3,dZ2,dZm,dZ1F,dZg
      real*8 d1,d2,d3,d4,d5,d6,d7,d8,d9,d10

      zero=dcmplx(0d0,0d0)
      do i=1,6
         fv1(i)=zero
         fv2(i)=zero
      enddo
      s3=4d0*mt**2+mh**2-s-t1-t2-u1-u2
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

      call ccfunc(mt**2,s,-p1pq,0d0,mt,mt,
     $     c0,c11,c12,c21,c22,c20,c23)
c      c0=C0_(mt**2,s,mt**2+s-2d0*p1pq,0d0,mt,mt,1)
      call bfunc(s,mt,mt,b0s,b1s,b20)

      fvV21(1)=1d0-2d0*b0s+4d0*c20+4d0*(p1pq-mt**2)*c11
     $     +2d0*(2d0*p1pq-s)*c12+4d0*p1pq*c0
      fvV21(2)=2d0*mt*(c0+c11)
      fvV21(3)=4d0*mt*(c0-c21)
      fvV21(4)=4d0*mt*(-c12-c0+c23)
      fvV21(5)=4d0*(-c12-c0-c11-c23)
      fvV21(6)=4d0*(c12+c22)

      call ccfunc(mt**2,s,-p1pq,mt,0d0,0d0,
     $     c0,c11,c12,c21,c22,c20,c23)
c      c0=C0_(mt**2,s,mt**2+s-2d0*p1pq,mt,0d0,0d0,1)
      call bfunc(s,0d0,0d0,b0s,b1s,b20)

      fvV31(1)=1d0-2d0*b0s-4d0*c20+2d0*(2d0*p1pq-mt**2)*c11
     $     +(2d0*p1pq-s)*c12-2d0*mt**2*c0
      fvV31(2)=3d0*mt*(c0+c11)
      fvV31(3)=mt*(6d0*c0+4d0*c21+10d0*c11)
      fvV31(4)=-2d0*mt*(2d0*c11+4d0*c12+3d0*c0+2d0*c23)
      fvV31(5)=2d0*(c12-c11+2d0*c23)
      fvV31(6)=-2d0*(c12+2d0*c22)

      call ccfunc(mt**2,s,-p2pq,0d0,mt,mt,
     $     c0,c11,c12,c21,c22,c20,c23)
c      c0=C0_(mt**2,s,mt**2+s-2d0*p2pq,0d0,mt,mt,1)
      call bfunc(s,mt,mt,b0s,b1s,b20)

      fvV22(1)=1d0-2d0*b0s+4d0*c20+4d0*(p2pq-mt**2)*c11
     $     +2d0*(2d0*p2pq-s)*c12+4d0*p2pq*c0
      fvV22(2)=2d0*mt*(c0+c11)
      fvV22(3)=4d0*mt*(-c0+c21)
      fvV22(4)=4d0*mt*(c12-c11-c23)
      fvV22(5)=4d0*(-c12-c0-c11-c23)
      fvV22(6)=4d0*(c12+c22)

      call ccfunc(mt**2,s,-p2pq,mt,0d0,0d0,
     $     c0,c11,c12,c21,c22,c20,c23)
c      c0=C0_(mt**2,s,mt**2+s-2d0*p2pq,mt,0d0,0d0,1)
      call bfunc(s,0d0,0d0,b0s,b1s,b20)

      fvV32(1)=1d0-2d0*b0s-4d0*c20+2d0*(2d0*p2pq-mt**2)*c11
     $     +(2d0*p2pq-s)*c12-2d0*mt**2*c0
      fvV32(2)=3d0*mt*(c0+c11)
      fvV32(3)=mt*(-6d0*c0-4d0*c21-10d0*c11)
      fvV32(4)=2d0*mt*(-c11+4d0*c12+2d0*c23)
      fvV32(5)=2d0*(c12-c11+2d0*c23)
      fvV32(6)=-2d0*(c12+2d0*c22)

      do i=1,6
         fv1(i)=-(cf-cg/2d0)*fvV21(i)-cg/2d0*fvV31(i)
         fv2(i)=-(cf-cg/2d0)*fvV22(i)-cg/2d0*fvV32(i)
      enddo
c add counterterm to fv1(1),fv2(1):
      fv1(1)=fv1(1)+dZ1F
      fv2(1)=fv2(1)+dZ1F
c check leftover mu-dep.:
c      fv1(1)=fv1(1)-cg*dlog(muedr**2/mt**2)
c      fv2(1)=fv2(1)-cg*dlog(muedr**2/mt**2)
c add CT for comparison with Laura:
c      fv1(1)=fv1(1)-4d0*cf+cg*dlog(mt**2/muedr**2)
c      fv2(1)=fv2(1)-4d0*cf+cg*dlog(mt**2/muedr**2)
c      fv1(1)=fv1(1)+(cf+cg)*dlog(mt**2/muedr**2)
c      fv2(1)=fv2(1)+(cf+cg)*dlog(mt**2/muedr**2)
c coefficients for comparison with Laura:
c V2(1)
c      d1 =  -6.293045578405626E-003
c      d2 =   8.078588451298933E-003
c      d3 =   -13.3903724723193     
c      d4 =   1.761213070774138E-006
c      d5 =  -3.252345094524684E-005
c      d6 =   2.653852826612369E-005
c      d7 =  -8.753819152900792E-007
c      d8 =   6.550349616522915E-003
c      d9 =  -6.550349616522915E-003
c      d10 =  -2.953098960568526E-005
c      write(6,*)fv1(1),
c     $     d3+mt*(d8+d9)+(mt**2-2d0*p1pq)*d10-4d0*cf
c      write(6,*)fv1(2),
c     $     d9-mt*d10
c      write(6,*)fv1(3),
c     $     d1+d2+mt*(d4+d5+d6+d7)
c      write(6,*)fv1(4),
c     $     d2+mt*(d5+d7)
c      write(6,*)fv1(5),
c     $     d6+d7
c      write(6,*)fv1(6),
c     $     d7
c V2(2)
c      d1 =  -2.809858992913099E-003
c      d2 =   5.148024225744337E-003
c      d3 =   -10.4641773459415     
c      d4 =  -6.234146418132681E-007
c      d5 =  -2.468211219457296E-005
c      d6 =   1.837721464074102E-005
c      d7 =   2.253766215857979E-006
c      d8 =   5.159381859819261E-003
c      d9 =  -5.159381859819261E-003
c      d10 =  -2.152966341765698E-005
c      write(6,*)fv2(1),
c     $     d3+mt*(d8+d9)+(mt**2-2d0*p2pq)*d10-4d0*cf
c      write(6,*)fv2(2),
c     $     -d8-mt*d10
c      write(6,*)fv2(3),
c     $     d1+d2+mt*(d4+d5+d6+d7)
c      write(6,*)fv2(4),
c     $     -d1-mt*(d4+d6)
c      write(6,*)fv2(5),
c     $     -d4-d5+2d0*d10
c      write(6,*)fv2(6),
c     $     d4
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine formh(s,t1,t2,u1,u2,mt,mh,fh1,fh2,fh3)
      implicit none
      integer i
      real*8 ncf,cg,cf,mt,mh
      real*8 s,t1,t2,u1,u2,s3,p3p4,p1p2,p1p3,p1p4,p2p3,p2p4,p1pq,p2pq
      real*8 c1123(0:2),c2123(0:3)
      complex*16 zero,fh1(3),fh2(3),fh3(3)
      complex*16 c11,c12,c21,c22,c20,c23,b0s,b1s,b20
      complex*16 c0,C0_
      real*8 delta,muedr
      common/dimreg/delta,muedr
      real*8 dZ3,dZ2,dZm,dZ1F,dZg
      common/renconst/dZ3,dZ2,dZm,dZ1F,dZg

      zero=dcmplx(0d0,0d0)
      do i=1,3
         fh1(i)=zero
         fh2(i)=zero
         fh3(i)=zero
      enddo
      s3=4d0*mt**2+mh**2-s-t1-t2-u1-u2
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

      call ccfunc(s+mt**2-2d0*p1pq,mh**2,
     $     p2pq-p1p2-mt**2-s+2d0*p1pq,0d0,mt,mt,
     $     c0,c11,c12,c21,c22,c20,c23)
      call bfunc(mh**2,mt,mt,b0s,b1s,b20)

      fh1(1)=4d0*b0s-2d0+(4d0*mt**2+4d0*(p2pq-p1p2))*c0
     $     +4d0*(s+mt**2-2d0*p1pq+p2pq-p1p2)*c11
     $     +4d0*(-s+2d0*p1pq)*c12
      fh1(2)=2d0*mt*(-c0-2d0*c11+2d0*c12)
      fh1(3)=2d0*mt*(-c0-2d0*c12)

      call ccfunc(mt**2,mh**2,p1pq-mt**2-p1p2,0d0,mt,mt,
     $     c0,c11,c12,c21,c22,c20,c23)

      fh2(1)=4d0*b0s-2d0+(4d0*mt**2+4d0*(p1pq-p1p2))*c0
     $     +4d0*(mt**2+p1pq-p1p2)*c11+4d0*(s-2d0*p2pq)*c12
      fh2(2)=2d0*mt*(-c0-2d0*c11+2d0*c12)
      fh2(3)=2d0*mt*(-c0-2d0*c12)

      call ccfunc(mt**2-2d0*p1p3,mh**2,
     $     2d0*p1p3-mt**2+p2p3-p3p4-p1p2+p1p4,0d0,mt,mt,
     $     c0,c11,c12,c21,c22,c20,c23)

      fh3(1)=4d0*b0s-2d0+(4d0*mt**2+4d0*(p2p3-p3p4-p1p2+p1p4))*c0
     $     +4d0*(mt**2-2d0*p1p3+p2p3-p3p4-p1p2+p1p4)*c11
     $     +4d0*(2d0*p1p3-2d0*p2p4)*c12
      fh3(2)=2d0*mt*(-c0-2d0*c11+2d0*c12)
      fh3(3)=2d0*mt*(-c0-2d0*c12)

      do i=1,3
         fh1(i)=cf*fh1(i)
         fh2(i)=cf*fh2(i)
         fh3(i)=cf*fh3(i)
      enddo
c add counterterm to fh1(1),fh2(1):
      fh1(1)=fh1(1)+dZ2+dZm
      fh2(1)=fh2(1)+dZ2+dZm
      fh3(1)=fh3(1)+dZ2+dZm
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine forms(s,t1,t2,u1,u2,mt,mh,fs1,fs2,fs3,fs4)
      implicit none
      integer i
      real*8 ncf,cg,cf,mt,mh
      real*8 s,t1,t2,u1,u2,s3,p3p4,p1p2,p1p3,p1p4,p2p3,p2p4,p1pq,p2pq
      complex*16 zero,fs1(2),fs2(2),fs3(2),fs4(2),sigmav,sigmas
      complex*16 b0s,b1s,b20
      real*8 delta,muedr
      common/dimreg/delta,muedr
      real*8 dZ3,dZ2,dZm,dZ1F,dZg
      common/renconst/dZ3,dZ2,dZm,dZ1F,dZg

      zero=dcmplx(0d0,0d0)
      do i=1,2
         fs1(i)=zero
         fs2(i)=zero
         fs3(i)=zero
         fs4(i)=zero
      enddo
      s3=4d0*mt**2+mh**2-s-t1-t2-u1-u2
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

      call bfunc(s+mt**2-2d0*p1pq,mt,0d0,b0s,b1s,b20)
      sigmav=-cf*(2d0*b1s+1d0)
      sigmas=-cf*(4d0*b0s-2d0)
c add counterterm:
      sigmav=sigmav+dZ2
      sigmas=sigmas-dZ2-dZm

      fs1(1)=(s+mt**2-2d0*p1pq)*sigmav+mt**2*sigmas
      fs1(2)=mt*(sigmav+sigmas)

      call bfunc(s+mt**2-2d0*p2pq,mt,0d0,b0s,b1s,b20)
      sigmav=-cf*(2d0*b1s+1d0)
      sigmas=-cf*(4d0*b0s-2d0)
c add counterterm:
      sigmav=sigmav+dZ2
      sigmas=sigmas-dZ2-dZm

      fs2(1)=(s+mt**2-2d0*p2pq)*sigmav+mt**2*sigmas
      fs2(2)=mt*(sigmav+sigmas)

c fs3 and fs4 only contribute to the t,u-channel
      call bfunc(mt**2-2d0*p1p3,mt,0d0,b0s,b1s,b20)
      sigmav=-cf*(2d0*b1s+1d0)
      sigmas=-cf*(4d0*b0s-2d0)
c add counterterm:
      sigmav=sigmav+dZ2
      sigmas=sigmas-dZ2-dZm

      fs3(1)=(mt**2-2d0*p1p3)*sigmav+mt**2*sigmas
      fs3(2)=mt*(sigmav+sigmas)

      call bfunc(mt**2-2d0*p2p4,mt,0d0,b0s,b1s,b20)
      sigmav=-cf*(2d0*b1s+1d0)
      sigmas=-cf*(4d0*b0s-2d0)
c add counterterm:
      sigmav=sigmav+dZ2
      sigmas=sigmas-dZ2-dZm

      fs4(1)=(mt**2-2d0*p2p4)*sigmav+mt**2*sigmas
      fs4(2)=mt*(sigmav+sigmas)

      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine ggself(s,mt,mh,sigma)
      implicit none
      real*8 ncf,cg,cf,nlf,mt,mh,s
      complex*16 zero,sigma
      complex*16 b0s,b1s,b20
      real*8 delta,muedr
      common/dimreg/delta,muedr
      real*8 dZ3,dZ2,dZm,dZ1F,dZg
      common/renconst/dZ3,dZ2,dZm,dZ1F,dZg

      zero=dcmplx(0d0,0d0)
      sigma=zero
c color factors:
      ncf=3d0
      cg=ncf
      cf=(ncf**2-1d0)/ncf/2d0
c number of light flavors:
      nlf=5d0

      call bfunc(s,mt,mt,b0s,b1s,b20)
      sigma=(2d0/3d0*nlf*(-dlog(s/muedr**2)+5d0/3d0)+
     $     2d0/3d0*(2d0*mt**2/s*dlog(mt**2/muedr**2)+
     $     (1d0+2d0*mt**2/s)*b0s-1d0/3d0))-
     $     cg*(-5d0/3d0*dlog(s/muedr**2)+31d0/9d0)

c add counterterm:
      sigma=sigma+dZ3
c check of leftover mu-dep.:
c      sigma=sigma-(-2d0/3d0*nlf+5d0/3d0*cg)*dlog(s/muedr**2)
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine ggvert(s,mt,mh,vertex)
      implicit none
      real*8 ncf,cg,cf,nlf,mt,mh,s
      complex*16 zero,vertex
      complex*16 b0s,b1s,b20,c0,c11,c12,c21,c22,c20,c23,b00,c00
      real*8 pi,pi2,w2
      common/para/pi,pi2,w2
      real*8 delta,muedr
      common/dimreg/delta,muedr
      real*8 dZ3,dZ2,dZm,dZ1F,dZg
      common/renconst/dZ3,dZ2,dZm,dZ1F,dZg

      zero=dcmplx(0d0,0d0)
      vertex=zero
c color factors:
      ncf=3d0
      cg=ncf
      cf=(ncf**2-1d0)/ncf/2d0
c number of light flavors:
      nlf=5d0

      call ccfunc(0d0,s,-s/2d0,mt,mt,mt,
     $     c0,c11,c12,c21,c22,c20,c23)
      call bfunc(s,mt,mt,b0s,b1s,b20)
c b00 and C00 are IR divergent:
c b00 is also UV divergent:
      b00=dcmplx(-dlog(mt**2/muedr**2),0d0)
      c00=1d0/s*dcmplx(dlog(mt**2/s)**2/2d0-2d0/3d0*pi**2,0d0)

c G[1+eps] (4 pi mu^2/mt^2)^eps x
c 3/4 (-1/eps_IR^2-1/eps_IR ln(mt2/s)-2/eps_IR) has been subtracted:
      vertex=2d0/3d0*nlf*(-dlog(s/muedr**2)+7d0/6d0)+
     $     (5d0/3d0*b0s-4d0*c20+4d0/9d0-4d0*mt**2*c12+
     $     4d0/3d0*mt**2/s*(b0s+dlog(mt**2/muedr**2)))+
     $     cg/4d0*((51d0/6d0*dlog(s/muedr**2)-95d0/6d0-
     $     3d0*s*c00+15d0*b00)+(dlog(s/muedr**2)/6d0-11d0/18d0)-
     $     9d0*b00)

c add contribution from the insertions into the external 
c gluon legs: 
c G[1+eps] (4 pi mu^2/mt^2)^eps 1/eps_IR (2/3 nlf-5/3 cg)
c has been subtracted
      vertex=vertex+(-2d0/3d0*nlf+5d0/3d0*cg)*b00

c add counterterm:
      vertex=vertex+dZg+3d0/2d0*dZ3
c check of leftover mu-dep.:
c      vertex=vertex+cg*dlog(s/muedr**2)
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine renorm(mt,mh)
      implicit none
      real*8 mt,mh,dZ3,dZ2,dZm,dZ1F,dZg
      real*8 db0,ncf,cg,cf
      common/renconst/dZ3,dZ2,dZm,dZ1F,dZg
      real*8 delta,muedr
      common/dimreg/delta,muedr
c color factors:
      ncf=3d0
      cg=ncf
      cf=(ncf**2-1d0)/ncf/2d0
c db0 is IR divergent:
c -G[1+eps] (4 pi mu^2/mt^2)^eps 1/2/eps_IR has been subtracted
      db0=-1d0/mt**2
c dZ3 in modified MSbar:
c UV divergences have been subtracted: N_eps=2/eps-gamma_E+ln(4pi)
      dZ3=2d0/3d0*dlog(mt**2/muedr**2)
      dZ2=cf*(dlog(mt**2/muedr**2)+4d0*mt**2*db0)
      dZm=-cf*(-3d0*dlog(mt**2/muedr**2)+4d0)
      dZ1F=dZ2
      dZg=-1d0/3d0*dlog(mt**2/muedr**2)
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine gluform_t(s,t1,t2,u1,u2,mt,mh,kmv1,kmv2,kmv3,
     $     kmv4,kmv5)
      implicit none
      integer i,j
      real*8 s,s3,p3p4,p1p2,p1p3,p1p4,p2p3,p2p4
      real*8 t1,t2,u1,u2,mt,mh,s1,s2
      complex*16 zero,kmv1(68),kmv2(68),kmv3(68)
      complex*16 kmvf1(5,68),kmvf2(5,68),kmvf3(5,68)
      complex*16 kmvf4(6,68),kmv4(68),kmvf5(2,68),kmv5(68)
      complex*16 fvf1(6),fvf2(6),fh1(3),fh2(3),fh3(3),
     $     fs1(2),fs2(2),fs3(2),fs4(2),kmvt(15,68)
      complex*16 fv3,fv4,fv8,fv10
      
      zero=dcmplx(0d0,0d0)
      do i=1,68
         kmv1(i)=zero
         kmv2(i)=zero
         kmv3(i)=zero
         kmv4(i)=zero
         kmv5(i)=zero
         do j=1,5
            kmvf1(j,i)=zero
            kmvf2(j,i)=zero
            kmvf3(j,i)=zero
         enddo
         do j=1,6
            kmvf4(j,i)=zero
         enddo
         do j=1,2
            kmvf5(j,i)=zero
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

c contribution from V7(1), V7(2), V7(3):

      call formfvf_t(s,t1,t2,u1,u2,mt,mh,fvf1,fvf2)
c
c non abelian:
c      write(6,*)dreal(fvf1(3)),-(3.963644095016382E-003
c     $     -8.347535652522810E-003)
c      write(6,*)dreal(fvf1(1)),3.36820473431370
c     $     +(mt**2+2d0*p3p4-2d0*p1p4-2d0*p1p3)
c     $     *4.427137414317813E-005
c      write(6,*)dreal(fvf1(2)),9.895164879666209E-003
c      write(6,*)dreal(fvf1(4)),-(2.875522775256243E-006
c     $     +5.361881469548196E-005-3.492393359087430E-005
c     $     -1.373362874139274E-005)
c      write(6,*)dreal(fvf1(5)),-(-3.492393359087430E-005
c     $     -1.373362874139274E-005)
c      write(6,*)dreal(fvf1(6)),-4.427137414317813E-005
c      write(6,*)'sum:', 
c     $     -2.195204695664224E-007*(3.963644095016382E-003+
c     $     -8.347535652522810E-003)+7.885922844466861E-009*
c     $     3.36820473431370-2.198200922432374E-003*(
c     $     2.875522775256243E-006+5.361881469548196E-005)+
c     $     1.614654495887201E-005*(-3.492393359087430E-005
c     $     -1.373362874139274E-005)-(-3.415844932398249E-007
c     $     +3.272555046983402E-006)*9.895164879666209E-003
c     $     + 8.828700346744467E-004*4.427137414317813E-005
c abelian: 
c      write(6,*)dreal(fvf1(3)),-(5.161698927677813E-004
c     $     + 1.332955231590983E-004)
c      write(6,*)dreal(fvf1(1)),0.362057014051139 
c     $     +(mt**2+2d0*p3p4-2d0*p1p4-2d0*p1p3)
c     $     *9.331399654121832E-007
c      write(6,*)dreal(fvf1(2)),0d0
c      write(6,*)dreal(fvf1(4)),-(1.140369175162632E-006
c     $     -6.903663086612742E-007-2.556646239485641E-006
c     $     +1.235890579820817E-006)
c      write(6,*)dreal(fvf1(5)),-(-2.556646239485641E-006
c     $     +1.235890579820817E-006)
c      write(6,*)dreal(fvf1(6)),- 9.331399654121832E-007

c V7(1)
      kmvf1(1,5)=2d0*p1p3*(fvf1(1)-2d0*(p1p4-p3p4)*fvf1(6))
      kmvf1(1,6)=-fvf1(1)+2d0*(mt*fvf1(2)+(mt**2-p1p3)*fvf1(6))
      kmvf1(1,7)=-2d0*p1p3*(fvf1(2)+mt*fvf1(6))
      kmvf1(1,8)=2d0*(p1p4-p3p4)*(fvf1(2)+mt*fvf1(6))
      kmvf1(1,14)=-2d0*fvf1(2)+fvf1(3)+mt*fvf1(4)-2d0*mt*fvf1(6)
      kmvf1(1,17)=2d0*(fvf1(1)-2d0*mt*fvf1(2)-2d0*(mt**2-p1p3)*fvf1(6))
      kmvf1(1,20)=2d0*p1p3*(-fvf1(4)+fvf1(5)+2d0*fvf1(6))
      kmvf1(1,26)=2d0*(-fvf1(1)+mt*fvf1(3)+(mt**2-p1p3)*fvf1(4)
     $     -(p1p4-p3p4)*fvf1(5))
      kmvf1(1,29)=-4d0*(p1p4-p3p4)*(fvf1(2)+mt*fvf1(6))
      kmvf1(1,32)=-2d0*p1p3*(fvf1(3)+mt*fvf1(4))
      kmvf1(1,50)=2d0*(2d0*fvf1(2)-fvf1(3)-mt*fvf1(4)+2d0*mt*fvf1(6))
      kmvf1(1,68)=4d0*(fvf1(1)-mt*fvf1(3)+(p1p3-mt**2)*fvf1(4)
     $     +(p1p4-p3p4)*fvf1(5))

c contribution for comparison with Laura:
c same as above for fv3=0,fv4=0,fv8=-fvf1(4),fv10=-fvf1(6)      
c      fvf1(1)=dcmplx(0d0,0d0)
c      fvf1(2)=dcmplx(1d0,0d0)
c      fvf1(3)=dcmplx(0d0,0d0)
c      fvf1(4)=dcmplx(0d0,0d0)
c      fvf1(5)=dcmplx(0d0,0d0)
c      fvf1(6)=dcmplx(0d0,0d0)
c      fv3=dcmplx(0d0,0d0)
c      fv4=dcmplx(0d0,0d0)
c      fv8=dcmplx(0d0,0d0)
c      fv10=dcmplx(0d0,0d0)
c      fv8=-fvf1(4)
c      fv10=-fvf1(6)
c      kmvf1(1,5)=2d0*p1p3*(fvf1(1)+2d0*(p1p4-p3p4)*fv10+mt*fv4)
c      kmvf1(1,6)=-fvf1(1)+2d0*(mt*fvf1(2)-(mt**2-p1p3)*fv10)+mt*fv4
c      kmvf1(1,7)=-2d0*p1p3*(fvf1(2)-mt*fv10+fv4)
c      kmvf1(1,8)=2d0*(p1p4-p3p4)*(fvf1(2)-mt*fv10)-2d0*p1p3*fv4
c      kmvf1(1,14)=-2d0*fvf1(2)+fvf1(3)-mt*fv8+2d0*mt*fv10
c      kmvf1(1,17)=2d0*(fvf1(1)-2d0*mt*fvf1(2)+2d0*p1p3*fvf1(6)+
c     $     2d0*mt**2*fv10-mt*fv4)
c      kmvf1(1,20)=2d0*p1p3*(fv8+fvf1(5)-2d0*fv10)
c      kmvf1(1,26)=2d0*(-fvf1(1)+mt*fvf1(3)-(mt**2-p1p3)*fv8+mt*fv4
c     $     -(p1p4-p3p4)*fvf1(5))
c      kmvf1(1,29)=-4d0*(p1p4-p3p4)*(fvf1(2)-mt*fv10)-4d0*p1p3*fv3
c      kmvf1(1,32)=-2d0*p1p3*(fvf1(3)-mt*fv8+2d0*fv4)
c      kmvf1(1,41)=2d0*(2d0*(fvf1(6)+fv10)-fvf1(4)-fv8)
c      kmvf1(1,50)=2d0*(2d0*fvf1(2)-fvf1(3)+mt*fv8-2d0*mt*fv10)
c      kmvf1(1,68)=4d0*(fvf1(1)-mt*fvf1(3)+p1p3*fvf1(4)-mt*fv4+mt**2*fv8
c     $     +(p1p4-p3p4)*fvf1(5))
c      kmvf1(1,59)=4d0*(fv3+fv4)
c      kmvf1(1,11)=2d0*(fv3+fv4)
c      kmvf1(1,23)=-4d0*mt*(fv3+fv4)+4d0*(p3p4-p1p4)*(fv10+fvf1(6))

c V7(3)
      kmvf2(1,5)=2d0*p2p4*fvf2(1)
      kmvf2(1,6)=-fvf2(1)+2d0*mt*fvf2(2)
      kmvf2(1,7)=2d0*p2p4*fvf2(2)
      kmvf2(1,8)=2d0*p2p4*fvf2(2)
      kmvf2(1,13)=-2d0*fvf2(2)-fvf2(3)
      kmvf2(1,16)=-2d0*fvf2(1)+4d0*mt*fvf2(2)
      kmvf2(1,19)=2d0*p2p4*fvf2(4)
      kmvf2(1,25)=2d0*(fvf2(1)+mt*fvf2(3)+p2p4*fvf2(4))
      kmvf2(1,28)=4d0*p2p4*fvf2(2)
      kmvf2(1,31)=2d0*p2p4*fvf2(3)
      kmvf2(1,46)=-4d0*fvf2(2)-2d0*fvf2(3)
      kmvf2(1,64)=4d0*(fvf2(1)+mt*fvf2(3)+p2p4*fvf2(4))
c V7(2)
      kmvf3(1,6)=-fvf2(1)+2d0*mt*fvf2(2)
      kmvf3(1,8)=2d0*p2p4*fvf2(2)
      kmvf3(1,13)=-2d0*fvf2(2)-fvf2(3)
      kmvf3(1,17)=2d0*fvf2(1)-4d0*mt*fvf2(2)
      kmvf3(1,25)=2d0*(fvf2(1)+mt*fvf2(3)+p2p4*fvf2(4))
      kmvf3(1,29)=-4d0*p2p4*fvf2(2)
      kmvf3(1,47)=4d0*fvf2(2)+2d0*fvf2(3)
      kmvf3(1,65)=-4d0*(fvf2(1)+mt*fvf2(3)+p2p4*fvf2(4))

c contribution from V8(1), V8(2), V8(3):
      call formfvf_t(s,t2,t1,u2,u1,mt,mh,fvf1,fvf2)
      fvf1(2)=-fvf1(2)
      fvf1(3)=-fvf1(3)
      fvf2(2)=-fvf2(2)
      fvf2(3)=-fvf2(3)
c V8(3)
      kmvf2(2,5)=2d0*p2p4*(fvf1(1)-2d0*(p2p3-p3p4)*fvf1(6))
      kmvf2(2,6)=-fvf1(1)+2d0*(-mt*fvf1(2)+(mt**2-p2p4)*fvf1(6))
      kmvf2(2,7)=2d0*(p2p3-p3p4)*(fvf1(2)-mt*fvf1(6))
      kmvf2(2,8)=-2d0*p2p4*(fvf1(2)-mt*fvf1(6))
      kmvf2(2,10)=-2d0*fvf1(2)+fvf1(3)-mt*fvf1(4)+2d0*mt*fvf1(6)
      kmvf2(2,16)=2d0*(-fvf1(1)-mt*fvf1(3)+(mt**2-p2p4)*fvf1(4)
     $     +(p3p4-p2p3)*fvf1(5))
      kmvf2(2,22)=2d0*p2p4*(-fvf1(4)+fvf1(5)+2d0*fvf1(6))
      kmvf2(2,25)=2d0*(fvf1(1)+2d0*mt*fvf1(2)+2d0*(p2p4-mt**2)*fvf1(6))
      kmvf2(2,28)=2d0*p2p4*(-fvf1(3)+mt*fvf1(4))
      kmvf2(2,31)=4d0*(p2p3-p3p4)*(-fvf1(2)+mt*fvf1(6))
      kmvf2(2,55)=2d0*(2d0*fvf1(2)-fvf1(3)+mt*fvf1(4)-2d0*mt*fvf1(6))
      kmvf2(2,64)=4d0*(fvf1(1)+mt*fvf1(3)-(mt**2-p2p4)*fvf1(4)
     $     -(p3p4-p2p3)*fvf1(5))
c V8(1)
      kmvf1(2,5)=2d0*p1p3*fvf2(1)
      kmvf1(2,6)=-fvf2(1)-2d0*mt*fvf2(2)
      kmvf1(2,7)=2d0*p1p3*fvf2(2)
      kmvf1(2,8)=2d0*p1p3*fvf2(2)
      kmvf1(2,11)=-2d0*fvf2(2)-fvf2(3)
      kmvf1(2,17)=2d0*(fvf2(1)-mt*fvf2(3)+p1p3*fvf2(4))
      kmvf1(2,23)=2d0*p1p3*fvf2(4)
      kmvf1(2,26)=-2d0*fvf2(1)-4d0*mt*fvf2(2)
      kmvf1(2,29)=2d0*p1p3*fvf2(3)
      kmvf1(2,32)=4d0*p1p3*fvf2(2)
      kmvf1(2,59)=-4d0*fvf2(2)-2d0*fvf2(3)
      kmvf1(2,68)=4d0*(fvf2(1)-mt*fvf2(3)+p1p3*fvf2(4))
c V8(2)
      kmvf3(2,6)=-fvf2(1)-2d0*mt*fvf2(2)
      kmvf3(2,7)=2d0*p1p3*fvf2(2)
      kmvf3(2,11)=-2d0*fvf2(2)-fvf2(3)
      kmvf3(2,17)=2d0*(fvf2(1)-mt*fvf2(3)+p1p3*fvf2(4))
      kmvf3(2,25)=2d0*fvf2(1)+4d0*mt*fvf2(2)
      kmvf3(2,31)=-4d0*p1p3*fvf2(2)
      kmvf3(2,56)=4d0*fvf2(2)+2d0*fvf2(3)
      kmvf3(2,65)=-4d0*(fvf2(1)-mt*fvf2(3)+p1p3*fvf2(4))

c contribution from V9(1), V9(2), V9(3):
      call formh(s,t1,t2,u1,u2,mt,mh,fh1,fh2,fh3)
c V9(1)
      kmvf1(3,5)=2d0*p1p3*(fh1(1)+mt*fh1(2)+mt*fh1(3))
      kmvf1(3,6)=-fh1(1)-mt*fh1(2)-mt*fh1(3)
      kmvf1(3,8)=(s-2d0*(p1p3+p1p4))*fh1(2)
      kmvf1(3,17)=2d0*(fh1(1)+mt*fh1(3)+mt*fh1(2))
      kmvf1(3,26)=-2d0*(fh1(1)+mt*fh1(3)+mt*fh1(2))
      kmvf1(3,29)=-2d0*(s-2d0*(p1p3+p1p4))*fh1(2)
      kmvf1(3,68)=4d0*(fh1(1)+mt*fh1(3)+mt*fh1(2))
c V9(3)
      kmvf2(3,5)=2d0*p2p4*(fh2(1)+mt*fh2(2)+mt*fh2(3))
      kmvf2(3,6)=-fh2(1)-mt*fh2(2)-mt*fh2(3)
      kmvf2(3,7)=-(s-2d0*(p2p4+p2p3))*fh2(3)
      kmvf2(3,16)=-2d0*(fh2(1)+mt*fh2(3)+mt*fh2(2))
      kmvf2(3,25)=2d0*(fh2(1)+mt*fh2(3)+mt*fh2(2))
      kmvf2(3,31)=2d0*(s-2d0*(p2p4+p2p3))*fh2(3)
      kmvf2(3,64)=4d0*(fh2(1)+mt*fh2(3)+mt*fh2(2))
c V9(2)
      kmvf3(3,6)=-(fh3(1)+mt*fh3(2)+mt*fh3(3))
      kmvf3(3,7)=2d0*p1p3*fh3(2)
      kmvf3(3,8)=-2d0*p2p4*fh3(3)
      kmvf3(3,17)=2d0*(fh3(1)+mt*fh3(3)+mt*fh3(2))
      kmvf3(3,25)=2d0*(fh3(1)+mt*fh3(3)+mt*fh3(2))
      kmvf3(3,29)=4d0*p2p4*fh3(3)
      kmvf3(3,31)=-4d0*p1p3*fh3(2)
      kmvf3(3,65)=-4d0*(fh3(1)+mt*fh3(3)+mt*fh3(2))

c contribution from S3(1), S3(2) and S4(2)
      call forms(s,t1,t2,u1,u2,mt,mh,fs1,fs2,fs3,fs4)
c S3(1)      
      kmvf1(4,5)=2d0*p1p3*(fs1(1)+mt*fs1(2))
      kmvf1(4,6)=-(fs1(1)+mt*fs1(2))
      kmvf1(4,8)=fs1(2)*(s-2d0*(p1p3+p1p4))
      kmvf1(4,17)=2d0*(fs1(1)+mt*fs1(2))
      kmvf1(4,26)=-2d0*(fs1(1)+mt*fs1(2))
      kmvf1(4,29)=-2d0*fs1(2)*(s-2d0*(p1p3+p1p4))
      kmvf1(4,68)=4d0*(fs1(1)+mt*fs1(2))
c S3(2)
      kmvf2(4,5)=2d0*p2p4*(fs2(1)+mt*fs2(2))
      kmvf2(4,6)=-(fs2(1)+mt*fs2(2))
      kmvf2(4,7)=-fs2(2)*(s-2d0*(p2p3+p2p4))
      kmvf2(4,16)=-2d0*(fs2(1)+mt*fs2(2))
      kmvf2(4,25)=2d0*(fs2(1)+mt*fs2(2))
      kmvf2(4,31)=2d0*fs2(2)*(s-2d0*(p2p3+p2p4))
      kmvf2(4,64)=4d0*(fs2(1)+mt*fs2(2))
c S4(2)
      kmvf3(4,6)=-(fs3(1)+mt*fs3(2))/2d0/p1p3
      kmvf3(4,7)=fs3(2)
      kmvf3(4,17)=(fs3(1)+mt*fs3(2))/p1p3
      kmvf3(4,25)=(fs3(1)+mt*fs3(2))/p1p3
      kmvf3(4,31)=-2d0*fs3(2)
      kmvf3(4,65)=-2d0*(fs3(1)+mt*fs3(2))/p1p3
      do i=1,68
         kmvf1(4,i)=-1d0/(s+t2+u2-2d0*mt**2)*kmvf1(4,i)
         kmvf2(4,i)=-1d0/(s+t1+u1-2d0*mt**2)*kmvf2(4,i)
      enddo
c contribution from S4(1), S4(3), S4(4)
c S4(1) 
      kmvf1(5,5)=(fs3(1)-mt*fs3(2))
      kmvf1(5,6)=-(fs3(1)+mt*fs3(2))/2d0/p1p3
      kmvf1(5,7)=fs3(2)
      kmvf1(5,8)=fs3(2)
      kmvf1(5,17)=(fs3(1)+mt*fs3(2))/p1p3
      kmvf1(5,26)=-(fs3(1)+mt*fs3(2))/p1p3
      kmvf1(5,29)=-2d0*fs3(2)
      kmvf1(5,32)=2d0*fs3(2)
      kmvf1(5,68)=2d0*(fs3(1)+mt*fs3(2))/p1p3
c S4(4)
      kmvf2(5,5)=(fs4(1)-mt*fs4(2))
      kmvf2(5,6)=-(fs4(1)+mt*fs4(2))/2d0/p2p4
      kmvf2(5,7)=-fs4(2)
      kmvf2(5,8)=-fs4(2)
      kmvf2(5,16)=-(fs4(1)+mt*fs4(2))/p2p4
      kmvf2(5,25)=(fs4(1)+mt*fs4(2))/p2p4
      kmvf2(5,28)=-2d0*fs4(2)
      kmvf2(5,31)=2d0*fs4(2)
      kmvf2(5,64)=2d0*(fs4(1)+mt*fs4(2))/p2p4
c S4(3)
      kmvf3(5,6)=-(fs4(1)+mt*fs4(2))/2d0/p2p4
      kmvf3(5,8)=-fs4(2)
      kmvf3(5,17)=(fs4(1)+mt*fs4(2))/p2p4
      kmvf3(5,25)=(fs4(1)+mt*fs4(2))/p2p4
      kmvf3(5,29)=2d0*fs4(2)
      kmvf3(5,65)=-2d0*(fs4(1)+mt*fs4(2))/p2p4
c IR+UV finite contributions:
      call kmvt_finite(s,t1,t2,u1,u2,mt,mh,kmvt)
c contribution from V10(1+2),V10(3+4),B5(1),B5(2),B6(1),B6(2)
      do i=1,68
         kmvf4(1,i)=2d0*kmvt(1,i)/(t2-mt**2)/
     $        (2d0*mt**2+2d0*p1p2-2d0*p1p3-2d0*p2p3)
         kmvf4(2,i)=2d0*kmvt(2,i)/(t1-mt**2)/
     $        (2d0*mt**2+2d0*p1p2-2d0*p1p4-2d0*p2p4)
         kmvf4(3,i)=kmvt(3,i)/(t2-mt**2)
         kmvf4(4,i)=kmvt(4,i)/(t2-mt**2)
         kmvf4(5,i)=kmvt(5,i)/(t1-mt**2)
         kmvf4(6,i)=kmvt(6,i)/(t1-mt**2)
      enddo
c contribution from B9(1),B9(2):
      do i=1,68
         kmvf5(1,i)=kmvt(7,i)/(s1-mt**2)
         kmvf5(2,i)=kmvt(8,i)/(s2-mt**2)
      enddo
c sum all contributions:
      do i=1,68
         do j=1,5
            kmv1(i)=kmv1(i)+kmvf1(j,i)
            kmv2(i)=kmv2(i)+kmvf2(j,i)
            kmv3(i)=kmv3(i)+kmvf3(j,i)
         enddo
         do j=1,6
            kmv4(i)=kmv4(i)+kmvf4(j,i)
         enddo
         do j=1,2
            kmv5(i)=kmv5(i)+kmvf5(j,i)
         enddo
      enddo
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine formfvf_t(s,t1,t2,u1,u2,mt,mh,fv1,fv2)
c form factors for the t-channel UV divergent contribution
c vertex corrections, V7,V8
      implicit none
      integer i
      real*8 ncf,cg,cf,mt,mh,nlf
      real*8 s,t1,t2,u1,u2,s3,p3p4,p1p2,p1p3,p1p4,p2p3,p2p4,p1pq,p2pq
      real*8 c1123(0:2),c2123(0:3)
      complex*16 zero,fv1(6),fvV21(6),fvV31(6)
      complex*16 fv2(6),fvV22(6),fvV32(6)
      complex*16 c11,c12,c21,c22,c20,c23,b0s,b1s,b20
      complex*16 c0,C0_,C0_ir,b00
      complex*16 spence,ieps
      real*8 delta,muedr
      common/dimreg/delta,muedr
      real*8 dZ3,dZ2,dZm,dZ1F,dZg
      common/renconst/dZ3,dZ2,dZm,dZ1F,dZg
      real*8 d1,d2,d3,d4,d5,d6,d7,d8,d9,d10
      real*8 tau,omega

      zero=dcmplx(0d0,0d0)
      ieps=dcmplx(0d0,1d-8)
      do i=1,6
         fv1(i)=zero
         fv2(i)=zero
      enddo
      s3=4d0*mt**2+mh**2-s-t1-t2-u1-u2
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
c V7(1)
      call ccfunc(mt**2-2d0*p1pq+s,0d0,p1p4-p3p4,0d0,mt,mt,
     $     c0,c11,c12,c21,c22,c20,c23)
c      c0=C0_(mt**2,s,mt**2+s-2d0*p1pq,0d0,mt,mt,1)
      call bfunc(0d0,mt,mt,b0s,b1s,b20)

      fvV21(1)=-1d0+2d0*b0s-4d0*c20
     $     +4d0*(mt**2-2d0*(p1p3+p1p4-p3p4))*c11
     $     +4d0*(p1p4-p3p4)*c12+4d0*(p3p4-p1p3-p1p4)*c0
      fvV21(2)=zero
      fvV21(3)=-8d0*mt*(c0+c11)
      fvV21(4)=4d0*(c0+2d0*c11+c21)
      fvV21(5)=4d0*(c0+2d0*c11-c12+c21-c23)
      fvV21(6)=-2d0*(c0+c11)

      c0=C0_ir(mt,0d0,s+mt**2-2d0*p1pq,
     $     mt**2-2d0*p1p3,0d0,0d0,mt)
      call ccfunc_ir(mt,0d0,mt**2-2d0*p1p3,p1p4-p3p4,0d0,0d0,mt,
     $     c0,c11,c12,c21,c22,c20,c23)
c compare with Laura:
c      tau=2d0*p1p3
c      omega=s-2d0*p1pq
c      write(6,*)c0,1d0/(tau+omega)*(-dlog(tau/mt**2)**2/2d0+
c     $     dlog(omega/mt**2)**2/2d0+dlog(1d0+mt**2/omega)**2/2d0+
c     $     spence(1d0-mt**2/tau+ieps)+spence(omega/(omega+mt**2)+ieps)-
c     $     5d0/6d0*(4d0*datan(1d0))**2)
      call bfunc(mt**2-2d0*p1pq+s,0d0,mt,b0s,b1s,b20)

      fvV31(1)=1d0-2d0*b0s-4d0*c20+4d0*(p3p4-p1p4)*c0+
     $     2d0*(2d0*(p1p3-p1p4+p3p4)-mt**2)*c12+2d0*(p3p4-p1p4)*c11
      fvV31(2)=3d0*mt*c0
      fvV31(3)=-6d0*mt*c12
      fvV31(4)=4d0*(c12+c22)
      fvV31(5)=2d0*(c0-c11+c12-2d0*c23)
      fvV31(6)=-3d0*(c0+c12)
c V7(2),V7(3)
      call ccfunc(mt**2,0d0,-p2p4,0d0,mt,mt,
     $     c0,c11,c12,c21,c22,c20,c23)
c      c0=C0_(mt**2,s,mt**2+s-2d0*p2pq,0d0,mt,mt,1)
      call bfunc(0d0,mt,mt,b0s,b1s,b20)

      fvV22(1)=-1d0+2d0*b0s-4d0*c20-4d0*(p2p4-mt**2)*c11
     $     -4d0*p2p4*c12-4d0*p2p4*c0
      fvV22(2)=-2d0*mt*(c0+c11)
      fvV22(3)=-4d0*mt*(-c0+c21)
      fvV22(4)=4d0*(c0+c11+c12+c23)
      fvV22(5)=zero
      fvV22(6)=zero

      c0=C0_ir(mt,0d0,mt**2-2d0*p2p4,mt**2,0d0,0d0,mt)
      call ccfunc_ir(mt,0d0,mt**2-2d0*p2p4,-p2p4,0d0,0d0,mt,
     $     c0,c11,c12,c21,c22,c20,c23)
c compare with Laura:
c      write(6,*)c0,1d0/2d0/p2p4*(-dlog(2d0*p2p4/mt**2)**2/2d0
c     $     +spence(1d0-mt**2/2d0/(p2p4+ieps)))
      call bfunc(mt**2,0d0,mt,b0s,b1s,b20)

      fvV32(1)=1d0-2d0*b0s-4d0*c20+2d0*p2p4*c11
     $     +(2d0*p2p4-2d0*mt**2)*c12-2d0*p2p4*c0
      fvV32(2)=-3d0*mt*c12
      fvV32(3)=2d0*mt*(c12-2d0*c22)
      fvV32(4)=2d0*(2d0*c0+c11+4d0*c12+2d0*c22+2d0*c23)
      fvV32(5)=zero
      fvV32(6)=zero

      do i=1,6
         fv1(i)=(cf-cg/2d0)*fvV21(i)-cg/2d0*fvV31(i)
         fv2(i)=(cf-cg/2d0)*fvV22(i)-cg/2d0*fvV32(i)
      enddo
c add counterterm to fv1(1),fv2(1):
      fv1(1)=fv1(1)+dZ1F
      fv2(1)=fv2(1)+dZ1F
c add CT for comparison with Laura:
c      fv1(1)=fv1(1)-4d0*cf+cg*dlog(mt**2/muedr**2)
c      fv2(1)=fv2(1)-4d0*cf+cg*dlog(mt**2/muedr**2)
c      fv1(1)=fv1(1)+(cf+cg)*dlog(mt**2/muedr**2)
c      fv2(1)=fv2(1)+(cf+cg)*dlog(mt**2/muedr**2)
c add contribution from the insertions into the external 
c gluon legs: 
c G[1+eps] (4 pi mu^2/mt^2)^eps 1/eps_IR (2/3 nlf-5/3 cg)
c has been subtracted
      nlf=5d0
      b00=dcmplx(-dlog(mt**2/muedr**2),0d0)
      fv1(1)=fv1(1)+((-2d0/3d0*nlf+5d0/3d0*cg)*b00)/2d0
      fv2(1)=fv2(1)+((-2d0/3d0*nlf+5d0/3d0*cg)*b00)/2d0
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine ccfunc(xpi1,xpi2,pi12,m0,mi1,mi2,
     $     c0,c11,c12,c21,c22,c20,c23)
*---------------------------------------------------------------------
c entwicklungskoeff. fuer das vektorielle und tensorielle dreipkt.-integral
c d0=k^2-m0^2, d1=(k+pi1)^2-m1^2, d2=(k+pi1+pi2)^2-m2^2
c xpi1=pi1^2,xpi2=pi2^2,pi12=pi1.pi2
c only IR finite C-functions !
*---------------------------------------------------------------------
      implicit none
      real*8 xpi1,xpi2,pi12,mi1,mi2,m0,xmi1,xmi2,xm0,f1,f2,q12,det
c      real*8 ch1(0:2)
      complex*16 c11,c12,c21,c22,c20,c23,
     &     b0q20,b0110,b0221,b1110,b1q20,b1221,b20
      complex*16 c0,C0_
      complex*16 a,b,c,d
      xm0 = m0*m0
      xmi1 = mi1*mi1
      xmi2 = mi2*mi2
      q12=xpi1+xpi2+2d0*pi12
      f1=xmi1-xm0-xpi1
      f2=xpi1-xmi1+xmi2-q12
      det=xpi1*xpi2-pi12**2
c      if(dabs(det).lt.1d-1)then
c         write(6,*)'warning: potential instability in ccfunc'
c         stop
c      endif
*
      if(dabs(q12).ne.0d0.and.dabs(q12).lt.1d-6)q12=0d0
      if(dabs(xpi1).ne.0d0.and.dabs(xpi1).lt.1d-6)xpi1=0d0
      if(dabs(xpi2).ne.0d0.and.dabs(xpi2).lt.1d-6)xpi2=0d0
      c0=C0_(xpi1,xpi2,q12,m0,mi1,mi2,0)
c      write(6,*)c0
c      call cmue(xpi1,xpi2,q12,m0,mi1,mi2,ch1)
c      c0=ch1(0)
c      write(6,*)c0
      call bfunc(q12,m0,mi2,b0q20,b1q20,b20)
      call bfunc(xpi1,m0,mi1,b0110,b1110,b20)
      call bfunc(xpi2,mi1,mi2,b0221,b1221,b20)
*---------------------------------------------------------------------
*--entwicklungskoeff. des vektorintegrals:----------------------------
      c11=1d0/2d0/det*(pi12*(b0q20-b0110)+xpi2*(b0q20-b0221)
     $     +(xpi2*f1-pi12*f2)*c0)

      c12=-1d0/2d0/det*(pi12*(b0q20-b0221)+xpi1*(b0q20-b0110)
     $     +(-xpi1*f2+pi12*f1)*c0)
*---------------------------------------------------------------------
*------------entwicklungskoeff. des tensorintegrals--------------
      c20 = (b0221+1d0)/4d0+xm0/2d0*c0-(f1*c11+f2*c12)/4d0
     $     -pi12/4d0/det*(xpi1*b1110+xpi2*b1221-q12*b1q20
     $     -pi12*b0221+(xpi1*f2-pi12*f1)*c11
     $     +(pi12*f2-xpi2*f1)*c12)

      a=1d0/2d0*((xpi1+pi12)*b1q20-pi12*b1221+xpi1*b0221
     $     +f1*(xpi1*c11+pi12*c12))
      b=1d0/2d0*(-(xpi2+pi12)*b1q20+pi12*b1110
     $     +f2*(pi12*c11+xpi2*c12))
      c=b0221+xm0*c0
      d=1d0/2d0*((xpi2+pi12)*b1q20+pi12*b0221-xpi2*b1221
     $     +f1*(pi12*c11+xpi2*c12))

      c21 = 1d0/det*(-b+xpi2*(c-3d0*c20+1d0/2d0)) 

      c22 = 1d0/det*(-a+xpi1*(c-3d0*c20+1d0/2d0)) 

      c23 = 1d0/det*(d+pi12*(-c+3d0*c20-1d0/2d0)) 
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine ccfunc2(xpi1,xpi2,pi12,m0,mi1,mi2,
     $     c0,c11,c12,c21,c22,c20,c23)
*---------------------------------------------------------------------
c entwicklungskoeff. fuer das vektorielle und tensorielle dreipkt.-integral
c d0=k^2-m0^2, d1=(k+pi1)^2-m1^2, d2=(k+pi2)^2-m2^2
c xpi1=pi1^2,xpi2=pi2^2,pi12=pi1.pi2
*---------------------------------------------------------------------
      implicit none
      real*8 xpi1,xpi2,pi12,mi1,mi2,m0,xmi1,xmi2,xm0,f1,f2,q12,det
      complex*16 c0,c11,c12,c21,c22,c20,c23,
     &     b01,b02,b03,b11,b12,b13,b20,C0_
      xm0 = m0*m0
      xmi1 = mi1*mi1
      xmi2 = mi2*mi2
      q12=xpi1+xpi2-2d0*pi12
      f1=xpi1-xmi1+xm0
      f2=xpi2-xmi2+xm0
      det=xpi1*xpi2-pi12**2
c      if(dabs(det).lt.1d-1)then
c         write(6,*)'warning: potential instability in ccfunc2'
c         stop
c      endif
      if(dabs(q12).ne.0d0.and.dabs(q12).lt.1d-6)q12=0d0
      if(dabs(xpi1).ne.0d0.and.dabs(xpi1).lt.1d-6)xpi1=0d0
      if(dabs(xpi2).ne.0d0.and.dabs(xpi2).lt.1d-6)xpi2=0d0
      c0=C0_(xpi1,q12,xpi2,m0,mi1,mi2,0)
c
c UV divergent B0 and B1 functions
c UV divergences have been subtracted: N_eps=2/eps-gamma_E+ln(4pi)
c
      call bfunc(q12,mi1,mi2,b01,b11,b20)
      call bfunc(xpi2,m0,mi2,b02,b12,b20)
      call bfunc(xpi1,m0,mi1,b03,b13,b20)
*---------------------------------------------------------------------
*--entwicklungskoeff. des vektorintegrals:----------------------------
      c11=1d0/2d0/det*(pi12*(b01-b03)+xpi2*(b02-b01)
     $     +(-xpi2*f1+pi12*f2)*c0)

      c12=1d0/2d0/det*(pi12*(b01-b02)+xpi1*(b03-b01)
     $     +(-xpi1*f2+pi12*f1)*c0)
*---------------------------------------------------------------------
*------------entwicklungskoeff. des tensorintegrals--------------
      c20 = (b01+1d0)/4d0+xm0/2d0*c0+(f1*c11+f2*c12)/4d0

      c21 = 1d0/2d0/det*(pi12*(-b01-b11-b13+f2*c11)+
     $     xpi2*(b01+b11-2d0*c20-f1*c11))

      c22 = 1d0/2d0/det*(pi12*(b11-b12+f1*c12)+
     $     xpi1*(-b11-2d0*c20-f2*c12))

      c23 = 1d0/2d0/det*(pi12*(b11+2d0*c20+f2*c12)+
     $     xpi2*(-b11+b12-f1*c12))
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine ccfunc_ir(xmt,xpi1,xpi2,pi12,m0,mi1,mi2,
     $     c0,c11,c12,c21,c22,c20,c23)
*---------------------------------------------------------------------
c entwicklungskoeff. fuer das vektorielle und tensorielle dreipkt.-integral
c d0=k^2-m0^2, d1=(k+pi1)^2-m1^2, d2=(k+pi2)^2-m2^2
c xpi1=pi1^2,xpi2=pi2^2,pi12=pi1.pi2
c only for IR divergent C0 ! C0 is input !
*---------------------------------------------------------------------
      implicit none
      real*8 xpi1,xpi2,pi12,mi1,mi2,m0,xmi1,xmi2,xm0,f1,f2,q12,det
      real*8 xmt
      complex*16 c11,c12,c21,c22,c20,c23,b00,b10,
     &     b01,b02,b03,b11,b12,b13,b20
      complex*16 c0
      real*8 delta,muedr
      common/dimreg/delta,muedr
      xm0 = m0*m0
      xmi1 = mi1*mi1
      xmi2 = mi2*mi2
      q12=xpi1+xpi2-2d0*pi12
      f1=xpi1-xmi1+xm0
      f2=xpi2-xmi2+xm0
      det=xpi1*xpi2-pi12**2
c      if(dabs(det).lt.1d-1)then
c         write(6,*)'warning: potential instability in ccfunc_ir'
c         stop
c      endif
c
c IR and UV divergent B0 and B1 functions
c 
c -G[1+eps] (4 pi mu^2/mt^2)^eps 1/eps_IR has been subtracted
c UV divergences have been subtracted: N_eps=2/eps-gamma_E+ln(4pi)
      b00=dcmplx(-dlog(xmt**2/muedr**2),0d0)
      b10=-b00/2d0
      if(dabs(q12).lt.1d-5.and.dabs(mi1).lt.1d-5.and.
     $     dabs(mi2).lt.1d-5) then
         b01=b00
         b11=b10
      else
         call bfunc(q12,mi1,mi2,b01,b11,b20)
      endif
      if(dabs(xpi2).lt.1d-5.and.dabs(m0).lt.1d-5.and.
     $     dabs(mi2).lt.1d-5) then
         b02=b00
         b12=b10
      else
         call bfunc(xpi2,m0,mi2,b02,b12,b20)
      endif
      if(dabs(xpi1).lt.1d-5.and.dabs(m0).lt.1d-5.and.
     $     dabs(mi1).lt.1d-5) then
         b03=b00
         b13=b10
      else
         call bfunc(xpi1,m0,mi1,b03,b13,b20)
      endif
*---------------------------------------------------------------------
*--entwicklungskoeff. des vektorintegrals:----------------------------
      c11=1d0/2d0/det*(pi12*(b01-b03)+xpi2*(b02-b01)
     $     +(-xpi2*f1+pi12*f2)*c0)

      c12=1d0/2d0/det*(pi12*(b01-b02)+xpi1*(b03-b01)
     $     +(-xpi1*f2+pi12*f1)*c0)
*---------------------------------------------------------------------
*------------entwicklungskoeff. des tensorintegrals--------------
      c20 = (b01+1d0)/4d0+xm0/2d0*c0+(f1*c11+f2*c12)/4d0

      c21 = 1d0/2d0/det*(pi12*(-b01-b11-b13+f2*c11)+
     $     xpi2*(b01+b11-2d0*c20-f1*c11))

      c22 = 1d0/2d0/det*(pi12*(b11-b12+f1*c12)+
     $     xpi1*(-b11-2d0*c20-f2*c12))

      c23 = 1d0/2d0/det*(pi12*(b11+2d0*c20+f2*c12)+
     $     xpi2*(-b11+b12-f1*c12))
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      function C0_ir(xmt,q01,q12,q20,m0,m1,m2)
************************************************************************
* IR singular scalar 3-point function
* Dimensional regularization in d=4-2 eps_IR
* G[1+eps] (4 pi mu^2/mt^2)^eps_IR x poles is subtracted
*-----------------------------------------------------------------------
* d0=k^2-m0^2, d1=(k+q1)^2-m1^2, d2=(k+q2)^2-m2^2
* q01=q1^2, q20=q2^2, q12=(q1-q2)^2
* m0=0,q1^2=m1^2, and q2^2=m2^2 with m1,m2 arbitrary or
* q2^2 =/= m2^2 with m1=0.
*-----------------------------------------------------------------------
*	April 03, 2002, DW
************************************************************************
      implicit none
      real*8 xmt,q01,q12,q20,m0,m1,m2,m12,m22
      real*8 pi,eps,mm1,mm4,p12
      complex*16 C0_ir,ieps,spence,sc,sq,xs,xq
      C0_ir=dcmplx(0d0,0d0)
      m12=m1*m1
      m22=m2*m2
      p12=-q12+q01+q20
      pi=4d0*datan(1d0)                                  
      eps=1d-17
      ieps=dcmplx(0d0,eps)
      if (dabs(q20-m22).lt.1d-5) then
         if (dabs(m1*m2).gt.1d-5) then
            sc  = q12+dabs(q12)*ieps
            mm1 = dsqrt(m12)
            mm4 = dsqrt(m22)
            xs = -4d0*mm1*mm4/(sc-(mm1-mm4)**2)/
     &           (cdsqrt(1d0-4d0*mm1*mm4/(sc-(mm1-mm4)**2))+1d0)**2 
            C0_ir = xs/mm1/mm4/(1d0-xs**2)*(
     &           cdlog(xs)*(-cdlog(xs)/2d0+2d0*cdlog(1d0-xs**2)
     $           - dlog(xmt**2/mm1/mm4))
     &           - pi**2/6d0+spence(xs**2)+dlog(mm1/mm4)**2/2d0
     &           + spence(1d0-xs*mm1/mm4) + spence(1d0-xs*mm4/mm1))
         else if (dabs(m1).lt.1d-5.and.dabs(m2).lt.1d-5) then
            sc=p12-dabs(p12)*ieps
            xs=xmt**2/sc
            C0_ir = 1d0/sc*(-cdlog(xs)**2/2d0+pi**2/6d0)
         else if (dabs(m1).lt.1d-5.and.dabs(m2).gt.1d-5) then
            sc=p12-dabs(p12)*ieps
            xs=xmt**2/sc
            xq=q20/sc
            C0_ir = -1d0/sc*((cdlog(xs)**2-cdlog(xq)**2)/4d0
     $           +cdlog(xs)*cdlog(xq)/2d0-spence(1d0-xq))
         endif
      else if (dabs(m1).lt.1d-5.and.dabs(q20-m22).gt.1d-5.and.
     $        dabs(q12-m22).gt.1d-5) then
         sc=q12+dabs(q12)*ieps-m22
         xs=-sc/xmt**2
         sq=m22-q20-dabs(q12)*ieps
         xq=sq/xmt**2
         C0_ir=1d0/(sc+sq)*(cdlog(xs)**2/2d0
     $        -cdlog(xq)**2/2d0+spence(-q20/sq)-spence(q12/sc))
      else if (dabs(m1).lt.1d-5.and.dabs(q20-m22).gt.1d-5.and.
     $        dabs(q12-m22).lt.1d-5) then
         sq=m22-q20-ieps
         xq=sq/xmt**2
         C0_ir=1d0/sq*(-cdlog(xq)**2/2d0+spence(1d0-1d0/xq))
      endif
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine ddfunc(xpi1,xpi2,xpi3,pi12,pi13,pi23,
     $     m0,mi1,mi2,mi3,d1,d2,d3,c1,c2)
*---------------------------------------------------------------------
c entwicklungskoeff. fuer das vektorielle und tensorielle vierpkt.-integral
c d0=k^2-m0^2, d1=(k+pi1)^2-m1^2, d2=(k+pi2)^2-m2^2, d3=(k+pi3)^2-m3^2
c xpi1=pi1^2,xpi2=pi2^2,xpi3=pi3^2,pi12=pi1.pi2,pi13=pi1.pi3,pi23=pi2.pi3
c Up to rank 3 !
c Only IR finite D-functions ! 
*---------------------------------------------------------------------
      implicit none
      integer i
      real*8 xpi1,xpi2,xpi3,pi12,pi13,pi23,q12,q13,q23,
     $     mi1,mi2,mi3,m0,xmi1,xmi2,xmi3,xm0,cf1,cf2,cf3,det
      real*8 mat11,mat12,mat13,mat21,mat22,mat23,mat31,mat32,mat33
      complex*16 d1(0:3),d2(0:6),d3(0:3,0:3),c1(0:2,4),c2(0:3,4)
      complex*16 D0_
      complex*16 c1234(0:2),c2234(0:3),c1134(0:2),c2134(0:3),
     $     c1124(0:2),c2124(0:3),c1123(0:2),c2123(0:3),ccc
      complex*16 s11,s12,s13,s20,s211,s212,s213,s221,s222,s223,s231,
     $     s232,s233,s311,s312,s313,s321,s322,s323,s331,s332,s333,
     $     s3113,s3213,s3313,s310,s320,s330
      xm0 = m0*m0
      xmi1 = mi1*mi1
      xmi2 = mi2*mi2
      xmi3 = mi3*mi3
      q12=xpi1+xpi2-2d0*pi12
      q13=xpi1+xpi3-2d0*pi13
      q23=xpi2+xpi3-2d0*pi23
      cf1=xpi1-xmi1+xm0
      cf2=xpi2-xmi2+xm0
      cf3=xpi3-xmi3+xm0
c Gram determinant
      det=xpi1*xpi2*xpi3-xpi1*pi23**2-xpi2*pi13**2-xpi3*pi12**2+
     $     2d0*pi12*pi13*pi23
c      if(dabs(det).lt.1d-1)then
c         write(6,*)'warning: potential instability in ddfunc'
c         write(6,*)xpi1,xpi2,xpi3,pi23,pi13,pi12
c         stop
c      endif
c matrix elements of the inverse Gram matrix x det
      mat11=xpi2*xpi3-pi23**2
      mat12=pi13*pi23-pi12*xpi3
      mat13=pi12*pi23-pi13*xpi2
      mat21=mat12
      mat22=xpi1*xpi3-pi13**2
      mat23=pi12*pi13-pi23*xpi1
      mat31=mat13
      mat32=mat23
      mat33=xpi1*xpi2-pi12**2
*
      d1(0)=D0_(xpi1,q12,q23,xpi3,xpi2,q13,
     $     m0,mi1,mi2,mi3,0)
      call ccfunc2(q12,q13,pi23-pi12-pi13+xpi1,mi1,mi2,mi3,
     $     c1234(0),c1234(1),c1234(2),
     $     c2234(1),c2234(2),c2234(0),c2234(3))
      call ccfunc2(xpi2,xpi3,pi23,m0,mi2,mi3,
     $     c1134(0),c1134(1),c1134(2),
     $     c2134(1),c2134(2),c2134(0),c2134(3))
      call ccfunc2(xpi1,xpi3,pi13,m0,mi1,mi3,
     $     c1124(0),c1124(1),c1124(2),
     $     c2124(1),c2124(2),c2124(0),c2124(3))
      call ccfunc2(xpi1,xpi2,pi12,m0,mi1,mi2,
     $     c1123(0),c1123(1),c1123(2),
     $     c2123(1),c2123(2),c2123(0),c2123(3))
      do i=0,2
         c1(i,1)=c1234(i)
         c1(i,2)=c1134(i)
         c1(i,3)=c1124(i)
         c1(i,4)=c1123(i)
      end do
      do i=0,3
         c2(i,1)=c2234(i)
         c2(i,2)=c2134(i)
         c2(i,3)=c2124(i)
         c2(i,4)=c2123(i)
      end do
*---------------------------------------------------------------------
*--entwicklungskoeff. des vektorintegrals:----------------------------
      s11 = (c1134(0)-c1234(0)-cf1*d1(0))/2d0
      s12 = (c1124(0)-c1234(0)-cf2*d1(0))/2d0
      s13 = (c1123(0)-c1234(0)-cf3*d1(0))/2d0
      d1(1) = (mat11*s11+mat12*s12+mat13*s13)/det
      d1(2) = (mat21*s11+mat22*s12+mat23*s13)/det
      d1(3) = (mat31*s11+mat32*s12+mat33*s13)/det
*---------------------------------------------------------------------
*--entwicklungskoeff. des tensorintegrals (Rank 2):-------------------
      ccc = c1234(1)+c1234(2)+c1234(0)
      s20 = c1234(0)+xm0*d1(0)
      s211 = (ccc-cf1*d1(1))/2d0
      s212 = (c1134(1)-c1234(1)-cf1*d1(2))/2d0
      s213 = (c1134(2)-c1234(2)-cf1*d1(3))/2d0
      s221 = (c1124(1)+ccc-cf2*d1(1))/2d0
      s222 = -(c1234(1)+cf2*d1(2))/2d0
      s223 = (c1124(2)-c1234(2)-cf2*d1(3))/2d0
      s231 = (c1123(1)+ccc-cf3*d1(1))/2d0
      s232 = (c1123(2)-c1234(1)-cf3*d1(2))/2d0
      s233 = -(c1234(2)+cf3*d1(3))/2d0
      d2(0) = s20-s211-s222-s233
      d2(1) = (mat11*(s211-d2(0))+mat12*s221+mat13*s231)/det
      d2(2) = (mat21*s212+mat22*(s222-d2(0))+mat23*s232)/det
      d2(3) = (mat31*s213+mat32*s223+mat33*(s233-d2(0)))/det
c---  >   d2(1,2)
      d2(4) = (mat11*s212+mat12*(s222-d2(0))+mat13*s232)/det
c---  >   d2(1,3)
      d2(5) = (mat11*s213+mat12*s223+mat13*(s233-d2(0)))/det
c---  >   d2(2,3)=d2(2,1) with m2<-->m4 for square boxes
      d2(6) = (mat21*s213+mat22*s223+mat23*(s233-d2(0)))/det
*---------------------------------------------------------------------
*--entwicklungskoeff. des tensorintegrals (Rank 3):-------------------
      s310 = (c2134(0)-c2234(0)-cf1*d2(0))/2d0
      s320 = (c2124(0)-c2234(0)-cf2*d2(0))/2d0
      s330 = (c2123(0)-c2234(0)-cf3*d2(0))/2d0
      ccc = ccc+(c2234(1)+c2234(2)-c1234(0))/2d0+c2234(3)
      s311 = -ccc-cf1*d2(1)/2d0
      s321 = (c2124(1)-cf2*d2(1))/2d0-ccc
      s331 = (c2123(1)-cf3*d2(1))/2d0-ccc
      s312 = (c2134(1)-c2234(1)-cf1*d2(2))/2d0
      s322 = -(c2234(1)+cf2*d2(2))/2d0
      s332 = (c2123(2)-c2234(1)-cf3*d2(2))/2d0
      s313 = (c2134(2)-c2234(2)-cf1*d2(3))/2d0
      s323 = (c2124(2)-c2234(2)-cf2*d2(3))/2d0
      s333 = -(c2234(2)+cf3*d2(3))/2d0
      s3113 = (c2234(3)+c2234(2)+c1234(2)-cf1*d2(5))/2d0
      s3213 = (c2124(3)+c2234(3)+c2234(2)+c1234(2)-cf2*d2(5))/2d0
      s3313 = (c2234(3)+c2234(2)+c1234(2)-cf3*d2(5))/2d0
      d3(0,1) = (mat11*s310+mat12*s320+mat13*s330)/det
      d3(0,2) = (mat21*s310+mat22*s320+mat23*s330)/det
      d3(0,3) = (mat31*s310+mat32*s320+mat33*s330)/det
      d3(1,1) = (mat11*(s311-2d0*d3(0,1))+mat12*s321+mat13*s331)/det
      d3(3,3) = (mat31*s313+mat32*s323+mat33*(s333-2d0*d3(0,3)))/det
      d3(1,2) = (mat21*(s311-2d0*d3(0,1))+mat22*s321+mat23*s331)/det
      d3(1,3) = (mat31*(s311-2d0*d3(0,1))+mat32*s321+mat33*s331)/det
      d3(2,1) = (mat11*s312+mat12*(s322-2d0*d3(0,2))+mat13*s332)/det
      d3(2,2) = (mat21*s312+mat22*(s322-2d0*d3(0,2))+mat23*s332)/det
      d3(2,3) = (mat31*s312+mat32*(s322-2d0*d3(0,2))+mat33*s332)/det
      d3(3,1) = (mat11*s313+mat12*s323+mat13*(s333-2d0*d3(0,3)))/det
      d3(3,2) = (mat21*s313+mat22*s323+mat23*(s333-2d0*d3(0,3)))/det
      d3(0,0) = (mat21*(s3113-d3(0,3))+mat22*s3213+mat23*
     1     (s3313-d3(0,1)))/det
 99   continue
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine ddfunc_ir(xmt,xpi1,xpi2,xpi3,pi12,pi13,pi23,
     $     m0,mi1,mi2,mi3,d1,d2,d3,c1,c2)
*---------------------------------------------------------------------
c entwicklungskoeff. fuer das vektorielle und tensorielle vierpkt.-integral
c d0=k^2-m0^2, d1=(k+pi1)^2-m1^2, d2=(k+pi2)^2-m2^2, d3=(k+pi3)^2-m3^2
c xpi1=pi1^2,xpi2=pi2^2,xpi3=pi3^2,pi12=pi1.pi2,pi13=pi1.pi3,pi23=pi2.pi3
c Up to rank 3, only IR div. box diagrams !
c always m0=0 and xpi1=mi1^2
c d1(0) is IR div. and input ! 
*---------------------------------------------------------------------
      implicit none
      integer i
      real*8 xpi1,xpi2,xpi3,pi12,pi13,pi23,q12,q13,q23,
     $     mi1,mi2,mi3,m0,xmi1,xmi2,xmi3,xm0,cf1,cf2,cf3,det,xmt
      real*8 mat11,mat12,mat13,mat21,mat22,mat23,mat31,mat32,mat33
      complex*16 d1(0:3),d2(0:6),d3(0:3,0:3),c1(0:2,4),c2(0:3,4),c0_ir
      complex*16 c1234(0:2),c2234(0:3),c1134(0:2),c2134(0:3),
     $     c1124(0:2),c2124(0:3),c1123(0:2),c2123(0:3),ccc
      complex*16 s11,s12,s13,s20,s211,s212,s213,s221,s222,s223,s231,
     $     s232,s233,s311,s312,s313,s321,s322,s323,s331,s332,s333,
     $     s3113,s3213,s3313,s310,s320,s330
      xm0 = m0*m0
      xmi1 = mi1*mi1
      xmi2 = mi2*mi2
      xmi3 = mi3*mi3
      q12=xpi1+xpi2-2d0*pi12
      q13=xpi1+xpi3-2d0*pi13
      q23=xpi2+xpi3-2d0*pi23
      cf1=xpi1-xmi1+xm0
      cf2=xpi2-xmi2+xm0
      cf3=xpi3-xmi3+xm0
c Gram determinant
      det=xpi1*xpi2*xpi3-xpi1*pi23**2-xpi2*pi13**2-xpi3*pi12**2+
     $     2d0*pi12*pi13*pi23
c      if(dabs(det).lt.1d-1)then
c         write(6,*)'warning: potential instability in ddfunc_ir'
c         write(6,*)xpi1,xpi2,xpi3,pi23,pi13,pi12
c         stop
c      endif
c matrix elements of the inverse Gram matrix x det
      mat11=xpi2*xpi3-pi23**2
      mat12=pi13*pi23-pi12*xpi3
      mat13=pi12*pi23-pi13*xpi2
      mat21=mat12
      mat22=xpi1*xpi3-pi13**2
      mat23=pi12*pi13-pi23*xpi1
      mat31=mat13
      mat32=mat23
      mat33=xpi1*xpi2-pi12**2
*
      call ccfunc2(q12,q13,pi23-pi12-pi13+xpi1,mi1,mi2,mi3,
     $     c1234(0),c1234(1),c1234(2),
     $     c2234(1),c2234(2),c2234(0),c2234(3))
c check for IR divergences:
      if (dabs(xpi2-xmi2).lt.1d-5.and.dabs(xpi3-xmi3).lt.1d-5) then
         c1134(0)=C0_ir(xmt,xpi2,q23,xpi3,m0,mi2,mi3)
         call ccfunc_ir(xmt,xpi2,xpi3,pi23,m0,mi2,mi3,
     $        c1134(0),c1134(1),c1134(2),
     $        c2134(1),c2134(2),c2134(0),c2134(3))
      else if (dabs(xpi2-xmi2).lt.1d-5.and.dabs(xmi3).lt.1d-5) then
         c1134(0)=C0_ir(xmt,xpi2,q23,xpi3,m0,mi2,mi3)
         call ccfunc_ir(xmt,xpi2,xpi3,pi23,m0,mi2,mi3,
     $        c1134(0),c1134(1),c1134(2),
     $        c2134(1),c2134(2),c2134(0),c2134(3))
      else if (dabs(xpi3-xmi3).lt.1d-5.and.dabs(xmi2).lt.1d-5) then
         c1134(0)=C0_ir(xmt,xpi2,q23,xpi3,m0,mi2,mi3)
         call ccfunc_ir(xmt,xpi2,xpi3,pi23,m0,mi2,mi3,
     $        c1134(0),c1134(1),c1134(2),
     $        c2134(1),c2134(2),c2134(0),c2134(3))
      else if (dabs(xpi2).lt.1d-5.and.dabs(xmi2).lt.1d-5) then
         c1134(0)=C0_ir(xmt,xpi2,q23,xpi3,m0,mi2,mi3)
         call ccfunc_ir(xmt,xpi2,xpi3,pi23,m0,mi2,mi3,
     $        c1134(0),c1134(1),c1134(2),
     $        c2134(1),c2134(2),c2134(0),c2134(3))
      else
         call ccfunc2(xpi2,xpi3,pi23,m0,mi2,mi3,
     $        c1134(0),c1134(1),c1134(2),
     $        c2134(1),c2134(2),c2134(0),c2134(3))
      endif
      if (dabs(xpi3-xmi3).lt.1d-5) then
         c1124(0)=C0_ir(xmt,xpi1,q13,xpi3,m0,mi1,mi3)
         call ccfunc_ir(xmt,xpi1,xpi3,pi13,m0,mi1,mi3,
     $     c1124(0),c1124(1),c1124(2),
     $     c2124(1),c2124(2),c2124(0),c2124(3))
      else if (dabs(xmi1).lt.1d-5) then
         c1124(0)=C0_ir(xmt,xpi1,q13,xpi3,m0,mi1,mi3)
         call ccfunc_ir(xmt,xpi1,xpi3,pi13,m0,mi1,mi3,
     $        c1124(0),c1124(1),c1124(2),
     $        c2124(1),c2124(2),c2124(0),c2124(3))
      else
         call ccfunc2(xpi1,xpi3,pi13,m0,mi1,mi3,
     $        c1124(0),c1124(1),c1124(2),
     $        c2124(1),c2124(2),c2124(0),c2124(3))
      end if
      if (dabs(xpi2-xmi2).lt.1d-5) then
         c1123(0)=C0_ir(xmt,xpi1,q12,xpi2,m0,mi1,mi2)
         call ccfunc_ir(xmt,xpi1,xpi2,pi12,m0,mi1,mi2,
     $        c1123(0),c1123(1),c1123(2),
     $        c2123(1),c2123(2),c2123(0),c2123(3))
      else if (dabs(xmi1).lt.1d-5) then
         c1123(0)=C0_ir(xmt,xpi1,q12,xpi2,m0,mi1,mi2)
         call ccfunc_ir(xmt,xpi1,xpi2,pi12,m0,mi1,mi2,
     $        c1123(0),c1123(1),c1123(2),
     $        c2123(1),c2123(2),c2123(0),c2123(3))
      else 
         call ccfunc2(xpi1,xpi2,pi12,m0,mi1,mi2,
     $        c1123(0),c1123(1),c1123(2),
     $        c2123(1),c2123(2),c2123(0),c2123(3))
      end if
      do i=0,2
         c1(i,1)=c1234(i)
         c1(i,2)=c1134(i)
         c1(i,3)=c1124(i)
         c1(i,4)=c1123(i)
      end do
      do i=0,3
         c2(i,1)=c2234(i)
         c2(i,2)=c2134(i)
         c2(i,3)=c2124(i)
         c2(i,4)=c2123(i)
      end do
*---------------------------------------------------------------------
*--entwicklungskoeff. des vektorintegrals:----------------------------
      s11 = (c1134(0)-c1234(0)-cf1*d1(0))/2d0
      s12 = (c1124(0)-c1234(0)-cf2*d1(0))/2d0
      s13 = (c1123(0)-c1234(0)-cf3*d1(0))/2d0
      d1(1) = (mat11*s11+mat12*s12+mat13*s13)/det
      d1(2) = (mat21*s11+mat22*s12+mat23*s13)/det
      d1(3) = (mat31*s11+mat32*s12+mat33*s13)/det
*---------------------------------------------------------------------
*--entwicklungskoeff. des tensorintegrals (Rank 2):-------------------
      ccc = c1234(1)+c1234(2)+c1234(0)
      s20 = c1234(0)+xm0*d1(0)
      s211 = (ccc-cf1*d1(1))/2d0
      s212 = (c1134(1)-c1234(1)-cf1*d1(2))/2d0
      s213 = (c1134(2)-c1234(2)-cf1*d1(3))/2d0
      s221 = (c1124(1)+ccc-cf2*d1(1))/2d0
      s222 = -(c1234(1)+cf2*d1(2))/2d0
      s223 = (c1124(2)-c1234(2)-cf2*d1(3))/2d0
      s231 = (c1123(1)+ccc-cf3*d1(1))/2d0
      s232 = (c1123(2)-c1234(1)-cf3*d1(2))/2d0
      s233 = -(c1234(2)+cf3*d1(3))/2d0
      d2(0) = s20-s211-s222-s233
      d2(1) = (mat11*(s211-d2(0))+mat12*s221+mat13*s231)/det
      d2(2) = (mat21*s212+mat22*(s222-d2(0))+mat23*s232)/det
      d2(3) = (mat31*s213+mat32*s223+mat33*(s233-d2(0)))/det
c---  >   d2(1,2)
      d2(4) = (mat11*s212+mat12*(s222-d2(0))+mat13*s232)/det
c---  >   d2(1,3)
      d2(5) = (mat11*s213+mat12*s223+mat13*(s233-d2(0)))/det
c---  >   d2(2,3)=d2(2,1) with m2<-->m4 for square boxes
      d2(6) = (mat21*s213+mat22*s223+mat23*(s233-d2(0)))/det
*---------------------------------------------------------------------
*--entwicklungskoeff. des tensorintegrals (Rank 3):-------------------
      s310 = (c2134(0)-c2234(0)-cf1*d2(0))/2d0
      s320 = (c2124(0)-c2234(0)-cf2*d2(0))/2d0
      s330 = (c2123(0)-c2234(0)-cf3*d2(0))/2d0
      ccc = ccc+(c2234(1)+c2234(2)-c1234(0))/2d0+c2234(3)
      s311 = -ccc-cf1*d2(1)/2d0
      s321 = (c2124(1)-cf2*d2(1))/2d0-ccc
      s331 = (c2123(1)-cf3*d2(1))/2d0-ccc
      s312 = (c2134(1)-c2234(1)-cf1*d2(2))/2d0
      s322 = -(c2234(1)+cf2*d2(2))/2d0
      s332 = (c2123(2)-c2234(1)-cf3*d2(2))/2d0
      s313 = (c2134(2)-c2234(2)-cf1*d2(3))/2d0
      s323 = (c2124(2)-c2234(2)-cf2*d2(3))/2d0
      s333 = -(c2234(2)+cf3*d2(3))/2d0
      s3113 = (c2234(3)+c2234(2)+c1234(2)-cf1*d2(5))/2d0
      s3213 = (c2124(3)+c2234(3)+c2234(2)+c1234(2)-cf2*d2(5))/2d0
      s3313 = (c2234(3)+c2234(2)+c1234(2)-cf3*d2(5))/2d0
      d3(0,1) = (mat11*s310+mat12*s320+mat13*s330)/det
      d3(0,2) = (mat21*s310+mat22*s320+mat23*s330)/det
      d3(0,3) = (mat31*s310+mat32*s320+mat33*s330)/det
      d3(1,1) = (mat11*(s311-2d0*d3(0,1))+mat12*s321+mat13*s331)/det
      d3(3,3) = (mat31*s313+mat32*s323+mat33*(s333-2d0*d3(0,3)))/det
      d3(1,2) = (mat21*(s311-2d0*d3(0,1))+mat22*s321+mat23*s331)/det
      d3(1,3) = (mat31*(s311-2d0*d3(0,1))+mat32*s321+mat33*s331)/det
      d3(2,1) = (mat11*s312+mat12*(s322-2d0*d3(0,2))+mat13*s332)/det
      d3(2,2) = (mat21*s312+mat22*(s322-2d0*d3(0,2))+mat23*s332)/det
      d3(2,3) = (mat31*s312+mat32*(s322-2d0*d3(0,2))+mat33*s332)/det
      d3(3,1) = (mat11*s313+mat12*s323+mat13*(s333-2d0*d3(0,3)))/det
      d3(3,2) = (mat21*s313+mat22*s323+mat23*(s333-2d0*d3(0,3)))/det
      d3(0,0) = (mat21*(s3113-d3(0,3))+mat22*s3213+mat23*
     1     (s3313-d3(0,1)))/det
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      function D0_ir(xmt,p1,p2,p3,p4,p12,p23,m2,m3,m4)
************************************************************************
* IR singular scalar 4-point function
* Dimensional regularization in d=4-2 eps_IR
* G[1+eps] (4 pi mu^2/mt^2)^eps_IR x poles is subtracted
*-----------------------------------------------------------------------
* d0=k^2, d1=(k+p1)^2-m2^2, d2=(k+p1+p2)^2-m3^2, d3=(k-p4)^2-m4^2
* p12=(p1+p2)^2,p23=(p2+p3)^2
* always m1=0 and p1^2=m2^2 and p4^2=m4^2. 
*-----------------------------------------------------------------------
*	April 16, 2002, DW
************************************************************************
      implicit none
      real*8 xmt,m2,m3,m4,m22,m32,m42,mm2,mm3,mm4
      real*8 p1,p2,p3,p4,p12,p23,q1,q2,q3,q4,q12,q23
      real*8 pi,eps
      complex*16 D0_ir,ieps,spence,sc,tc,q2c,q3c,xs,x2,x3
      complex*16 t1,t2,t4,w1,w2,lap,lam,tau,omega
      D0_ir=dcmplx(0d0,0d0)
      mm2=m2    
      mm3=m3    
      mm4=m4    
      m22=m2*m2 
      m32=m3*m3 
      m42=m4*m4 
      q1=p1 
      q2=p2   
      q3=p3   
      q4=p4
      q12=p12   
      q23=p23
      pi=4d0*datan(1d0)                                  
      eps=1d-17
      ieps=dcmplx(0d0,eps)
      sc  = q23+dabs(q23)*ieps
      tc  = q12+dabs(q12)*ieps
      q2c = q2+dabs(q2)*ieps
      q3c = q3+dabs(q3)*ieps
      if (dabs(mm2*mm3*mm4).ne.0d0.and.dabs(q12-m32).gt.1d-5) then
         xs = -4d0*mm2*mm4/(sc-(mm2-mm4)**2) /
     &        (cdsqrt(1d0-4d0*mm2*mm4/( sc-(mm2-mm4)**2))+1d0 )**2 
         if (q2.ne.(mm2-mm3)**2) then
            x2 = -4d0*mm2*mm3/(q2c-(mm2-mm3)**2) /
     &           (cdsqrt(1d0-4d0*mm2*mm3/( q2c-(mm2-mm3)**2))+1d0 )**2 
         else
            x2 = 1d0
         endif
         if (q3.ne.(mm4-mm3)**2) then
            x3 = -4d0*mm4*mm3/(q3c-(mm4-mm3)**2) /
     &           (cdsqrt(1d0-4d0*mm4*mm3/( q3c-(mm4-mm3)**2))+1d0 )**2 
         else
            x3 = 1d0
         endif
         D0_ir = xs/mm2/mm4/(q12-m32)/(1d0-xs**2)*(
     &        2d0*cdlog(xs)*(cdlog(1d0-xs**2)
     $        -cdlog(xmt*mm3/(m32-tc)))
     &        +pi**2/2d0+spence(xs**2)+cdlog(x2)**2+cdlog(x3)**2
     &        -spence(xs*x2*x3)-(cdlog(xs)+cdlog(x2)+cdlog(x3))
     $        *cdlog(1d0-xs*x2*x3)
     &        -spence(xs*x2/x3)-(cdlog(xs)+cdlog(x2)-cdlog(x3))
     $        *cdlog(1d0-xs*x2/x3)
     &        -spence(xs/x2*x3)-(cdlog(xs)-cdlog(x2)+cdlog(x3))
     $        *cdlog(1d0-xs/x2*x3)
     &        -spence(xs/x2/x3)-(cdlog(xs)-cdlog(x2)-cdlog(x3))
     $        *cdlog(1d0-xs/x2/x3) )
      elseif (dabs(mm2).eq.0d0.and.dabs(q12-m32).gt.1d-5.and.
     $        dabs(m42-q23).gt.1d-5.and.dabs(q2-m32).gt.1d-5.and.
     $        dabs(mm3*mm4).ne.0d0) then
         if (dabs(q3).gt.1d-8) then
            t1=m32-q2c
            t2=m42-sc
            w2=tc-m32
            lap=(1d0+cdsqrt(1d0-4d0*mm3*mm4/q3c))/2d0
            lam=(1d0-cdsqrt(1d0-4d0*mm3*mm4/q3c))/2d0
            D0_ir=-1d0/(q12-m32)/(m42-q23)*(-5d0/6d0*pi**2
     $           +cdlog(w2/xmt**2)**2+cdlog(t2/xmt**2)**2
     $           -cdlog(t1/xmt**2)**2
     $           +2d0*cdlog((w2+t1)/t2)*cdlog(t1/w2)
     $           +2d0*cdlog((t2-t1)/w2)*cdlog(t1/t2)
     $           -2d0*spence((t2-t1-w2)/t2)-2d0*spence((w2+t1-t2)/w2)
     $           +2d0*spence(t1*(w2+t1-t2)/w2/t2)
     $           -cdlog(t2/t1)*cdlog(q3c/xmt**2)+spence(1d0/lap)
     $           -cdlog(t2/t1)*cdlog((-t1-lap*(t2-t1))/(t2-t1))
     $           +spence(t2/(lap*(t2-t1)+t1))
     $           -spence(t1/(lap*(t2-t1)+t1))
     $           +spence(1d0/lam)
     $           -cdlog(t2/t1)*cdlog((-t1-lam*(t2-t1))/(t2-t1))
     $           +spence(t2/(lam*(t2-t1)+t1))
     $           -spence(t1/(lam*(t2-t1)+t1)))
         else
            t2=m32-tc
            t4=m42-sc
            w1=q2c-m32
            D0_ir=1d0/(m32-q12)/(m42-q23)*(cdlog(t2/xmt**2)**2
     $           +cdlog(t4/xmt**2)**2-cdlog(w1/xmt**2)**2
     $           +3d0/2d0*pi**2+2d0*cdlog((t2+w1)/t4)
     $           *cdlog(t4/(t2+t4+w1))+2d0*cdlog((t4+w1)/t2)
     $           *cdlog(t2/(t2+t4+w1))-2d0*spence((t2+t4+w1)/t4)
     $           -2d0*spence((t2+t4+w1)/t2)
     $           -2d0*spence((t2+w1)*(t4+w1)/t2/t4))
         endif
      elseif (dabs(mm2).eq.0d0.and.dabs(q12-m32).gt.1d-5.and.
     $        dabs(m42-q23).gt.1d-5.and.dabs(q2-m32).lt.1d-5.and.
     $        dabs(mm3*mm4).ne.0d0) then
         if (dabs(q3).ne.0d0) then
            t1=m32-tc
            t4=m42-sc
            lap=(1d0+cdsqrt(1d0-4d0*mm3*mm4/q3c))/2d0
            lam=(1d0-cdsqrt(1d0-4d0*mm3*mm4/q3c))/2d0
            D0_ir=1d0/(m32-q12)/(m42-q23)*(cdlog(t1/xmt**2)**2
     $           +cdlog(t4/xmt**2)**2-cdlog(t4/t1)**2
     $           -2d0/3d0*pi**2+2d0*spence(1d0/lap)
     $           +2d0*spence(1d0/lam))
         endif
      elseif (dabs(mm2).eq.0d0.and.dabs(mm4).eq.0d0.and.
     $        dabs(q12-m32).gt.1d-5.and.dabs(q23).gt.1d-5.and.
     $        dabs(mm3).ne.0d0) then
         tau=m32-tc
         if (dabs(q2-m32).gt.1d-5) then
            omega=q2c-m32
         else if (dabs(q3-m32).gt.1d-5) then
            omega=q3c-m32
         end if
         D0_ir=-1d0/q23/(m32-q12)*(2d0*cdlog(tau/xmt**2)
     $        *cdlog(sc/xmt**2)-cdlog(omega/xmt**2)**2+pi**2/3d0
     $        -2d0*spence(1d0+omega/tau))
      endif
      return
      end            
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine eefunc_ir(np,xmt,xpi1,xpi2,xpi3,xpi4,pi12,pi13,
     $     pi14,pi23,pi24,pi34,m0,mi1,mi2,mi3,mi4,e1,e2,e3,d1,d2,d3)
*---------------------------------------------------------------------
c entwicklungskoeff. fuer das vektorielle und tensorielle fuenfpkt.-integral
c d0=k^2-m0^2, d1=(k+pi1)^2-m1^2, d2=(k+pi2)^2-m2^2, d3=(k+pi3)^2-m3^2,
c d4=(k+pi4)^2-m4^2
c xpi1=pi1^2,xpi2=pi2^2,xpi3=pi3^2,xpi4=pi4^2,pi12=pi1.pi2,pi13=pi1.pi3,
c pi14=pi1.pi4,pi23=pi2.pi3,pi24=pi2.pi4,pi34=pi3.pi4
c Up to rank 3, only IR div. pentagon diagrams !
c always m0=0 and xpi1=mi1^2 and xpi2=mi2^2 and mi3=mt and mi4=mt
c e1(0) is IR div. and input ! 
c D.W., June 18, 2002
*---------------------------------------------------------------------
      implicit none
      integer i,j,np
      real*8 xpi1,xpi2,xpi3,xpi4,pi12,pi13,pi14,pi23,pi24,pi34,
     $     q12,q13,q14,q23,q24,q34,mi1,mi2,mi3,mi4,m0,
     $     xmi1,xmi2,xmi3,xmi4,xm0,cf3,cf4,det,xmt
      real*8 mat11,mat12,mat13,mat14,mat21,mat22,mat23,mat24,
     $     mat31,mat32,mat33,mat34,mat41,mat42,mat43,mat44
      complex*16 e1(0:4),e2(4,4),e3(0:4,4)
      complex*16 d1(0:3),d2(0:6),d3(0:3,0:3)
      complex*16 d11(0:3),d21(0:6),d31(0:3,0:3),D0_ir
      complex*16 d12(0:3),d22(0:6),d32(0:3,0:3)
      complex*16 d13(0:3),d23(0:6),d33(0:3,0:3)
      complex*16 d14(0:3),d24(0:6),d34(0:3,0:3)
      complex*16 d15(0:3),d25(0:6),d35(0:3,0:3)
      complex*16 c1(0:2,4),c2(0:3,4),ddd,diff(4)
      complex*16 s11,s12,s13,s14,s211,s212,s213,s214,s221,s222,s223,
     $     s224,s231,s232,s233,s234,s241,s242,s243,s244,s311,s312,s313,
     $     s314,s321,s322,s323,s324,s331,s332,s333,s334,s341,s342,s343,
     $     s344,s3121,s3122,s3341,s3342,s3343,s3344,s3123,s3124
      integer nlim
      real*8 dets

      xm0 = m0*m0
      xmi1 = mi1*mi1
      xmi2 = mi2*mi2
      xmi3 = mi3*mi3
      xmi4 = mi4*mi4
      q12=xpi1+xpi2-2d0*pi12
      q13=xpi1+xpi3-2d0*pi13
      q14=xpi1+xpi4-2d0*pi14
      q23=xpi2+xpi3-2d0*pi23
      q24=xpi2+xpi4-2d0*pi24
      q34=xpi3+xpi4-2d0*pi34
      cf3=xpi3-xmi3+xm0
      cf4=xpi4-xmi4+xm0
c Gram determinant
      det=(xpi1*xpi2*xpi3*xpi4-xpi3*xpi4*pi12**2-
     $     xpi2*xpi4*pi13**2-xpi2*xpi3*pi14**2+
     $     2d0*xpi4*pi12*pi13*pi23-xpi1*xpi4*pi23**2+
     $     pi14**2*pi23**2+2d0*xpi3*pi12*pi14*pi24-
     $     2d0*pi13*pi14*pi23*pi24-xpi1*xpi3*pi24**2+
     $     pi13**2*pi24**2+2d0*xpi2*pi13*pi14*pi34-
     $     2d0*pi12*pi14*pi23*pi34-
     $     2d0*pi12*pi13*pi24*pi34+
     $     2d0*xpi1*pi23*pi24*pi34-xpi1*xpi2*pi34**2+
     $     pi12**2*pi34**2)
      nlim=0
c      if(nlim.eq.0) then
cc         det=-kappas*kappai(1)*kappai(2)*
cc     $        kappai(3)*kappai(4)*kappai(5)
         dets=1d0
c      else
c         det=-kappas
c         dets=kappai(1)*kappai(2)*kappai(3)*kappai(4)*kappai(5)
c      endif
c
c      if(dabs(det).lt.1d-1)then
c         write(6,*)'warning: potential instability in eefunc_ir'
c         stop
c      endif
c matrix elements of the inverse Gram matrix x det
      mat11=(xpi2*xpi3*xpi4-xpi4*pi23**2-xpi3*pi24**2+
     $     2d0*pi23*pi24*pi34-xpi2*pi34**2)
      mat12=(-(xpi3*xpi4*pi12)+xpi4*pi13*pi23+
     $     xpi3*pi14*pi24-pi14*pi23*pi34-
     $     pi13*pi24*pi34+pi12*pi34**2)
      mat13=(-(xpi2*xpi4*pi13)+xpi4*pi12*pi23-
     $     pi14*pi23*pi24+pi13*pi24**2+
     $     xpi2*pi14*pi34-pi12*pi24*pi34)
      mat14=(-(xpi2*xpi3*pi14)+pi14*pi23**2+
     $     xpi3*pi12*pi24-pi13*pi23*pi24+
     $     xpi2*pi13*pi34-pi12*pi23*pi34)
      mat21=mat12
      mat22=(xpi1*xpi3*xpi4-xpi4*pi13**2-xpi3*pi14**2+
     $     2d0*pi13*pi14*pi34-xpi1*pi34**2)
      mat23=(xpi4*pi12*pi13-xpi1*xpi4*pi23+
     $     pi14**2*pi23-pi13*pi14*pi24-
     $     pi12*pi14*pi34+xpi1*pi24*pi34)
      mat24=(xpi3*pi12*pi14-pi13*pi14*pi23-
     $     xpi1*xpi3*pi24+pi13**2*pi24-
     $     pi12*pi13*pi34+xpi1*pi23*pi34)
      mat31=mat13
      mat32=mat23
      mat33=(xpi1*xpi2*xpi4-xpi4*pi12**2-xpi2*pi14**2+
     $     2d0*pi12*pi14*pi24-xpi1*pi24**2)
      mat34=(xpi2*pi13*pi14-pi12*pi14*pi23-
     $     pi12*pi13*pi24+xpi1*pi23*pi24-
     $     xpi1*xpi2*pi34+pi12**2*pi34)
      mat41=mat14
      mat42=mat24
      mat43=mat34
      mat44=(xpi1*xpi2*xpi3-xpi3*pi12**2-xpi2*pi13**2+
     $     2d0*pi12*pi13*pi23-xpi1*pi23**2)
*
      call ddfunc(q12,q13,q14,pi23-pi12-pi13+xpi1,
     $     pi24-pi12-pi14+xpi1,xpi1+pi34-pi14-pi13,
     $     mi1,mi2,mi3,mi4,d11,d21,d31,c1,c2)
      do i=0,3
         d1(i)=d11(i)
      enddo
      do i=0,6
         d2(i)=d21(i)
      enddo
      do i=0,3
         do j=0,3
            d3(i,j)=d31(i,j)
         enddo
      enddo
c only D0(2) of P1 is IR divergent:
      if(np.eq.1)then
         d12(0)=D0_ir(xmt,xpi2,xpi3,q34,q24,q23,xpi4,mi2,mi3,mi4)
         call ddfunc_ir(xmt,xpi2,xpi3,xpi4,pi23,pi24,pi34,
     $        m0,mi2,mi3,mi4,d12,d22,d32,c1,c2)
      else
         call ddfunc(xpi2,xpi3,xpi4,pi23,pi24,pi34,
     $        m0,mi2,mi3,mi4,d12,d22,d32,c1,c2)
      endif
c only D0(3) of P1,P5,P6 are IR divergent:
      if(np.eq.1.or.np.eq.5.or.np.eq.6)then
         d13(0)=D0_ir(xmt,xpi1,xpi4,q34,q13,q14,xpi3,mi1,mi4,mi3)
         call ddfunc_ir(xmt,xpi1,xpi3,xpi4,pi13,pi14,pi34,
     $        m0,mi1,mi3,mi4,d13,d23,d33,c1,c2)
      else
         call ddfunc(xpi1,xpi3,xpi4,pi13,pi14,pi34,
     $        m0,mi1,mi3,mi4,d13,d23,d33,c1,c2)
      endif
c all D0(4) are IR divergent:
      d14(0)=D0_ir(xmt,xpi1,q14,q24,xpi2,xpi4,q12,mi1,mi4,mi2)
      call ddfunc_ir(xmt,xpi1,xpi2,xpi4,pi12,pi14,pi24,
     $     m0,mi1,mi2,mi4,d14,d24,d34,c1,c2)
c all D0(5) are IR divergent:
      d15(0)=D0_ir(xmt,xpi1,q13,q23,xpi2,xpi3,q12,mi1,mi3,mi2)
      call ddfunc_ir(xmt,xpi1,xpi2,xpi3,pi12,pi13,pi23,
     $     m0,mi1,mi2,mi3,d15,d25,d35,c1,c2)
c      write(6,*)'d0(1)=',d11(0)
c      write(6,*)'d0(2)=',d12(0)
c      write(6,*)'d0(3)=',d13(0)
c      write(6,*)'d0(4)=',d14(0)
c      write(6,*)'d0(5)=',d15(0)
c
*--entwicklungskoeff. des vektorintegrals:----------------------------
      s11 = (d12(0)-d11(0))/2d0
      s12 = (d13(0)-d11(0))/2d0
      s13 = (d14(0)-d11(0)-cf3*e1(0))/2d0
      s14 = (d15(0)-d11(0)-cf4*e1(0))/2d0
      e1(1) = (mat11*s11+mat12*s12+mat13*s13+mat14*s14)/det
      e1(2) = (mat21*s11+mat22*s12+mat23*s13+mat24*s14)/det
      e1(3) = (mat31*s11+mat32*s12+mat33*s13+mat34*s14)/det
      e1(4) = (mat41*s11+mat42*s12+mat43*s13+mat44*s14)/det
*---------------------------------------------------------------------
*--entwicklungskoeff. des tensorintegrals (Rank 2):-------------------
      do i=0,3
         d11(i)=d11(i)*dets
         d12(i)=d12(i)*dets
         d13(i)=d13(i)*dets
         d14(i)=d14(i)*dets
         d15(i)=d15(i)*dets
      enddo
      ddd = d11(1)+d11(2)+d11(3)+d11(0)
      s211 = ddd/2d0
      s212 = (d12(1)-d11(1))/2d0
      s213 = (d12(2)-d11(2))/2d0
      s214 = (d12(3)-d11(3))/2d0
      s221 = (d13(1)+ddd)/2d0
      s222 = -d11(1)/2d0
      s223 = (d13(2)-d11(2))/2d0
      s224 = (d13(3)-d11(3))/2d0
      s231 = (d14(1)+ddd-cf3*e1(1))/2d0
      s232 = (d14(2)-d11(1)-cf3*e1(2))/2d0
      s233 = -(d11(2)+cf3*e1(3))/2d0
      s234 = (d14(3)-d11(3)-cf3*e1(4))/2d0
      s241 = (d15(1)+ddd-cf4*e1(1))/2d0
      s242 = (d15(2)-d11(1)-cf4*e1(2))/2d0
      s243 = (d15(3)-d11(2)-cf4*e1(3))/2d0
      s244 = -(d11(3)+cf4*e1(4))/2d0
      e2(1,1) = (mat11*s211+mat12*s221+mat13*s231+mat14*s241)/det
      e2(2,2) = (mat21*s212+mat22*s222+mat23*s232+mat24*s242)/det
      e2(3,3) = (mat31*s213+mat32*s223+mat33*s233+mat34*s243)/det
      e2(4,4) = (mat41*s214+mat42*s224+mat43*s234+mat44*s244)/det
      e2(1,2) = (mat11*s212+mat12*s222+mat13*s232+mat14*s242)/det
      e2(1,3) = (mat11*s213+mat12*s223+mat13*s233+mat14*s243)/det
      e2(1,4) = (mat11*s214+mat12*s224+mat13*s234+mat14*s244)/det
      e2(2,1) = (mat21*s211+mat22*s221+mat23*s231+mat24*s241)/det
      e2(2,3) = (mat21*s213+mat22*s223+mat23*s233+mat24*s243)/det
      e2(2,4) = (mat21*s214+mat22*s224+mat23*s234+mat24*s244)/det
      e2(3,1) = (mat31*s211+mat32*s221+mat33*s231+mat34*s241)/det
      e2(3,2) = (mat31*s212+mat32*s222+mat33*s232+mat34*s242)/det
      e2(3,4) = (mat31*s214+mat32*s224+mat33*s234+mat34*s244)/det
      e2(4,1) = (mat41*s211+mat42*s221+mat43*s231+mat44*s241)/det
      e2(4,2) = (mat41*s212+mat42*s222+mat43*s232+mat44*s242)/det
      e2(4,3) = (mat41*s213+mat42*s223+mat43*s233+mat44*s243)/det
*---------------------------------------------------------------------
*--entwicklungskoeff. des tensorintegrals (Rank 3):-------------------
      diff(1)=d22(0)-d21(0)
      diff(2)=d23(0)-d21(0)
      diff(3)=d24(0)-d21(0)
      diff(4)=d25(0)-d21(0)
      do i=1,4
         diff(i)=diff(i)/det*dets
      enddo
      do i=0,3
         d11(i)=d11(i)*dets
         d12(i)=d12(i)*dets
         d13(i)=d13(i)*dets
         d14(i)=d14(i)*dets
         d15(i)=d15(i)*dets
      enddo
      do i=0,6
         d21(i)=d21(i)*dets**2
         d22(i)=d22(i)*dets**2
         d23(i)=d23(i)*dets**2
         d24(i)=d24(i)*dets**2
         d25(i)=d25(i)*dets**2
      enddo
      ddd = ddd*dets+
     $     (d21(1)+d21(2)+d21(3)-d11(0))/2d0+d21(4)+d21(5)+d21(6)
      s311 = -ddd+mat11*diff(1)/2d0
      s321 = (d22(1)-d21(1))/2d0+mat22*diff(1)/2d0
      s331 = (d22(2)-d21(2))/2d0+mat33*diff(1)/2d0
      s341 = (d22(3)-d21(3))/2d0+mat44*diff(1)/2d0
      s312 = -ddd+d23(1)/2d0+mat11*diff(2)/2d0
      s322 = -d21(1)/2d0+mat22*diff(2)/2d0
      s332 = (d23(2)-d21(2))/2d0+mat33*diff(2)/2d0
      s342 = (d23(3)-d21(3))/2d0+mat44*diff(2)/2d0
      s313 = -ddd+(d24(1)-cf3*e2(1,1))/2d0+mat11*diff(3)/2d0
      s323 = (d24(2)-d21(1)-cf3*e2(2,2))/2d0+mat22*diff(3)/2d0
      s333 = (-d21(2)-cf3*e2(3,3))/2d0+mat33*diff(3)/2d0
      s343 = (d24(3)-d21(3)-cf3*e2(4,4))/2d0+mat44*diff(3)/2d0
      s314 = -ddd+(d25(1)-cf4*e2(1,1))/2d0+mat11*diff(4)/2d0
      s324 = (d25(2)-d21(1)-cf4*e2(2,2))/2d0+mat22*diff(4)/2d0
      s334 = (d25(3)-d21(2)-cf4*e2(3,3))/2d0+mat33*diff(4)/2d0
      s344 = (-d21(3)-cf4*e2(4,4))/2d0+mat44*diff(4)/2d0
      s3121 = (d11(1)+d21(1)+d21(4)+d21(5))/2d0+mat12*diff(1)/2d0
      s3122 = (d11(1)+d21(1)+d21(4)+d21(5))/2d0+mat12*diff(2)/2d0
      s3123 = (d11(1)+d21(1)+d21(4)+d24(4)+d21(5)-cf3*e2(1,2))/2d0+
     $     mat12*diff(3)/2d0
      s3124 = (d11(1)+d21(1)+d21(4)+d25(4)+d21(5)-cf4*e2(1,2))/2d0+
     $     mat12*diff(4)/2d0
      s3341 = (d22(6)-d21(6))/2d0+mat34*diff(1)/2d0
      s3342 = (d23(6)-d21(6))/2d0+mat34*diff(2)/2d0
      s3343 = (-d21(6)-cf3*e2(3,4))/2d0+mat34*diff(3)/2d0
      s3344 = (-d21(6)-cf4*e2(3,4))/2d0+mat34*diff(4)/2d0
c e3123
      e3(0,1) = (mat31*s3121+mat32*s3122+mat33*s3123+mat34*s3124)/det
c e3124
      e3(0,2) = (mat41*s3121+mat42*s3122+mat43*s3123+mat44*s3124)/det
c e3134
      e3(0,3) = (mat11*s3341+mat12*s3342+mat13*s3343+mat14*s3344)/det
c e3234
      e3(0,4) = (mat21*s3341+mat22*s3342+mat23*s3343+mat24*s3344)/det
c
      e3(1,1) = (mat11*s311+mat12*s312+mat13*s313+mat14*s314)/det
      e3(2,2) = (mat21*s321+mat22*s322+mat23*s323+mat24*s324)/det
      e3(3,3) = (mat31*s331+mat32*s332+mat33*s333+mat34*s334)/det
      e3(4,4) = (mat41*s341+mat42*s342+mat43*s343+mat44*s344)/det
      e3(1,2) = (mat21*s311+mat22*s312+mat23*s313+mat24*s314)/det
      e3(1,3) = (mat31*s311+mat32*s312+mat33*s313+mat34*s314)/det
      e3(1,4) = (mat41*s311+mat42*s312+mat43*s313+mat44*s314)/det
      e3(2,1) = (mat11*s321+mat12*s322+mat13*s323+mat14*s324)/det
      e3(2,3) = (mat31*s321+mat32*s322+mat33*s323+mat34*s324)/det
      e3(2,4) = (mat41*s321+mat42*s322+mat43*s323+mat44*s324)/det
      e3(3,1) = (mat11*s331+mat12*s332+mat13*s333+mat14*s334)/det
      e3(3,2) = (mat21*s331+mat22*s332+mat23*s333+mat24*s334)/det
      e3(3,4) = (mat41*s331+mat42*s332+mat43*s333+mat44*s334)/det
      e3(4,1) = (mat11*s341+mat12*s342+mat13*s343+mat14*s344)/det
      e3(4,2) = (mat21*s341+mat22*s342+mat23*s343+mat24*s344)/det
      e3(4,3) = (mat31*s341+mat32*s342+mat33*s343+mat34*s344)/det
c
      if(nlim.gt.0)then
         do i=0,3
            d1(i)=d1(i)*dets**3
         enddo
         do i=0,6
            d2(i)=d2(i)*dets**3
         enddo
         e1(0)=e1(0)*dets**3
         do i=1,4
            e1(i)=e1(i)*dets**2
         enddo
         do i=1,4
            do j=1,4
               e2(i,j)=e2(i,j)*dets
            enddo
         enddo
      endif
c      do i=1,4
c         e1(i)=dcmplx(0d0,0d0)
c      enddo
c      do i=1,4
c         do j=1,4
c            e2(i,j)=dcmplx(0d0,0d0)
c         enddo
c      enddo
c      do i=0,4
c         do j=1,4
c            e3(i,j)=dcmplx(0d0,0d0)
c         enddo
c      enddo
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine ddpenta_ir(np,xmt,xpi1,xpi2,xpi3,xpi4,pi12,pi13,
     $     pi14,pi23,pi24,pi34,m0,mi1,mi2,mi3,mi4,d1,d2,d3)
*---------------------------------------------------------------------
c entwicklungskoeff. fuer das vektorielle und tensorielle fuenfpkt.-integral
c d0=k^2-m0^2, d1=(k+pi1)^2-m1^2, d2=(k+pi2)^2-m2^2, d3=(k+pi3)^2-m3^2,
c d4=(k+pi4)^2-m4^2
c xpi1=pi1^2,xpi2=pi2^2,xpi3=pi3^2,xpi4=pi4^2,pi12=pi1.pi2,pi13=pi1.pi3,
c pi14=pi1.pi4,pi23=pi2.pi3,pi24=pi2.pi4,pi34=pi3.pi4
c Up to rank 3, only IR div. pentagon diagrams !
c always m0=0 and xpi1=mi1^2 and xpi2=mi2^2 and mi3=mt and mi4=mt
c e1(0) is IR div. and input ! 
c D.W., June 18, 2002
*---------------------------------------------------------------------
      implicit none
      integer i,j,np
      real*8 xpi1,xpi2,xpi3,xpi4,pi12,pi13,pi14,pi23,pi24,pi34,
     $     q12,q13,q14,q23,q24,q34,mi1,mi2,mi3,mi4,m0,
     $     xmi1,xmi2,xmi3,xmi4,xm0,xmt
      complex*16 d11(0:3),d21(0:6),d31(0:3,0:3),D0_ir
      complex*16 d12(0:3),d22(0:6),d32(0:3,0:3)
      complex*16 d13(0:3),d23(0:6),d33(0:3,0:3)
      complex*16 d14(0:3),d24(0:6),d34(0:3,0:3)
      complex*16 d15(0:3),d25(0:6),d35(0:3,0:3)
      complex*16 d1(0:3,5),d2(0:6,5),d3(0:3,0:3,5)
      complex*16 c1(0:2,4),c2(0:3,4)
      xm0 = m0*m0
      xmi1 = mi1*mi1
      xmi2 = mi2*mi2
      xmi3 = mi3*mi3
      xmi4 = mi4*mi4
      q12=xpi1+xpi2-2d0*pi12
      q13=xpi1+xpi3-2d0*pi13
      q14=xpi1+xpi4-2d0*pi14
      q23=xpi2+xpi3-2d0*pi23
      q24=xpi2+xpi4-2d0*pi24
      q34=xpi3+xpi4-2d0*pi34
*
      call ddfunc(q12,q13,q14,pi23-pi12-pi13+xpi1,
     $     pi24-pi12-pi14+xpi1,xpi1+pi34-pi14-pi13,
     $     mi1,mi2,mi3,mi4,d11,d21,d31,c1,c2)
c only D0(2) of P1 is IR divergent:
      if(np.eq.1)then
         d12(0)=D0_ir(xmt,xpi2,xpi3,q34,q24,q23,xpi4,mi2,mi3,mi4)
         call ddfunc_ir(xmt,xpi2,xpi3,xpi4,pi23,pi24,pi34,
     $        m0,mi2,mi3,mi4,d12,d22,d32,c1,c2)
      else
         call ddfunc(xpi2,xpi3,xpi4,pi23,pi24,pi34,
     $        m0,mi2,mi3,mi4,d12,d22,d32,c1,c2)
      endif
c only D0(3) of P1,P5,P6 are IR divergent:
      if(np.eq.1.or.np.eq.5.or.np.eq.6)then
         d13(0)=D0_ir(xmt,xpi1,xpi4,q34,q13,q14,xpi3,mi1,mi4,mi3)
         call ddfunc_ir(xmt,xpi1,xpi3,xpi4,pi13,pi14,pi34,
     $        m0,mi1,mi3,mi4,d13,d23,d33,c1,c2)
      else
         call ddfunc(xpi1,xpi3,xpi4,pi13,pi14,pi34,
     $        m0,mi1,mi3,mi4,d13,d23,d33,c1,c2)
      endif
c all D0(4) are IR divergent:
      d14(0)=D0_ir(xmt,xpi1,q14,q24,xpi2,xpi4,q12,mi1,mi4,mi2)
      call ddfunc_ir(xmt,xpi1,xpi2,xpi4,pi12,pi14,pi24,
     $     m0,mi1,mi2,mi4,d14,d24,d34,c1,c2)
c all D0(5) are IR divergent:
      d15(0)=D0_ir(xmt,xpi1,q13,q23,xpi2,xpi3,q12,mi1,mi3,mi2)
      call ddfunc_ir(xmt,xpi1,xpi2,xpi3,pi12,pi13,pi23,
     $     m0,mi1,mi2,mi3,d15,d25,d35,c1,c2)
      do i=0,3
         d1(i,1)=d11(i)
         d1(i,2)=d12(i)
         d1(i,3)=d13(i)
         d1(i,4)=d14(i)
         d1(i,5)=d15(i)
      enddo
      do i=0,6
         d2(i,1)=d21(i)
         d2(i,2)=d22(i)
         d2(i,3)=d23(i)
         d2(i,4)=d24(i)
         d2(i,5)=d25(i)
      enddo
      do i=0,3
         do j=0,3
            d3(i,j,1)=d31(i,j)
            d3(i,j,2)=d32(i,j)
            d3(i,j,3)=d33(i,j)
            d3(i,j,4)=d34(i,j)
            d3(i,j,5)=d35(i,j)
         enddo
      enddo
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      function E0_ir(np,xmt,p1,p2,p3,p4,p5,p12,p23,p34,p45,p15,
     $     m0,m1,m2,m3,m4)
*********************************************************************
* IR singular scalar 5-point functions for the six Pentagons
* Dimensional regularization in d=4-2 eps_IR
* G[1+eps] (4 pi mu^2/mt^2)^eps_IR x poles are subtracted
*--------------------------------------------------------------------
* d0=k^2-m0^2, d1=(k+p1)^2-m1^2, d2=(k+p1+p2)^2-m2^2, 
* d3=(k+p3+p4)^2-m3^2, d4=(k-p4)^2-m4^2
* p12=(p1+p2)^2,p23=(p2+p3)^2,p34=(p3+p4)^2,p45=(p4+p5)^2,
* p15=(p1+p5)^2
* always m0=0 and p1^2=m1^2 and p4^2=m4^2. 
*--------------------------------------------------------------------
* June 19, 2002, DW
*********************************************************************
      implicit none
      integer i,np
      real*8 xmt,p1,p2,p3,p4,p5,p12,p23,p34,p45,p15
      real*8 m0,m1,m2,m3,m4,m02,m12,m22,m32,m42
      real*8 detyy,y(0:4,0:4),dety(0:4)
      complex*16 E0_ir,D0_,D0_ir,d0(0:4)

      m02 = m0**2
      m12 = m1**2
      m22 = m2**2
      m32 = m3**2
      m42 = m4**2
      
      d0(0)=D0_(p2,p3,p4,p15,p23,p34,m1,m2,m3,m4,0)
      if (np.eq.1) then
         d0(1)=D0_ir(xmt,p5,p12,p3,p4,p34,p45,m4,m2,m3)
      else
         d0(1)=D0_(p12,p3,p4,p5,p45,p34,m0,m2,m3,m4,0)
      endif
      d0(2)=D0_ir(xmt,p1,p23,p4,p5,p45,p15,m1,m3,m4)
      d0(3)=D0_ir(xmt,p1,p2,p34,p5,p12,p15,m1,m2,m4)
      if (np.eq.1.or.np.eq.5.or.np.eq.6) then
         d0(4)=D0_ir(xmt,p1,p45,p3,p2,p23,p12,m1,m3,m2)
      elseif(np.eq.2.or.np.eq.4)then
         d0(4)=D0_(p1,p2,p3,p45,p12,p23,m0,m1,m2,m3,0)
      elseif(np.eq.3)then
         d0(4)=D0_(p1,p4,p3,p45,p12,p23,m0,m1,m2,m3,0)
      endif
c      write(6,*)'d0(0)=',d0(0)
c      write(6,*)'d0(1)=',d0(1)
c      write(6,*)'d0(2)=',d0(2)
c      write(6,*)'d0(3)=',d0(3)
c      write(6,*)'d0(4)=',d0(4)
c
      y(0,0) = 2d0*m02
      y(1,1) = 2d0*m12
      y(2,2) = 2d0*m22
      y(3,3) = 2d0*m32
      y(4,4) = 2d0*m42
      y(0,1) = m02+m12-p1
      y(1,2) = m12+m22-p2
      y(2,3) = m22+m32-p3
      y(3,4) = m32+m42-p4
      y(0,4) = m42+m02-p5
      y(0,2) = m02+m22-p12
      y(1,3) = m12+m32-p23
      y(2,4) = m22+m42-p34
      y(0,3) = m32+m02-p45
      y(1,4) = m42+m12-p15

      DETYY = Y(0,4)**2*Y(1,3)**2*Y(2,2)-2D0*Y(0,3)*Y(0,4)*Y(1,3)*
     &     Y(1,4)*Y(2,2)+Y(0,3)**2*Y(1,4)**2*Y(2,2)-2D0*Y(0,4)**2*Y(1,
     &     2)*Y(1,3)*Y(2,3)+2D0*Y(0,3)*Y(0,4)*Y(1,2)*Y(1,4)*Y(2,3)+2D0*
     &     Y(0,2)*Y(0,4)*Y(1,3)*Y(1,4)*Y(2,3)-2D0*Y(0,2)*Y(0,3)*Y(1,
     &     4)**2*Y(2,3)+Y(0,4)**2*Y(1,1)*Y(2,3)**2-2D0*Y(0,1)*Y(0,4)*
     &     Y(1,4)*Y(2,3)**2+Y(0,0)*Y(1,4)**2*Y(2,3)**2+2D0*Y(0,3)*Y(0,
     &     4)*Y(1,2)*Y(1,3)*Y(2,4)-2D0*Y(0,2)*Y(0,4)*Y(1,3)**2*Y(2,4)-
     &     2D0*Y(0,3)**2*Y(1,2)*Y(1,4)*Y(2,4)+2D0*Y(0,2)*Y(0,3)*Y(1,3)*
     &     Y(1,4)*Y(2,4)-2D0*Y(0,3)*Y(0,4)*Y(1,1)*Y(2,3)*Y(2,4)+2D0*
     &     Y(0,1)*Y(0,4)*Y(1,3)*Y(2,3)*Y(2,4)+2D0*Y(0,1)*Y(0,3)*Y(1,4)*
     &     Y(2,3)*Y(2,4)-2D0*Y(0,0)*Y(1,3)*Y(1,4)*Y(2,3)*Y(2,4)+Y(0,
     &     3)**2*Y(1,1)*Y(2,4)**2-2D0*Y(0,1)*Y(0,3)*Y(1,3)*Y(2,4)**2+
     &     Y(0,0)*Y(1,3)**2*Y(2,4)**2+Y(0,4)**2*Y(1,2)**2*Y(3,3)-2D0*
     &     Y(0,2)*Y(0,4)*Y(1,2)*Y(1,4)*Y(3,3)+Y(0,2)**2*Y(1,4)**2*Y(3,
     &     3)-Y(0,4)**2*Y(1,1)*Y(2,2)*Y(3,3)+2D0*Y(0,1)*Y(0,4)*Y(1,4)*
     &     Y(2,2)*Y(3,3)-Y(0,0)*Y(1,4)**2*Y(2,2)*Y(3,3)+2D0*Y(0,2)*Y(0,
     &     4)*Y(1,1)*Y(2,4)*Y(3,3)-2D0*Y(0,1)*Y(0,4)*Y(1,2)*Y(2,4)*Y(3,
     &     3)-2D0*Y(0,1)*Y(0,2)*Y(1,4)*Y(2,4)*Y(3,3)+2D0*Y(0,0)*Y(1,2)*
     &     Y(1,4)*Y(2,4)*Y(3,3)+Y(0,1)**2*Y(2,4)**2*Y(3,3)-Y(0,0)*Y(1,
     &     1)*Y(2,4)**2*Y(3,3)-2D0*Y(0,3)*Y(0,4)*Y(1,2)**2*Y(3,4)+2D0*
     &     Y(0,2)*Y(0,4)*Y(1,2)*Y(1,3)*Y(3,4)+2D0*Y(0,2)*Y(0,3)*Y(1,2)*
     &     Y(1,4)*Y(3,4)-2D0*Y(0,2)**2*Y(1,3)*Y(1,4)*Y(3,4)+2D0*Y(0,3)*
     &     Y(0,4)*Y(1,1)*Y(2,2)*Y(3,4)-2D0*Y(0,1)*Y(0,4)*Y(1,3)*Y(2,2)*
     &     Y(3,4)-2D0*Y(0,1)*Y(0,3)*Y(1,4)*Y(2,2)*Y(3,4)+2D0*Y(0,0)*
     &     Y(1,3)*Y(1,4)*Y(2,2)*Y(3,4)-2D0*Y(0,2)*Y(0,4)*Y(1,1)*Y(2,3)*
     &     Y(3,4)+2D0*Y(0,1)*Y(0,4)*Y(1,2)*Y(2,3)*Y(3,4)+2D0*Y(0,1)*
     &     Y(0,2)*Y(1,4)*Y(2,3)*Y(3,4)-2D0*Y(0,0)*Y(1,2)*Y(1,4)*Y(2,3)*
     &     Y(3,4)-2D0*Y(0,2)*Y(0,3)*Y(1,1)*Y(2,4)*Y(3,4)+2D0*Y(0,1)*
     &     Y(0,3)*Y(1,2)*Y(2,4)*Y(3,4)+2D0*Y(0,1)*Y(0,2)*Y(1,3)*Y(2,4)*
     &     Y(3,4)-2D0*Y(0,0)*Y(1,2)*Y(1,3)*Y(2,4)*Y(3,4)-2D0*Y(0,1)**2*
     &     Y(2,3)*Y(2,4)*Y(3,4)+2D0*Y(0,0)*Y(1,1)*Y(2,3)*Y(2,4)*Y(3,4)+
     &     Y(0,2)**2*Y(1,1)*Y(3,4)**2-2D0*Y(0,1)*Y(0,2)*Y(1,2)*Y(3,4)**
     &     2+Y(0,0)*Y(1,2)**2*Y(3,4)**2+Y(0,1)**2*Y(2,2)*Y(3,4)**2-Y(0,
     &     0)*Y(1,1)*Y(2,2)*Y(3,4)**2+Y(0,3)**2*Y(1,2)**2*Y(4,4)-2D0*
     &     Y(0,2)*Y(0,3)*Y(1,2)*Y(1,3)*Y(4,4)+Y(0,2)**2*Y(1,3)**2*Y(4,
     &     4)-Y(0,3)**2*Y(1,1)*Y(2,2)*Y(4,4)+2D0*Y(0,1)*Y(0,3)*Y(1,3)*
     &     Y(2,2)*Y(4,4)-Y(0,0)*Y(1,3)**2*Y(2,2)*Y(4,4)+2D0*Y(0,2)*Y(0,
     &     3)*Y(1,1)*Y(2,3)*Y(4,4)-2D0*Y(0,1)*Y(0,3)*Y(1,2)*Y(2,3)*Y(4,
     &     4)-2D0*Y(0,1)*Y(0,2)*Y(1,3)*Y(2,3)*Y(4,4)+2D0*Y(0,0)*Y(1,2)*
     &     Y(1,3)*Y(2,3)*Y(4,4)+Y(0,1)**2*Y(2,3)**2*Y(4,4)-Y(0,0)*Y(1,
     &     1)*Y(2,3)**2*Y(4,4)-Y(0,2)**2*Y(1,1)*Y(3,3)*Y(4,4)+2D0*Y(0,
     &     1)*Y(0,2)*Y(1,2)*Y(3,3)*Y(4,4)-Y(0,0)*Y(1,2)**2*Y(3,3)*Y(4,
     &     4)-Y(0,1)**2*Y(2,2)*Y(3,3)*Y(4,4)+Y(0,0)*Y(1,1)*Y(2,2)*Y(3,
     &     3)*Y(4,4)

c      if (dabs(detyy).lt.1d-1) then
c         write(6,*)'warning: potential instability in E0_ir'
c         stop
c      endif

      DETY(0) = -(Y(0,4)*Y(1,3)**2*Y(2,2))+Y(0,3)*Y(1,3)*Y(1,4)*
     &     Y(2,2)+Y(0,4)*Y(1,3)*Y(1,4)*Y(2,2)-Y(0,3)*Y(1,4)**2*Y(2,2)+
     &     2D0*Y(0,4)*Y(1,2)*Y(1,3)*Y(2,3)-Y(0,3)*Y(1,2)*Y(1,4)*Y(2,3)-
     &     Y(0,4)*Y(1,2)*Y(1,4)*Y(2,3)-Y(0,2)*Y(1,3)*Y(1,4)*Y(2,3)-Y(0,
     &     4)*Y(1,3)*Y(1,4)*Y(2,3)+Y(0,2)*Y(1,4)**2*Y(2,3)+Y(0,3)*Y(1,
     &     4)**2*Y(2,3)-Y(0,4)*Y(1,1)*Y(2,3)**2+Y(0,1)*Y(1,4)*Y(2,3)**
     &     2+Y(0,4)*Y(1,4)*Y(2,3)**2-Y(1,4)**2*Y(2,3)**2-Y(0,3)*Y(1,2)*
     &     Y(1,3)*Y(2,4)-Y(0,4)*Y(1,2)*Y(1,3)*Y(2,4)+Y(0,2)*Y(1,3)**2*
     &     Y(2,4)+Y(0,4)*Y(1,3)**2*Y(2,4)+2D0*Y(0,3)*Y(1,2)*Y(1,4)*Y(2,
     &     4)-Y(0,2)*Y(1,3)*Y(1,4)*Y(2,4)-Y(0,3)*Y(1,3)*Y(1,4)*Y(2,4)+
     &     Y(0,3)*Y(1,1)*Y(2,3)*Y(2,4)+Y(0,4)*Y(1,1)*Y(2,3)*Y(2,4)-Y(0,
     &     1)*Y(1,3)*Y(2,3)*Y(2,4)-Y(0,4)*Y(1,3)*Y(2,3)*Y(2,4)-Y(0,1)*
     &     Y(1,4)*Y(2,3)*Y(2,4)-Y(0,3)*Y(1,4)*Y(2,3)*Y(2,4)+2D0*Y(1,3)*
     &     Y(1,4)*Y(2,3)*Y(2,4)-Y(0,3)*Y(1,1)*Y(2,4)**2+Y(0,1)*Y(1,3)*
     &     Y(2,4)**2+Y(0,3)*Y(1,3)*Y(2,4)**2-Y(1,3)**2*Y(2,4)**2-Y(0,
     &     4)*Y(1,2)**2*Y(3,3)+Y(0,2)*Y(1,2)*Y(1,4)*Y(3,3)+Y(0,4)*Y(1,
     &     2)*Y(1,4)*Y(3,3)-Y(0,2)*Y(1,4)**2*Y(3,3)+Y(0,4)*Y(1,1)*Y(2,
     &     2)*Y(3,3)-Y(0,1)*Y(1,4)*Y(2,2)*Y(3,3)-Y(0,4)*Y(1,4)*Y(2,2)*
     &     Y(3,3)+Y(1,4)**2*Y(2,2)*Y(3,3)-Y(0,2)*Y(1,1)*Y(2,4)*Y(3,3)-
     &     Y(0,4)*Y(1,1)*Y(2,4)*Y(3,3)+Y(0,1)*Y(1,2)*Y(2,4)*Y(3,3)+Y(0,
     &     4)*Y(1,2)*Y(2,4)*Y(3,3)+Y(0,1)*Y(1,4)*Y(2,4)*Y(3,3)+Y(0,2)*
     &     Y(1,4)*Y(2,4)*Y(3,3)-2D0*Y(1,2)*Y(1,4)*Y(2,4)*Y(3,3)-Y(0,1)*
     &     Y(2,4)**2*Y(3,3)+Y(1,1)*Y(2,4)**2*Y(3,3)+Y(0,3)*Y(1,2)**2*
     &     Y(3,4)+Y(0,4)*Y(1,2)**2*Y(3,4)-Y(0,2)*Y(1,2)*Y(1,3)*Y(3,4)-
     &     Y(0,4)*Y(1,2)*Y(1,3)*Y(3,4)-Y(0,2)*Y(1,2)*Y(1,4)*Y(3,4)-Y(0,
     &     3)*Y(1,2)*Y(1,4)*Y(3,4)+2D0*Y(0,2)*Y(1,3)*Y(1,4)*Y(3,4)-Y(0,
     &     3)*Y(1,1)*Y(2,2)*Y(3,4)-Y(0,4)*Y(1,1)*Y(2,2)*Y(3,4)+Y(0,1)*
     &     Y(1,3)*Y(2,2)*Y(3,4)+Y(0,4)*Y(1,3)*Y(2,2)*Y(3,4)+Y(0,1)*Y(1,
     &     4)*Y(2,2)*Y(3,4)+Y(0,3)*Y(1,4)*Y(2,2)*Y(3,4)-2D0*Y(1,3)*Y(1,
     &     4)*Y(2,2)*Y(3,4)+Y(0,2)*Y(1,1)*Y(2,3)*Y(3,4)+Y(0,4)*Y(1,1)*
     &     Y(2,3)*Y(3,4)-Y(0,1)*Y(1,2)*Y(2,3)*Y(3,4)-Y(0,4)*Y(1,2)*Y(2,
     &     3)*Y(3,4)-Y(0,1)*Y(1,4)*Y(2,3)*Y(3,4)-Y(0,2)*Y(1,4)*Y(2,3)*
     &     Y(3,4)+2D0*Y(1,2)*Y(1,4)*Y(2,3)*Y(3,4)+Y(0,2)*Y(1,1)*Y(2,4)*
     &     Y(3,4)+Y(0,3)*Y(1,1)*Y(2,4)*Y(3,4)-Y(0,1)*Y(1,2)*Y(2,4)*Y(3,
     &     4)-Y(0,3)*Y(1,2)*Y(2,4)*Y(3,4)-Y(0,1)*Y(1,3)*Y(2,4)*Y(3,4)-
     &     Y(0,2)*Y(1,3)*Y(2,4)*Y(3,4)+2D0*Y(1,2)*Y(1,3)*Y(2,4)*Y(3,4)+
     &     2D0*Y(0,1)*Y(2,3)*Y(2,4)*Y(3,4)-2D0*Y(1,1)*Y(2,3)*Y(2,4)*
     &     Y(3,4)-Y(0,2)*Y(1,1)*Y(3,4)**2+Y(0,1)*Y(1,2)*Y(3,4)**2+Y(0,
     &     2)*Y(1,2)*Y(3,4)**2-Y(1,2)**2*Y(3,4)**2-Y(0,1)*Y(2,2)*Y(3,
     &     4)**2+Y(1,1)*Y(2,2)*Y(3,4)**2-Y(0,3)*Y(1,2)**2*Y(4,4)+Y(0,
     &     2)*Y(1,2)*Y(1,3)*Y(4,4)+Y(0,3)*Y(1,2)*Y(1,3)*Y(4,4)-Y(0,2)*
     &     Y(1,3)**2*Y(4,4)+Y(0,3)*Y(1,1)*Y(2,2)*Y(4,4)-Y(0,1)*Y(1,3)*
     &     Y(2,2)*Y(4,4)-Y(0,3)*Y(1,3)*Y(2,2)*Y(4,4)+Y(1,3)**2*Y(2,2)*
     &     Y(4,4)-Y(0,2)*Y(1,1)*Y(2,3)*Y(4,4)-Y(0,3)*Y(1,1)*Y(2,3)*Y(4,
     &     4)+Y(0,1)*Y(1,2)*Y(2,3)*Y(4,4)+Y(0,3)*Y(1,2)*Y(2,3)*Y(4,4)+
     &     Y(0,1)*Y(1,3)*Y(2,3)*Y(4,4)+Y(0,2)*Y(1,3)*Y(2,3)*Y(4,4)-2D0*
     &     Y(1,2)*Y(1,3)*Y(2,3)*Y(4,4)-Y(0,1)*Y(2,3)**2*Y(4,4)+Y(1,1)*
     &     Y(2,3)**2*Y(4,4)+Y(0,2)*Y(1,1)*Y(3,3)*Y(4,4)-Y(0,1)*Y(1,2)*
     &     Y(3,3)*Y(4,4)-Y(0,2)*Y(1,2)*Y(3,3)*Y(4,4)+Y(1,2)**2*Y(3,3)*
     &     Y(4,4)+Y(0,1)*Y(2,2)*Y(3,3)*Y(4,4)-Y(1,1)*Y(2,2)*Y(3,3)*Y(4,
     &     4)
      DETY(1) = Y(0,3)*Y(0,4)*Y(1,3)*Y(2,2)-Y(0,4)**2*Y(1,3)*Y(2,
     &     2)-Y(0,3)**2*Y(1,4)*Y(2,2)+Y(0,3)*Y(0,4)*Y(1,4)*Y(2,2)-Y(0,
     &     3)*Y(0,4)*Y(1,2)*Y(2,3)+Y(0,4)**2*Y(1,2)*Y(2,3)-Y(0,2)*Y(0,
     &     4)*Y(1,3)*Y(2,3)+Y(0,4)**2*Y(1,3)*Y(2,3)+2D0*Y(0,2)*Y(0,3)*
     &     Y(1,4)*Y(2,3)-Y(0,2)*Y(0,4)*Y(1,4)*Y(2,3)-Y(0,3)*Y(0,4)*Y(1,
     &     4)*Y(2,3)+Y(0,1)*Y(0,4)*Y(2,3)**2-Y(0,4)**2*Y(2,3)**2-Y(0,
     &     0)*Y(1,4)*Y(2,3)**2+Y(0,4)*Y(1,4)*Y(2,3)**2+Y(0,3)**2*Y(1,
     &     2)*Y(2,4)-Y(0,3)*Y(0,4)*Y(1,2)*Y(2,4)-Y(0,2)*Y(0,3)*Y(1,3)*
     &     Y(2,4)+2D0*Y(0,2)*Y(0,4)*Y(1,3)*Y(2,4)-Y(0,3)*Y(0,4)*Y(1,3)*
     &     Y(2,4)-Y(0,2)*Y(0,3)*Y(1,4)*Y(2,4)+Y(0,3)**2*Y(1,4)*Y(2,4)-
     &     Y(0,1)*Y(0,3)*Y(2,3)*Y(2,4)-Y(0,1)*Y(0,4)*Y(2,3)*Y(2,4)+2D0*
     &     Y(0,3)*Y(0,4)*Y(2,3)*Y(2,4)+Y(0,0)*Y(1,3)*Y(2,3)*Y(2,4)-Y(0,
     &     4)*Y(1,3)*Y(2,3)*Y(2,4)+Y(0,0)*Y(1,4)*Y(2,3)*Y(2,4)-Y(0,3)*
     &     Y(1,4)*Y(2,3)*Y(2,4)+Y(0,1)*Y(0,3)*Y(2,4)**2-Y(0,3)**2*Y(2,
     &     4)**2-Y(0,0)*Y(1,3)*Y(2,4)**2+Y(0,3)*Y(1,3)*Y(2,4)**2+Y(0,
     &     2)*Y(0,4)*Y(1,2)*Y(3,3)-Y(0,4)**2*Y(1,2)*Y(3,3)-Y(0,2)**2*
     &     Y(1,4)*Y(3,3)+Y(0,2)*Y(0,4)*Y(1,4)*Y(3,3)-Y(0,1)*Y(0,4)*Y(2,
     &     2)*Y(3,3)+Y(0,4)**2*Y(2,2)*Y(3,3)+Y(0,0)*Y(1,4)*Y(2,2)*Y(3,
     &     3)-Y(0,4)*Y(1,4)*Y(2,2)*Y(3,3)+Y(0,1)*Y(0,2)*Y(2,4)*Y(3,3)+
     &     Y(0,1)*Y(0,4)*Y(2,4)*Y(3,3)-2D0*Y(0,2)*Y(0,4)*Y(2,4)*Y(3,3)-
     &     Y(0,0)*Y(1,2)*Y(2,4)*Y(3,3)+Y(0,4)*Y(1,2)*Y(2,4)*Y(3,3)-Y(0,
     &     0)*Y(1,4)*Y(2,4)*Y(3,3)+Y(0,2)*Y(1,4)*Y(2,4)*Y(3,3)+Y(0,0)*
     &     Y(2,4)**2*Y(3,3)-Y(0,1)*Y(2,4)**2*Y(3,3)-Y(0,2)*Y(0,3)*Y(1,
     &     2)*Y(3,4)-Y(0,2)*Y(0,4)*Y(1,2)*Y(3,4)+2D0*Y(0,3)*Y(0,4)*Y(1,
     &     2)*Y(3,4)+Y(0,2)**2*Y(1,3)*Y(3,4)-Y(0,2)*Y(0,4)*Y(1,3)*Y(3,
     &     4)+Y(0,2)**2*Y(1,4)*Y(3,4)-Y(0,2)*Y(0,3)*Y(1,4)*Y(3,4)+Y(0,
     &     1)*Y(0,3)*Y(2,2)*Y(3,4)+Y(0,1)*Y(0,4)*Y(2,2)*Y(3,4)-2D0*Y(0,
     &     3)*Y(0,4)*Y(2,2)*Y(3,4)-Y(0,0)*Y(1,3)*Y(2,2)*Y(3,4)+Y(0,4)*
     &     Y(1,3)*Y(2,2)*Y(3,4)-Y(0,0)*Y(1,4)*Y(2,2)*Y(3,4)+Y(0,3)*Y(1,
     &     4)*Y(2,2)*Y(3,4)-Y(0,1)*Y(0,2)*Y(2,3)*Y(3,4)-Y(0,1)*Y(0,4)*
     &     Y(2,3)*Y(3,4)+2D0*Y(0,2)*Y(0,4)*Y(2,3)*Y(3,4)+Y(0,0)*Y(1,2)*
     &     Y(2,3)*Y(3,4)-Y(0,4)*Y(1,2)*Y(2,3)*Y(3,4)+Y(0,0)*Y(1,4)*Y(2,
     &     3)*Y(3,4)-Y(0,2)*Y(1,4)*Y(2,3)*Y(3,4)-Y(0,1)*Y(0,2)*Y(2,4)*
     &     Y(3,4)-Y(0,1)*Y(0,3)*Y(2,4)*Y(3,4)+2D0*Y(0,2)*Y(0,3)*Y(2,4)*
     &     Y(3,4)+Y(0,0)*Y(1,2)*Y(2,4)*Y(3,4)-Y(0,3)*Y(1,2)*Y(2,4)*Y(3,
     &     4)+Y(0,0)*Y(1,3)*Y(2,4)*Y(3,4)-Y(0,2)*Y(1,3)*Y(2,4)*Y(3,4)-
     &     2D0*Y(0,0)*Y(2,3)*Y(2,4)*Y(3,4)+2D0*Y(0,1)*Y(2,3)*Y(2,4)*
     &     Y(3,4)+Y(0,1)*Y(0,2)*Y(3,4)**2-Y(0,2)**2*Y(3,4)**2-Y(0,0)*
     &     Y(1,2)*Y(3,4)**2+Y(0,2)*Y(1,2)*Y(3,4)**2+Y(0,0)*Y(2,2)*Y(3,
     &     4)**2-Y(0,1)*Y(2,2)*Y(3,4)**2+Y(0,2)*Y(0,3)*Y(1,2)*Y(4,4)-
     &     Y(0,3)**2*Y(1,2)*Y(4,4)-Y(0,2)**2*Y(1,3)*Y(4,4)+Y(0,2)*Y(0,
     &     3)*Y(1,3)*Y(4,4)-Y(0,1)*Y(0,3)*Y(2,2)*Y(4,4)+Y(0,3)**2*Y(2,
     &     2)*Y(4,4)+Y(0,0)*Y(1,3)*Y(2,2)*Y(4,4)-Y(0,3)*Y(1,3)*Y(2,2)*
     &     Y(4,4)+Y(0,1)*Y(0,2)*Y(2,3)*Y(4,4)+Y(0,1)*Y(0,3)*Y(2,3)*Y(4,
     &     4)-2D0*Y(0,2)*Y(0,3)*Y(2,3)*Y(4,4)-Y(0,0)*Y(1,2)*Y(2,3)*Y(4,
     &     4)+Y(0,3)*Y(1,2)*Y(2,3)*Y(4,4)-Y(0,0)*Y(1,3)*Y(2,3)*Y(4,4)+
     &     Y(0,2)*Y(1,3)*Y(2,3)*Y(4,4)+Y(0,0)*Y(2,3)**2*Y(4,4)-Y(0,1)*
     &     Y(2,3)**2*Y(4,4)-Y(0,1)*Y(0,2)*Y(3,3)*Y(4,4)+Y(0,2)**2*Y(3,
     &     3)*Y(4,4)+Y(0,0)*Y(1,2)*Y(3,3)*Y(4,4)-Y(0,2)*Y(1,2)*Y(3,3)*
     &     Y(4,4)-Y(0,0)*Y(2,2)*Y(3,3)*Y(4,4)+Y(0,1)*Y(2,2)*Y(3,3)*Y(4,
     &     4)
      DETY(2) = -(Y(0,3)*Y(0,4)*Y(1,2)*Y(1,3))+Y(0,4)**2*Y(1,2)*
     &     Y(1,3)+Y(0,2)*Y(0,4)*Y(1,3)**2-Y(0,4)**2*Y(1,3)**2+Y(0,3)**
     &     2*Y(1,2)*Y(1,4)-Y(0,3)*Y(0,4)*Y(1,2)*Y(1,4)-Y(0,2)*Y(0,3)*
     &     Y(1,3)*Y(1,4)-Y(0,2)*Y(0,4)*Y(1,3)*Y(1,4)+2D0*Y(0,3)*Y(0,4)*
     &     Y(1,3)*Y(1,4)+Y(0,2)*Y(0,3)*Y(1,4)**2-Y(0,3)**2*Y(1,4)**2+
     &     Y(0,3)*Y(0,4)*Y(1,1)*Y(2,3)-Y(0,4)**2*Y(1,1)*Y(2,3)-Y(0,1)*
     &     Y(0,4)*Y(1,3)*Y(2,3)+Y(0,4)**2*Y(1,3)*Y(2,3)-Y(0,1)*Y(0,3)*
     &     Y(1,4)*Y(2,3)+2D0*Y(0,1)*Y(0,4)*Y(1,4)*Y(2,3)-Y(0,3)*Y(0,4)*
     &     Y(1,4)*Y(2,3)+Y(0,0)*Y(1,3)*Y(1,4)*Y(2,3)-Y(0,4)*Y(1,3)*Y(1,
     &     4)*Y(2,3)-Y(0,0)*Y(1,4)**2*Y(2,3)+Y(0,3)*Y(1,4)**2*Y(2,3)-
     &     Y(0,3)**2*Y(1,1)*Y(2,4)+Y(0,3)*Y(0,4)*Y(1,1)*Y(2,4)+2D0*Y(0,
     &     1)*Y(0,3)*Y(1,3)*Y(2,4)-Y(0,1)*Y(0,4)*Y(1,3)*Y(2,4)-Y(0,3)*
     &     Y(0,4)*Y(1,3)*Y(2,4)-Y(0,0)*Y(1,3)**2*Y(2,4)+Y(0,4)*Y(1,3)**
     &     2*Y(2,4)-Y(0,1)*Y(0,3)*Y(1,4)*Y(2,4)+Y(0,3)**2*Y(1,4)*Y(2,
     &     4)+Y(0,0)*Y(1,3)*Y(1,4)*Y(2,4)-Y(0,3)*Y(1,3)*Y(1,4)*Y(2,4)-
     &     Y(0,2)*Y(0,4)*Y(1,1)*Y(3,3)+Y(0,4)**2*Y(1,1)*Y(3,3)+Y(0,1)*
     &     Y(0,4)*Y(1,2)*Y(3,3)-Y(0,4)**2*Y(1,2)*Y(3,3)+Y(0,1)*Y(0,2)*
     &     Y(1,4)*Y(3,3)-2D0*Y(0,1)*Y(0,4)*Y(1,4)*Y(3,3)+Y(0,2)*Y(0,4)*
     &     Y(1,4)*Y(3,3)-Y(0,0)*Y(1,2)*Y(1,4)*Y(3,3)+Y(0,4)*Y(1,2)*Y(1,
     &     4)*Y(3,3)+Y(0,0)*Y(1,4)**2*Y(3,3)-Y(0,2)*Y(1,4)**2*Y(3,3)-
     &     Y(0,1)**2*Y(2,4)*Y(3,3)+Y(0,1)*Y(0,4)*Y(2,4)*Y(3,3)+Y(0,0)*
     &     Y(1,1)*Y(2,4)*Y(3,3)-Y(0,4)*Y(1,1)*Y(2,4)*Y(3,3)-Y(0,0)*Y(1,
     &     4)*Y(2,4)*Y(3,3)+Y(0,1)*Y(1,4)*Y(2,4)*Y(3,3)+Y(0,2)*Y(0,3)*
     &     Y(1,1)*Y(3,4)+Y(0,2)*Y(0,4)*Y(1,1)*Y(3,4)-2D0*Y(0,3)*Y(0,4)*
     &     Y(1,1)*Y(3,4)-Y(0,1)*Y(0,3)*Y(1,2)*Y(3,4)-Y(0,1)*Y(0,4)*Y(1,
     &     2)*Y(3,4)+2D0*Y(0,3)*Y(0,4)*Y(1,2)*Y(3,4)-Y(0,1)*Y(0,2)*Y(1,
     &     3)*Y(3,4)+2D0*Y(0,1)*Y(0,4)*Y(1,3)*Y(3,4)-Y(0,2)*Y(0,4)*Y(1,
     &     3)*Y(3,4)+Y(0,0)*Y(1,2)*Y(1,3)*Y(3,4)-Y(0,4)*Y(1,2)*Y(1,3)*
     &     Y(3,4)-Y(0,1)*Y(0,2)*Y(1,4)*Y(3,4)+2D0*Y(0,1)*Y(0,3)*Y(1,4)*
     &     Y(3,4)-Y(0,2)*Y(0,3)*Y(1,4)*Y(3,4)+Y(0,0)*Y(1,2)*Y(1,4)*Y(3,
     &     4)-Y(0,3)*Y(1,2)*Y(1,4)*Y(3,4)-2D0*Y(0,0)*Y(1,3)*Y(1,4)*Y(3,
     &     4)+2D0*Y(0,2)*Y(1,3)*Y(1,4)*Y(3,4)+Y(0,1)**2*Y(2,3)*Y(3,4)-
     &     Y(0,1)*Y(0,4)*Y(2,3)*Y(3,4)-Y(0,0)*Y(1,1)*Y(2,3)*Y(3,4)+Y(0,
     &     4)*Y(1,1)*Y(2,3)*Y(3,4)+Y(0,0)*Y(1,4)*Y(2,3)*Y(3,4)-Y(0,1)*
     &     Y(1,4)*Y(2,3)*Y(3,4)+Y(0,1)**2*Y(2,4)*Y(3,4)-Y(0,1)*Y(0,3)*
     &     Y(2,4)*Y(3,4)-Y(0,0)*Y(1,1)*Y(2,4)*Y(3,4)+Y(0,3)*Y(1,1)*Y(2,
     &     4)*Y(3,4)+Y(0,0)*Y(1,3)*Y(2,4)*Y(3,4)-Y(0,1)*Y(1,3)*Y(2,4)*
     &     Y(3,4)-Y(0,1)**2*Y(3,4)**2+Y(0,1)*Y(0,2)*Y(3,4)**2+Y(0,0)*
     &     Y(1,1)*Y(3,4)**2-Y(0,2)*Y(1,1)*Y(3,4)**2-Y(0,0)*Y(1,2)*Y(3,
     &     4)**2+Y(0,1)*Y(1,2)*Y(3,4)**2-Y(0,2)*Y(0,3)*Y(1,1)*Y(4,4)+
     &     Y(0,3)**2*Y(1,1)*Y(4,4)+Y(0,1)*Y(0,3)*Y(1,2)*Y(4,4)-Y(0,3)**
     &     2*Y(1,2)*Y(4,4)+Y(0,1)*Y(0,2)*Y(1,3)*Y(4,4)-2D0*Y(0,1)*Y(0,
     &     3)*Y(1,3)*Y(4,4)+Y(0,2)*Y(0,3)*Y(1,3)*Y(4,4)-Y(0,0)*Y(1,2)*
     &     Y(1,3)*Y(4,4)+Y(0,3)*Y(1,2)*Y(1,3)*Y(4,4)+Y(0,0)*Y(1,3)**2*
     &     Y(4,4)-Y(0,2)*Y(1,3)**2*Y(4,4)-Y(0,1)**2*Y(2,3)*Y(4,4)+Y(0,
     &     1)*Y(0,3)*Y(2,3)*Y(4,4)+Y(0,0)*Y(1,1)*Y(2,3)*Y(4,4)-Y(0,3)*
     &     Y(1,1)*Y(2,3)*Y(4,4)-Y(0,0)*Y(1,3)*Y(2,3)*Y(4,4)+Y(0,1)*Y(1,
     &     3)*Y(2,3)*Y(4,4)+Y(0,1)**2*Y(3,3)*Y(4,4)-Y(0,1)*Y(0,2)*Y(3,
     &     3)*Y(4,4)-Y(0,0)*Y(1,1)*Y(3,3)*Y(4,4)+Y(0,2)*Y(1,1)*Y(3,3)*
     &     Y(4,4)+Y(0,0)*Y(1,2)*Y(3,3)*Y(4,4)-Y(0,1)*Y(1,2)*Y(3,3)*Y(4,
     &     4)
      DETY(3) = Y(0,3)*Y(0,4)*Y(1,2)**2-Y(0,4)**2*Y(1,2)**2-Y(0,
     &     2)*Y(0,4)*Y(1,2)*Y(1,3)+Y(0,4)**2*Y(1,2)*Y(1,3)-Y(0,2)*Y(0,
     &     3)*Y(1,2)*Y(1,4)+2D0*Y(0,2)*Y(0,4)*Y(1,2)*Y(1,4)-Y(0,3)*Y(0,
     &     4)*Y(1,2)*Y(1,4)+Y(0,2)**2*Y(1,3)*Y(1,4)-Y(0,2)*Y(0,4)*Y(1,
     &     3)*Y(1,4)-Y(0,2)**2*Y(1,4)**2+Y(0,2)*Y(0,3)*Y(1,4)**2-Y(0,
     &     3)*Y(0,4)*Y(1,1)*Y(2,2)+Y(0,4)**2*Y(1,1)*Y(2,2)+Y(0,1)*Y(0,
     &     4)*Y(1,3)*Y(2,2)-Y(0,4)**2*Y(1,3)*Y(2,2)+Y(0,1)*Y(0,3)*Y(1,
     &     4)*Y(2,2)-2D0*Y(0,1)*Y(0,4)*Y(1,4)*Y(2,2)+Y(0,3)*Y(0,4)*Y(1,
     &     4)*Y(2,2)-Y(0,0)*Y(1,3)*Y(1,4)*Y(2,2)+Y(0,4)*Y(1,3)*Y(1,4)*
     &     Y(2,2)+Y(0,0)*Y(1,4)**2*Y(2,2)-Y(0,3)*Y(1,4)**2*Y(2,2)+Y(0,
     &     2)*Y(0,4)*Y(1,1)*Y(2,3)-Y(0,4)**2*Y(1,1)*Y(2,3)-Y(0,1)*Y(0,
     &     4)*Y(1,2)*Y(2,3)+Y(0,4)**2*Y(1,2)*Y(2,3)-Y(0,1)*Y(0,2)*Y(1,
     &     4)*Y(2,3)+2D0*Y(0,1)*Y(0,4)*Y(1,4)*Y(2,3)-Y(0,2)*Y(0,4)*Y(1,
     &     4)*Y(2,3)+Y(0,0)*Y(1,2)*Y(1,4)*Y(2,3)-Y(0,4)*Y(1,2)*Y(1,4)*
     &     Y(2,3)-Y(0,0)*Y(1,4)**2*Y(2,3)+Y(0,2)*Y(1,4)**2*Y(2,3)+Y(0,
     &     2)*Y(0,3)*Y(1,1)*Y(2,4)-2D0*Y(0,2)*Y(0,4)*Y(1,1)*Y(2,4)+Y(0,
     &     3)*Y(0,4)*Y(1,1)*Y(2,4)-Y(0,1)*Y(0,3)*Y(1,2)*Y(2,4)+2D0*Y(0,
     &     1)*Y(0,4)*Y(1,2)*Y(2,4)-Y(0,3)*Y(0,4)*Y(1,2)*Y(2,4)-Y(0,1)*
     &     Y(0,2)*Y(1,3)*Y(2,4)-Y(0,1)*Y(0,4)*Y(1,3)*Y(2,4)+2D0*Y(0,2)*
     &     Y(0,4)*Y(1,3)*Y(2,4)+Y(0,0)*Y(1,2)*Y(1,3)*Y(2,4)-Y(0,4)*Y(1,
     &     2)*Y(1,3)*Y(2,4)+2D0*Y(0,1)*Y(0,2)*Y(1,4)*Y(2,4)-Y(0,1)*Y(0,
     &     3)*Y(1,4)*Y(2,4)-Y(0,2)*Y(0,3)*Y(1,4)*Y(2,4)-2D0*Y(0,0)*Y(1,
     &     2)*Y(1,4)*Y(2,4)+2D0*Y(0,3)*Y(1,2)*Y(1,4)*Y(2,4)+Y(0,0)*Y(1,
     &     3)*Y(1,4)*Y(2,4)-Y(0,2)*Y(1,3)*Y(1,4)*Y(2,4)+Y(0,1)**2*Y(2,
     &     3)*Y(2,4)-Y(0,1)*Y(0,4)*Y(2,3)*Y(2,4)-Y(0,0)*Y(1,1)*Y(2,3)*
     &     Y(2,4)+Y(0,4)*Y(1,1)*Y(2,3)*Y(2,4)+Y(0,0)*Y(1,4)*Y(2,3)*Y(2,
     &     4)-Y(0,1)*Y(1,4)*Y(2,3)*Y(2,4)-Y(0,1)**2*Y(2,4)**2+Y(0,1)*
     &     Y(0,3)*Y(2,4)**2+Y(0,0)*Y(1,1)*Y(2,4)**2-Y(0,3)*Y(1,1)*Y(2,
     &     4)**2-Y(0,0)*Y(1,3)*Y(2,4)**2+Y(0,1)*Y(1,3)*Y(2,4)**2-Y(0,
     &     2)**2*Y(1,1)*Y(3,4)+Y(0,2)*Y(0,4)*Y(1,1)*Y(3,4)+2D0*Y(0,1)*
     &     Y(0,2)*Y(1,2)*Y(3,4)-Y(0,1)*Y(0,4)*Y(1,2)*Y(3,4)-Y(0,2)*Y(0,
     &     4)*Y(1,2)*Y(3,4)-Y(0,0)*Y(1,2)**2*Y(3,4)+Y(0,4)*Y(1,2)**2*
     &     Y(3,4)-Y(0,1)*Y(0,2)*Y(1,4)*Y(3,4)+Y(0,2)**2*Y(1,4)*Y(3,4)+
     &     Y(0,0)*Y(1,2)*Y(1,4)*Y(3,4)-Y(0,2)*Y(1,2)*Y(1,4)*Y(3,4)-Y(0,
     &     1)**2*Y(2,2)*Y(3,4)+Y(0,1)*Y(0,4)*Y(2,2)*Y(3,4)+Y(0,0)*Y(1,
     &     1)*Y(2,2)*Y(3,4)-Y(0,4)*Y(1,1)*Y(2,2)*Y(3,4)-Y(0,0)*Y(1,4)*
     &     Y(2,2)*Y(3,4)+Y(0,1)*Y(1,4)*Y(2,2)*Y(3,4)+Y(0,1)**2*Y(2,4)*
     &     Y(3,4)-Y(0,1)*Y(0,2)*Y(2,4)*Y(3,4)-Y(0,0)*Y(1,1)*Y(2,4)*Y(3,
     &     4)+Y(0,2)*Y(1,1)*Y(2,4)*Y(3,4)+Y(0,0)*Y(1,2)*Y(2,4)*Y(3,4)-
     &     Y(0,1)*Y(1,2)*Y(2,4)*Y(3,4)+Y(0,2)**2*Y(1,1)*Y(4,4)-Y(0,2)*
     &     Y(0,3)*Y(1,1)*Y(4,4)-2D0*Y(0,1)*Y(0,2)*Y(1,2)*Y(4,4)+Y(0,1)*
     &     Y(0,3)*Y(1,2)*Y(4,4)+Y(0,2)*Y(0,3)*Y(1,2)*Y(4,4)+Y(0,0)*Y(1,
     &     2)**2*Y(4,4)-Y(0,3)*Y(1,2)**2*Y(4,4)+Y(0,1)*Y(0,2)*Y(1,3)*
     &     Y(4,4)-Y(0,2)**2*Y(1,3)*Y(4,4)-Y(0,0)*Y(1,2)*Y(1,3)*Y(4,4)+
     &     Y(0,2)*Y(1,2)*Y(1,3)*Y(4,4)+Y(0,1)**2*Y(2,2)*Y(4,4)-Y(0,1)*
     &     Y(0,3)*Y(2,2)*Y(4,4)-Y(0,0)*Y(1,1)*Y(2,2)*Y(4,4)+Y(0,3)*Y(1,
     &     1)*Y(2,2)*Y(4,4)+Y(0,0)*Y(1,3)*Y(2,2)*Y(4,4)-Y(0,1)*Y(1,3)*
     &     Y(2,2)*Y(4,4)-Y(0,1)**2*Y(2,3)*Y(4,4)+Y(0,1)*Y(0,2)*Y(2,3)*
     &     Y(4,4)+Y(0,0)*Y(1,1)*Y(2,3)*Y(4,4)-Y(0,2)*Y(1,1)*Y(2,3)*Y(4,
     &     4)-Y(0,0)*Y(1,2)*Y(2,3)*Y(4,4)+Y(0,1)*Y(1,2)*Y(2,3)*Y(4,4)
      DETY(4) = -(Y(0,3)**2*Y(1,2)**2)+Y(0,3)*Y(0,4)*Y(1,2)**2+
     &     2D0*Y(0,2)*Y(0,3)*Y(1,2)*Y(1,3)-Y(0,2)*Y(0,4)*Y(1,2)*Y(1,3)-
     &     Y(0,3)*Y(0,4)*Y(1,2)*Y(1,3)-Y(0,2)**2*Y(1,3)**2+Y(0,2)*Y(0,
     &     4)*Y(1,3)**2-Y(0,2)*Y(0,3)*Y(1,2)*Y(1,4)+Y(0,3)**2*Y(1,2)*
     &     Y(1,4)+Y(0,2)**2*Y(1,3)*Y(1,4)-Y(0,2)*Y(0,3)*Y(1,3)*Y(1,4)+
     &     Y(0,3)**2*Y(1,1)*Y(2,2)-Y(0,3)*Y(0,4)*Y(1,1)*Y(2,2)-2D0*Y(0,
     &     1)*Y(0,3)*Y(1,3)*Y(2,2)+Y(0,1)*Y(0,4)*Y(1,3)*Y(2,2)+Y(0,3)*
     &     Y(0,4)*Y(1,3)*Y(2,2)+Y(0,0)*Y(1,3)**2*Y(2,2)-Y(0,4)*Y(1,3)**
     &     2*Y(2,2)+Y(0,1)*Y(0,3)*Y(1,4)*Y(2,2)-Y(0,3)**2*Y(1,4)*Y(2,
     &     2)-Y(0,0)*Y(1,3)*Y(1,4)*Y(2,2)+Y(0,3)*Y(1,3)*Y(1,4)*Y(2,2)-
     &     2D0*Y(0,2)*Y(0,3)*Y(1,1)*Y(2,3)+Y(0,2)*Y(0,4)*Y(1,1)*Y(2,3)+
     &     Y(0,3)*Y(0,4)*Y(1,1)*Y(2,3)+2D0*Y(0,1)*Y(0,3)*Y(1,2)*Y(2,3)-
     &     Y(0,1)*Y(0,4)*Y(1,2)*Y(2,3)-Y(0,3)*Y(0,4)*Y(1,2)*Y(2,3)+2D0*
     &     Y(0,1)*Y(0,2)*Y(1,3)*Y(2,3)-Y(0,1)*Y(0,4)*Y(1,3)*Y(2,3)-Y(0,
     &     2)*Y(0,4)*Y(1,3)*Y(2,3)-2D0*Y(0,0)*Y(1,2)*Y(1,3)*Y(2,3)+2D0*
     &     Y(0,4)*Y(1,2)*Y(1,3)*Y(2,3)-Y(0,1)*Y(0,2)*Y(1,4)*Y(2,3)-Y(0,
     &     1)*Y(0,3)*Y(1,4)*Y(2,3)+2D0*Y(0,2)*Y(0,3)*Y(1,4)*Y(2,3)+Y(0,
     &     0)*Y(1,2)*Y(1,4)*Y(2,3)-Y(0,3)*Y(1,2)*Y(1,4)*Y(2,3)+Y(0,0)*
     &     Y(1,3)*Y(1,4)*Y(2,3)-Y(0,2)*Y(1,3)*Y(1,4)*Y(2,3)-Y(0,1)**2*
     &     Y(2,3)**2+Y(0,1)*Y(0,4)*Y(2,3)**2+Y(0,0)*Y(1,1)*Y(2,3)**2-
     &     Y(0,4)*Y(1,1)*Y(2,3)**2-Y(0,0)*Y(1,4)*Y(2,3)**2+Y(0,1)*Y(1,
     &     4)*Y(2,3)**2+Y(0,2)*Y(0,3)*Y(1,1)*Y(2,4)-Y(0,3)**2*Y(1,1)*
     &     Y(2,4)-Y(0,1)*Y(0,3)*Y(1,2)*Y(2,4)+Y(0,3)**2*Y(1,2)*Y(2,4)-
     &     Y(0,1)*Y(0,2)*Y(1,3)*Y(2,4)+2D0*Y(0,1)*Y(0,3)*Y(1,3)*Y(2,4)-
     &     Y(0,2)*Y(0,3)*Y(1,3)*Y(2,4)+Y(0,0)*Y(1,2)*Y(1,3)*Y(2,4)-Y(0,
     &     3)*Y(1,2)*Y(1,3)*Y(2,4)-Y(0,0)*Y(1,3)**2*Y(2,4)+Y(0,2)*Y(1,
     &     3)**2*Y(2,4)+Y(0,1)**2*Y(2,3)*Y(2,4)-Y(0,1)*Y(0,3)*Y(2,3)*
     &     Y(2,4)-Y(0,0)*Y(1,1)*Y(2,3)*Y(2,4)+Y(0,3)*Y(1,1)*Y(2,3)*Y(2,
     &     4)+Y(0,0)*Y(1,3)*Y(2,3)*Y(2,4)-Y(0,1)*Y(1,3)*Y(2,3)*Y(2,4)+
     &     Y(0,2)**2*Y(1,1)*Y(3,3)-Y(0,2)*Y(0,4)*Y(1,1)*Y(3,3)-2D0*Y(0,
     &     1)*Y(0,2)*Y(1,2)*Y(3,3)+Y(0,1)*Y(0,4)*Y(1,2)*Y(3,3)+Y(0,2)*
     &     Y(0,4)*Y(1,2)*Y(3,3)+Y(0,0)*Y(1,2)**2*Y(3,3)-Y(0,4)*Y(1,2)**
     &     2*Y(3,3)+Y(0,1)*Y(0,2)*Y(1,4)*Y(3,3)-Y(0,2)**2*Y(1,4)*Y(3,
     &     3)-Y(0,0)*Y(1,2)*Y(1,4)*Y(3,3)+Y(0,2)*Y(1,2)*Y(1,4)*Y(3,3)+
     &     Y(0,1)**2*Y(2,2)*Y(3,3)-Y(0,1)*Y(0,4)*Y(2,2)*Y(3,3)-Y(0,0)*
     &     Y(1,1)*Y(2,2)*Y(3,3)+Y(0,4)*Y(1,1)*Y(2,2)*Y(3,3)+Y(0,0)*Y(1,
     &     4)*Y(2,2)*Y(3,3)-Y(0,1)*Y(1,4)*Y(2,2)*Y(3,3)-Y(0,1)**2*Y(2,
     &     4)*Y(3,3)+Y(0,1)*Y(0,2)*Y(2,4)*Y(3,3)+Y(0,0)*Y(1,1)*Y(2,4)*
     &     Y(3,3)-Y(0,2)*Y(1,1)*Y(2,4)*Y(3,3)-Y(0,0)*Y(1,2)*Y(2,4)*Y(3,
     &     3)+Y(0,1)*Y(1,2)*Y(2,4)*Y(3,3)-Y(0,2)**2*Y(1,1)*Y(3,4)+Y(0,
     &     2)*Y(0,3)*Y(1,1)*Y(3,4)+2D0*Y(0,1)*Y(0,2)*Y(1,2)*Y(3,4)-Y(0,
     &     1)*Y(0,3)*Y(1,2)*Y(3,4)-Y(0,2)*Y(0,3)*Y(1,2)*Y(3,4)-Y(0,0)*
     &     Y(1,2)**2*Y(3,4)+Y(0,3)*Y(1,2)**2*Y(3,4)-Y(0,1)*Y(0,2)*Y(1,
     &     3)*Y(3,4)+Y(0,2)**2*Y(1,3)*Y(3,4)+Y(0,0)*Y(1,2)*Y(1,3)*Y(3,
     &     4)-Y(0,2)*Y(1,2)*Y(1,3)*Y(3,4)-Y(0,1)**2*Y(2,2)*Y(3,4)+Y(0,
     &     1)*Y(0,3)*Y(2,2)*Y(3,4)+Y(0,0)*Y(1,1)*Y(2,2)*Y(3,4)-Y(0,3)*
     &     Y(1,1)*Y(2,2)*Y(3,4)-Y(0,0)*Y(1,3)*Y(2,2)*Y(3,4)+Y(0,1)*Y(1,
     &     3)*Y(2,2)*Y(3,4)+Y(0,1)**2*Y(2,3)*Y(3,4)-Y(0,1)*Y(0,2)*Y(2,
     &     3)*Y(3,4)-Y(0,0)*Y(1,1)*Y(2,3)*Y(3,4)+Y(0,2)*Y(1,1)*Y(2,3)*
     &     Y(3,4)+Y(0,0)*Y(1,2)*Y(2,3)*Y(3,4)-Y(0,1)*Y(1,2)*Y(2,3)*Y(3,
     &     4)

      E0_ir=dcmplx(0d0,0d0)
      do i=0,4
         E0_ir=E0_ir+d0(i)*dety(i)/detyy
      enddo
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine dmn(q1,q2,q3,q4,q5,q6,m1,m2,m3,m4,d1,d2,d3,ch1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       vierpunktfkt. mit bis zu 4 integrationsimpulsen in zaehler
c       realteil
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       05.05.87 sa
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(a-y)
      implicit complex*16(z)
      real*8 c1123(0:2),c2123(0:3),c1134(0:2),c2134(0:3),d2(0:6)
      real*8 c1234(0:2),c2234(0:3),c1124(0:2),c2124(0:3),d1(0:3)
      real*8 d3(0:3,0:3),ch1(0:2)
      complex*16 d0gen,d0_
*
      q7 = q1+q2+q3+q6-q4-q5
      call chutmn(q2,q3,q7,m2,m3,m4,c1234,c2234)
      call chutmn(q1,q2,q4,m1,m2,m3,c1123,c2123)
      call chutmn(q1,q7,q6,m1,m2,m4,c1124,c2124)
      call chutmn(q4,q3,q6,m1,m3,m4,c1134,c2134)
c      write(6,*)q2,q3,q7,m2,m3,m4,c1234(1)
      ch1(0) = c1234(0)
      ch1(1) = -c1234(0)-c1234(1)-c1234(2)
      ch1(2) = c1234(2)
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      p12 = (q1+q4-q2)/2d0
      p13 = (q1+q6-q7)/2d0
      p23 = (q4+q6-q3)/2d0
      det = q1*q4*q6+2d0*p12*p13*p23-q1*p23*p23-q4*p13*p13-q6*p12*p12
      mat11 = q4*q6-p23*p23
      mat12 = p13*p23-q6*p12
      mat13 = p12*p23-q4*p13
      mat21 = mat12
      mat22 = q1*q6-p13*p13
      mat23 = p12*p13-q1*p23
      mat31 = mat13
      mat32 = mat23
      mat33 = q1*q4-p12*p12
      cf1 = q1+m12-m22
      cf2 = q4+m12-m32
      cf3 = q6+m12-m42
*
c old version:
c      if ((m22.lt.q1).and.(m32.lt.q1).and.(m42.lt.q1)) goto 50
c      d1(0) = ggttf(q7,q4,q1,m12,m22)
c      write(6,*)m1,m2,m3,m4,q1,q2,q3,q4,q5,q6
c      if (m12.ge.4d0*m22) then
c         d1(0) = dreal(d0gen(q1,q2,q3,q4,q5,q6,m1,m2,m3,m4))
c      end if
c      goto 100
c--->   mass singular 3 fermion box
c 50   continue
c new version:
c      if (dabs(m1-dsqrt(q1)).lt.1d-2) then
c         m1 = m1+1d-2
c         m12 = m1**2
c      end if
c old version:
c      d1(0) = ggttb(q7,q4,q1,m12,m22)
c      goto 100
c 100  s11 = (c1134(0)-c1234(0)-cf1*d1(0))/2d0
c new version:
      goto 100
c 100  d0new = dreal(D0_(q1,q2,q3,q6,q4,q7,M1,M2,M3,M4,0))
 100  d0new = dreal(D0_(q1,q3,q7,q4,q6,q2,M1,M2,M3,M4,0))
      d1(0) = d0new
      s11 = (c1134(0)-c1234(0)-cf1*d1(0))/2d0
      s12 = (c1124(0)-c1234(0)-cf2*d1(0))/2d0
      s13 = (c1123(0)-c1234(0)-cf3*d1(0))/2d0
      d1(1) = (mat11*s11+mat12*s12+mat13*s13)/det
      d1(2) = (mat21*s11+mat22*s12+mat23*s13)/det
      d1(3) = (mat31*s11+mat32*s12+mat33*s13)/det
      ccc = c1234(1)+c1234(2)+c1234(0)
      s20 = c1234(0)+m12*d1(0)
      s211 = (ccc-cf1*d1(1))/2d0
      s212 = (c1134(1)-c1234(1)-cf1*d1(2))/2d0
      s213 = (c1134(2)-c1234(2)-cf1*d1(3))/2d0
      s221 = (c1124(1)+ccc-cf2*d1(1))/2d0
      s222 = -(c1234(1)+cf2*d1(2))/2d0
      s223 = (c1124(2)-c1234(2)-cf2*d1(3))/2d0
      s231 = (c1123(1)+ccc-cf3*d1(1))/2d0
      s232 = (c1123(2)-c1234(1)-cf3*d1(2))/2d0
      s233 = -(c1234(2)+cf3*d1(3))/2d0
      d2(0) = s20-s211-s222-s233
      d2(1) = (mat11*(s211-d2(0))+mat12*s221+mat13*s231)/det
      d2(2) = (mat21*s212+mat22*(s222-d2(0))+mat23*s232)/det
      d2(3) = (mat31*s213+mat32*s223+mat33*(s233-d2(0)))/det
c---  >   d2(1,2)
      d2(4) = (mat11*s212+mat12*(s222-d2(0))+mat13*s232)/det
c---  >   d2(1,3)
      d2(5) = (mat11*s213+mat12*s223+mat13*(s233-d2(0)))/det
c---  >   d2(2,3)=d2(2,1) with m2<-->m4 for square boxes
      d2(6) = (mat21*s213+mat22*s223+mat23*(s233-d2(0)))/det
      s310 = (c2134(0)-c2234(0)-cf1*d2(0))/2d0
      s320 = (c2124(0)-c2234(0)-cf2*d2(0))/2d0
      s330 = (c2123(0)-c2234(0)-cf3*d2(0))/2d0
      ccc = ccc+(c2234(1)+c2234(2)-c1234(0))/2d0+c2234(3)
      s311 = -ccc-cf1*d2(1)/2d0
      s321 = (c2124(1)-cf2*d2(1))/2d0-ccc
      s331 = (c2123(1)-cf3*d2(1))/2d0-ccc
      s312 = (c2134(1)-c2234(1)-cf1*d2(2))/2d0
      s322 = -(c2234(1)+cf2*d2(2))/2d0
      s332 = (c2123(2)-c2234(1)-cf3*d2(2))/2d0
      s313 = (c2134(2)-c2234(2)-cf1*d2(3))/2d0
      s323 = (c2124(2)-c2234(2)-cf2*d2(3))/2d0
      s333 = -(c2234(2)+cf3*d2(3))/2d0
      s3113 = (c2234(3)+c2234(2)+c1234(2)-cf1*d2(5))/2d0
      s3213 = (c2124(3)+c2234(3)+c2234(2)+c1234(2)-cf2*d2(5))/2d0
      s3313 = (c2234(3)+c2234(2)+c1234(2)-cf3*d2(5))/2d0
      d3(0,1) = (mat11*s310+mat12*s320+mat13*s330)/det
      d3(0,2) = (mat21*s310+mat22*s320+mat23*s330)/det
      d3(0,3) = (mat31*s310+mat32*s320+mat33*s330)/det
      d3(1,1) = (mat11*(s311-2d0*d3(0,1))+mat12*s321+mat13*s331)/det
      d3(3,3) = (mat31*s313+mat32*s323+mat33*(s333-2d0*d3(0,3)))/det
      d3(1,2) = (mat21*(s311-2d0*d3(0,1))+mat22*s321+mat23*s331)/det
      d3(1,3) = (mat31*(s311-2d0*d3(0,1))+mat32*s321+mat33*s331)/det
      d3(2,1) = (mat11*s312+mat12*(s322-2d0*d3(0,2))+mat13*s332)/det
      d3(2,2) = (mat21*s312+mat22*(s322-2d0*d3(0,2))+mat23*s332)/det
      d3(2,3) = (mat31*s312+mat32*(s322-2d0*d3(0,2))+mat33*s332)/det
      d3(3,1) = (mat11*s313+mat12*s323+mat13*(s333-2d0*d3(0,3)))/det
      d3(3,2) = (mat21*s313+mat22*s323+mat23*(s333-2d0*d3(0,3)))/det
      d3(0,0) = (mat21*(s3113-d3(0,3))+mat22*s3213+mat23*
     1     (s3313-d3(0,1)))/det
      return
      end
