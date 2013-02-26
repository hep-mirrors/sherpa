c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine matri(s,t1,t2,u1,u2,mt,mh,mat)
c----------------------------------------------------------------------
c mat(i,j) = (summe ueber polarisat. der gluonen)*trace(i-channel-
c            born-matrixelem.*(p1s-mt)*j-entwicklungskoeff.*(p2+mt))/4
c  i: s,t,u-channel, j: 1-68 t-channel SMEs
c----------------------------------------------------------------------
      implicit none
      integer i,j
      real*8 s,s1,s2,s3,t1,t2,u1,u2,mt,mh
      real*8 cs1,cs2,ct1,ct2,ct3,cu1,cu2,cu3
      real*8 p1p2,p1p3,p1p4,p2p3,p2p4
      real*8 mat(3,68)
c initialization
      do i=1,3
         do j=1,68
            mat(i,j)=0d0
         end do
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
      ct1=1d0/(s1-mt**2)/(t2-mt**2)
      ct2=1d0/(s2-mt**2)/(t1-mt**2)
      ct3=1d0/(t1-mt**2)/(t2-mt**2)
      cu1=1d0/(s1-mt**2)/(u2-mt**2)
      cu2=1d0/(s2-mt**2)/(u1-mt**2)
      cu3=1d0/(u1-mt**2)/(u2-mt**2)
*
      mat(1,1)=cs1*(-16*mt**2*p1p3 + 16*p1p2*p1p3 + 
     $     16*mt**2*p1p4 - 16*p1p2*p1p4 + 16*p1p4*p2p3 - 
     $     16*p1p3*p2p4) + 
     $     cs2*(16*mt**2*p2p3 - 16*p1p2*p2p3 + 16*p1p4*p2p3 - 
     $     16*mt**2*p2p4 + 16*p1p2*p2p4 - 16*p1p3*p2p4)
      mat(2,1)=ct1*(-16*mt**4 + 4*mt**2*s + 16*mt**2*p1p2 - 4*s*p1p2 - 
     $     16*mt**2*p1p3 + 24*p1p2*p1p3 + 8*p1p2*p1p4 + 
     $     (64*mt**2*p1p3*p1p4)/s - (64*p1p2*p1p3*p1p4)/s - 
     $     8*mt**2*p2p3 + 8*p1p4*p2p3 + 
     $     (16*p1p3*p1p4*p2p3)/s - (16*p1p4**2*p2p3)/s - 
     $     8*mt**2*p2p4 - 8*p1p3*p2p4 - 
     $     (16*p1p3**2*p2p4)/s + (16*p1p3*p1p4*p2p4)/s) + 
     $     ct3*(4*mt**2*s + 16*mt**2*p1p2 - 4*s*p1p2 - 16*p1p2**2 - 
     $     8*mt**2*p1p3 + 8*p1p2*p1p4 + 8*p1p2*p2p3 + 
     $     8*p1p4*p2p3 - (32*mt**2*p1p4*p2p3)/s + 
     $     (32*p1p2*p1p4*p2p3)/s - (16*p1p4**2*p2p3)/s - 
     $     (16*p1p4*p2p3**2)/s - 8*mt**2*p2p4 - 
     $     8*p1p3*p2p4 - (32*mt**2*p1p3*p2p4)/s + 
     $     (32*p1p2*p1p3*p2p4)/s + 
     $     (16*p1p3*p1p4*p2p4)/s + (16*p1p3*p2p3*p2p4)/s)+
     $     ct2*(-16*mt**4 + 4*mt**2*s + 16*mt**2*p1p2 - 4*s*p1p2 - 
     $     8*mt**2*p1p3 - 8*mt**2*p1p4 + 8*p1p2*p2p3 + 
     $     8*p1p4*p2p3 - (16*p1p4*p2p3**2)/s - 
     $     16*mt**2*p2p4 + 24*p1p2*p2p4 - 8*p1p3*p2p4 + 
     $     (64*mt**2*p2p3*p2p4)/s - (64*p1p2*p2p3*p2p4)/s + 
     $     (16*p1p3*p2p3*p2p4)/s + 
     $     (16*p1p4*p2p3*p2p4)/s - (16*p1p3*p2p4**2)/s)
      mat(3,1)=cu1*(-16*mt**4 + 4*mt**2*s + 16*mt**2*p1p2 - 4*s*p1p2 + 
     $     8*p1p2*p1p3 - 16*mt**2*p1p4 + 24*p1p2*p1p4 + 
     $     (64*mt**2*p1p3*p1p4)/s - (64*p1p2*p1p3*p1p4)/s - 
     $     8*mt**2*p2p3 - 8*p1p4*p2p3 + 
     $     (16*p1p3*p1p4*p2p3)/s - (16*p1p4**2*p2p3)/s - 
     $     8*mt**2*p2p4 + 8*p1p3*p2p4 - 
     $     (16*p1p3**2*p2p4)/s + (16*p1p3*p1p4*p2p4)/s) + 
     $     cu3*(4*mt**2*s + 16*mt**2*p1p2 - 4*s*p1p2 - 16*p1p2**2 + 
     $     8*p1p2*p1p3 - 8*mt**2*p1p4 - 8*mt**2*p2p3 - 
     $     8*p1p4*p2p3 - (32*mt**2*p1p4*p2p3)/s + 
     $     (32*p1p2*p1p4*p2p3)/s + 
     $     (16*p1p3*p1p4*p2p3)/s + 8*p1p2*p2p4 + 
     $     8*p1p3*p2p4 - (32*mt**2*p1p3*p2p4)/s + 
     $     (32*p1p2*p1p3*p2p4)/s - (16*p1p3**2*p2p4)/s + 
     $     (16*p1p4*p2p3*p2p4)/s - (16*p1p3*p2p4**2)/s) + 
     $     cu2*(-16*mt**4 + 4*mt**2*s + 16*mt**2*p1p2 - 4*s*p1p2 - 
     $     8*mt**2*p1p3 - 8*mt**2*p1p4 - 16*mt**2*p2p3 + 
     $     24*p1p2*p2p3 - 8*p1p4*p2p3 - 
     $     (16*p1p4*p2p3**2)/s + 8*p1p2*p2p4 + 
     $     8*p1p3*p2p4 + (64*mt**2*p2p3*p2p4)/s - 
     $     (64*p1p2*p2p3*p2p4)/s + 
     $     (16*p1p3*p2p3*p2p4)/s + 
     $     (16*p1p4*p2p3*p2p4)/s - (16*p1p3*p2p4**2)/s)

      mat(1,2)=cs1*(-4*mt**2*s**2 + 4*s**2*p1p2 - 8*mt**2*s*p1p3 + 
     $     8*s*p1p2*p1p3 + 8*mt**2*s*p1p4 - 8*s*p1p2*p1p4 - 
     $     8*s*p1p4*p2p3 - 16*p1p3*p1p4*p2p3 + 
     $     16*p1p4**2*p2p3 - 24*s*p1p3*p2p4 + 
     $     16*p1p3**2*p2p4 - 16*p1p3*p1p4*p2p4) + 
     $     cs2*(-4*mt**2*s**2 + 4*s**2*p1p2 + 8*mt**2*s*p2p3 - 
     $     8*s*p1p2*p2p3 - 8*s*p1p4*p2p3 + 
     $     16*p1p4*p2p3**2 - 8*mt**2*s*p2p4 + 
     $     8*s*p1p2*p2p4 - 24*s*p1p3*p2p4 - 
     $     16*p1p3*p2p3*p2p4 - 16*p1p4*p2p3*p2p4 + 
     $     16*p1p3*p2p4**2)
      mat(2,2)=ct1*(-8*mt**4*s + 8*mt**2*s*p1p2 - 8*mt**2*s*p1p3 + 
     $     16*s*p1p2*p1p3 + 32*mt**2*p1p3*p1p4 - 
     $     32*p1p2*p1p3*p1p4 - 16*mt**2*p1p4*p2p3 - 
     $     32*p1p3*p1p4*p2p3 + 
     $     (64*p1p3*p1p4**2*p2p3)/s - 8*mt**2*s*p2p4 + 
     $     16*mt**2*p1p3*p2p4 - 16*s*p1p3*p2p4 + 
     $     32*p1p3*p1p4*p2p4 - (64*p1p3**2*p1p4*p2p4)/s) +
     $     ct3*(8*mt**2*s*p1p2 - 8*s*p1p2**2 - 8*mt**2*s*p1p3 - 
     $     16*mt**2*p1p4*p2p3 + 32*p1p2*p1p4*p2p3 - 
     $     (32*p1p4**2*p2p3**2)/s - 8*mt**2*s*p2p4 - 
     $     16*mt**2*p1p3*p2p4 - 16*s*p1p3*p2p4 + 
     $     32*p1p3*p1p4*p2p4 + 32*p1p3*p2p3*p2p4 + 
     $     (32*p1p3**2*p2p4**2)/s) + 
     $     ct2*(-8*mt**4*s + 8*mt**2*s*p1p2 - 8*mt**2*s*p1p3 - 
     $     16*mt**2*p1p4*p2p3 - 8*mt**2*s*p2p4 + 
     $     16*s*p1p2*p2p4 + 16*mt**2*p1p3*p2p4 - 
     $     16*s*p1p3*p2p4 + 32*mt**2*p2p3*p2p4 - 
     $     32*p1p2*p2p3*p2p4 + 32*p1p3*p2p3*p2p4 - 
     $     32*p1p4*p2p3*p2p4 + 
     $     (64*p1p4*p2p3**2*p2p4)/s - 
     $     (64*p1p3*p2p3*p2p4**2)/s)
      mat(3,2)=cu1*(-8*mt**4*s + 4*mt**2*s**2 + 8*mt**2*s*p1p2 - 
     $     4*s**2*p1p2 + 8*s*p1p2*p1p3 - 8*mt**2*s*p1p4 + 
     $     8*s*p1p2*p1p4 + 32*mt**2*p1p3*p1p4 - 
     $     32*p1p2*p1p3*p1p4 - 16*mt**2*p1p4*p2p3 + 
     $     8*s*p1p4*p2p3 - 16*p1p3*p1p4*p2p3 - 
     $     16*p1p4**2*p2p3 + (64*p1p3*p1p4**2*p2p3)/s - 
     $     8*mt**2*s*p2p4 + 16*mt**2*p1p3*p2p4 + 
     $     8*s*p1p3*p2p4 - 16*p1p3**2*p2p4 + 
     $     48*p1p3*p1p4*p2p4 - (64*p1p3**2*p1p4*p2p4)/s) +
     $     cu3*(4*mt**2*s**2 + 8*mt**2*s*p1p2 - 4*s**2*p1p2 - 
     $     8*s*p1p2**2 + 8*s*p1p2*p1p3 - 
     $     16*mt**2*p1p4*p2p3 + 8*s*p1p4*p2p3 + 
     $     32*p1p2*p1p4*p2p3 - 16*p1p3*p1p4*p2p3 - 
     $     (32*p1p4**2*p2p3**2)/s + 8*s*p1p2*p2p4 - 
     $     16*mt**2*p1p3*p2p4 + 8*s*p1p3*p2p4 - 
     $     16*p1p3**2*p2p4 - 16*p1p4*p2p3*p2p4 - 
     $     16*p1p3*p2p4**2 + (32*p1p3**2*p2p4**2)/s) + 
     $     cu2*(-8*mt**4*s + 4*mt**2*s**2 +8*mt**2*s*p1p2-4*s**2*p1p2- 
     $     8*mt**2*s*p1p3 - 8*mt**2*s*p2p3 + 8*s*p1p2*p2p3 - 
     $     16*mt**2*p1p4*p2p3 + 8*s*p1p4*p2p3 - 
     $     16*p1p4*p2p3**2 + 8*s*p1p2*p2p4 + 
     $     16*mt**2*p1p3*p2p4 + 8*s*p1p3*p2p4 + 
     $     32*mt**2*p2p3*p2p4 - 32*p1p2*p2p3*p2p4 + 
     $     48*p1p3*p2p3*p2p4 - 16*p1p4*p2p3*p2p4 + 
     $     (64*p1p4*p2p3**2*p2p4)/s - 16*p1p3*p2p4**2 - 
     $     (64*p1p3*p2p3*p2p4**2)/s)

      mat(1,3)=cs1*(8*mt*s*p1p4 + 16*mt*p1p3*p1p4 - 
     $     16*mt*p1p4**2 + 8*mt*s*p2p4 - 16*mt*p1p3*p2p4 + 
     $     16*mt*p1p4*p2p4) + 
     $     cs2*(8*mt*s*p1p4 - 16*mt*p1p4*p2p3 + 8*mt*s*p2p4 + 
     $     16*mt*p1p4*p2p4 + 16*mt*p2p3*p2p4 - 16*mt*p2p4**2)
      mat(2,3)=ct1*(-4*mt**3*s - 4*mt*s*p1p2 + 16*mt**3*p1p4 + 
     $     32*mt*p1p3*p1p4 - (64*mt*p1p3*p1p4**2)/s + 
     $     8*mt*p1p4*p2p3 - 16*mt**3*p2p4 + 8*mt*s*p2p4 - 
     $     8*mt*p1p3*p2p4 + (64*mt*p1p3*p1p4*p2p4)/s) + 
     $     ct3*(4*mt**3*s + 4*mt*s*p1p2 - 16*mt*p1p2*p1p4 - 
     $     8*mt*p1p4*p2p3 + (32*mt*p1p4**2*p2p3)/s + 
     $     8*mt*s*p2p4 + 16*mt*p1p2*p2p4 - 8*mt*p1p3*p2p4 + 
     $     (32*mt*p1p3*p1p4*p2p4)/s - 16*mt*p2p3*p2p4 - 
     $     (32*mt*p1p4*p2p3*p2p4)/s - (32*mt*p1p3*p2p4**2)/s) + 
     $     ct2*(4*mt**3*s + 4*mt*s*p1p2 + 16*mt**3*p1p4 - 
     $     8*mt*p1p4*p2p3 - 16*mt**3*p2p4 + 8*mt*s*p2p4 - 
     $     8*mt*p1p3*p2p4 + 16*mt*p1p4*p2p4 - 
     $     16*mt*p2p3*p2p4 - (64*mt*p1p4*p2p3*p2p4)/s - 
     $     16*mt*p2p4**2 + (64*mt*p2p3*p2p4**2)/s)
      mat(3,3)=cu1*(-4*mt**3*s - 4*mt*s*p1p2 + 16*mt**3*p1p4 - 
     $     8*mt*s*p1p4 + 16*mt*p1p3*p1p4 + 16*mt*p1p4**2 - 
     $     (64*mt*p1p3*p1p4**2)/s + 8*mt*p1p4*p2p3 - 
     $     16*mt**3*p2p4 + 8*mt*p1p3*p2p4 - 
     $     16*mt*p1p4*p2p4 + (64*mt*p1p3*p1p4*p2p4)/s) + 
     $     cu3*(-4*mt**3*s - 4*mt*s*p1p2 - 8*mt*s*p1p4 - 
     $     16*mt*p1p2*p1p4 + 16*mt*p1p3*p1p4 + 
     $     8*mt*p1p4*p2p3 + (32*mt*p1p4**2*p2p3)/s + 
     $     16*mt*p1p2*p2p4 + 8*mt*p1p3*p2p4 + 
     $     (32*mt*p1p3*p1p4*p2p4)/s - 
     $     (32*mt*p1p4*p2p3*p2p4)/s - (32*mt*p1p3*p2p4**2)/s) + 
     $     cu2*(4*mt**3*s + 4*mt*s*p1p2 + 16*mt**3*p1p4 - 8*mt*s*p1p4+ 
     $     8*mt*p1p4*p2p3 - 16*mt**3*p2p4 - 8*mt*p1p3*p2p4 - 
     $     32*mt*p2p3*p2p4 - (64*mt*p1p4*p2p3*p2p4)/s + 
     $     (64*mt*p2p3*p2p4**2)/s)

      mat(1,4)=cs1*(-8*mt*s*p1p3 + 16*mt*p1p3**2 - 
     $     16*mt*p1p3*p1p4 - 8*mt*s*p2p3 - 
     $     16*mt*p1p3*p2p3 + 16*mt*p1p4*p2p3) + 
     $     cs2*(-8*mt*s*p1p3 - 8*mt*s*p2p3 - 16*mt*p1p3*p2p3 + 
     $     16*mt*p2p3**2 + 16*mt*p1p3*p2p4 - 16*mt*p2p3*p2p4)
      mat(2,4)=ct1*(-4*mt**3*s - 4*mt*s*p1p2 + 16*mt**3*p1p3 - 
     $     8*mt*s*p1p3 + 16*mt*p1p3**2 + 16*mt*p1p3*p1p4 - 
     $     (64*mt*p1p3**2*p1p4)/s - 16*mt**3*p2p3 - 
     $     16*mt*p1p3*p2p3 + 8*mt*p1p4*p2p3 + 
     $     (64*mt*p1p3*p1p4*p2p3)/s + 8*mt*p1p3*p2p4) + 
     $     ct3*(-4*mt**3*s - 4*mt*s*p1p2 - 8*mt*s*p1p3 - 
     $     16*mt*p1p2*p1p3 + 16*mt*p1p3*p1p4 + 
     $     16*mt*p1p2*p2p3 + 8*mt*p1p4*p2p3 + 
     $     (32*mt*p1p3*p1p4*p2p3)/s - 
     $     (32*mt*p1p4*p2p3**2)/s + 8*mt*p1p3*p2p4 + 
     $     (32*mt*p1p3**2*p2p4)/s - (32*mt*p1p3*p2p3*p2p4)/s) + 
     $     ct2*(4*mt**3*s + 4*mt*s*p1p2 + 16*mt**3*p1p3 - 8*mt*s*p1p3- 
     $     16*mt**3*p2p3 - 8*mt*p1p4*p2p3 + 8*mt*p1p3*p2p4 - 
     $     32*mt*p2p3*p2p4 - (64*mt*p1p3*p2p3*p2p4)/s + 
     $     (64*mt*p2p3**2*p2p4)/s)
      mat(3,4)=cu1*(-4*mt**3*s - 4*mt*s*p1p2 + 16*mt**3*p1p3 + 
     $     32*mt*p1p3*p1p4 - (64*mt*p1p3**2*p1p4)/s - 
     $     16*mt**3*p2p3 + 8*mt*s*p2p3 - 8*mt*p1p4*p2p3 + 
     $     (64*mt*p1p3*p1p4*p2p3)/s + 8*mt*p1p3*p2p4) + 
     $     cu3*(4*mt**3*s + 4*mt*s*p1p2 - 16*mt*p1p2*p1p3 + 
     $     8*mt*s*p2p3 + 16*mt*p1p2*p2p3 - 8*mt*p1p4*p2p3 + 
     $     (32*mt*p1p3*p1p4*p2p3)/s - 
     $     (32*mt*p1p4*p2p3**2)/s - 8*mt*p1p3*p2p4 + 
     $     (32*mt*p1p3**2*p2p4)/s - 16*mt*p2p3*p2p4 - 
     $     (32*mt*p1p3*p2p3*p2p4)/s) + 
     $     cu2*(4*mt**3*s + 4*mt*s*p1p2 + 16*mt**3*p1p3 - 
     $     16*mt**3*p2p3 + 8*mt*s*p2p3 + 16*mt*p1p3*p2p3 - 
     $     8*mt*p1p4*p2p3 - 16*mt*p2p3**2 - 8*mt*p1p3*p2p4 - 
     $     16*mt*p2p3*p2p4 - (64*mt*p1p3*p2p3*p2p4)/s + 
     $     (64*mt*p2p3**2*p2p4)/s)

      mat(1,5)=cs1*(-16*mt**2*p1p3 + 16*p1p2*p1p3 + 
     $     16*mt**2*p1p4 - 16*p1p2*p1p4 + 16*p1p4*p2p3 - 
     $     16*p1p3*p2p4) + 
     $     cs2*(16*mt**2*p2p3 - 16*p1p2*p2p3 + 16*p1p4*p2p3 - 
     $     16*mt**2*p2p4 + 16*p1p2*p2p4 - 16*p1p3*p2p4)
      mat(2,5)=ct1*(-16*mt**4 + 8*mt**2*s + 16*mt**2*p1p2 - 
     $     32*mt**2*p1p3 + 16*p1p2*p1p3 + 
     $     (64*mt**2*p1p3*p1p4)/s - (64*p1p2*p1p3*p1p4)/s + 
     $     (32*p1p3*p1p4*p2p3)/s - 16*mt**2*p2p4 - 
     $     32*p1p3*p2p4 + (32*p1p3**2*p2p4)/s + 
     $     (64*p1p3*p1p4*p2p4)/s) + 
     $     ct2*(-16*mt**4 + 8*mt**2*s + 16*mt**2*p1p2 - 16*mt**2*p1p3- 
     $     32*mt**2*p2p4 + 16*p1p2*p2p4 - 32*p1p3*p2p4 + 
     $     (64*mt**2*p2p3*p2p4)/s - (64*p1p2*p2p3*p2p4)/s + 
     $     (64*p1p3*p2p3*p2p4)/s + 
     $     (32*p1p4*p2p3*p2p4)/s + (32*p1p3*p2p4**2)/s) + 
     $     ct3*(-16*mt**4 + 8*mt**2*s + 16*mt**2*p1p2 - 16*mt**2*p1p3+ 
     $     (64*mt**2*p1p3*p1p4)/s - (32*mt**2*p1p4*p2p3)/s - 
     $     (32*p1p2*p1p4*p2p3)/s + 
     $     (64*p1p4**2*p2p3**2)/s**2 - 16*mt**2*p2p4 - 
     $     32*p1p3*p2p4 - (32*mt**2*p1p3*p2p4)/s - 
     $     (32*p1p2*p1p3*p2p4)/s + 
     $     (64*p1p3*p1p4*p2p4)/s + (64*mt**2*p2p3*p2p4)/s + 
     $     (64*p1p3*p2p3*p2p4)/s - 
     $     (128*p1p3*p1p4*p2p3*p2p4)/s**2 + 
     $     (64*p1p3**2*p2p4**2)/s**2)
      mat(3,5)=cu1*(-16*mt**4 + 16*mt**2*p1p2 - 8*s*p1p2 + 
     $     16*p1p2*p1p3 + 32*p1p2*p1p4 + 
     $     (64*mt**2*p1p3*p1p4)/s - (64*p1p2*p1p3*p1p4)/s + 
     $     16*p1p4*p2p3 - (32*p1p3*p1p4*p2p3)/s - 
     $     (64*p1p4**2*p2p3)/s - 16*mt**2*p2p4 + 
     $     16*p1p3*p2p4 - (32*p1p3**2*p2p4)/s) + 
     $     cu2*(-16*mt**4 + 16*mt**2*p1p2 - 8*s*p1p2 - 16*mt**2*p1p3 + 
     $     32*p1p2*p2p3 + 16*p1p4*p2p3 - 
     $     (64*p1p4*p2p3**2)/s + 16*p1p2*p2p4 + 
     $     16*p1p3*p2p4 + (64*mt**2*p2p3*p2p4)/s - 
     $     (64*p1p2*p2p3*p2p4)/s - 
     $     (32*p1p4*p2p3*p2p4)/s - (32*p1p3*p2p4**2)/s) + 
     $     cu3*(16*mt**4 + 16*mt**2*p1p2 - 8*s*p1p2 - 32*p1p2**2 + 
     $     16*p1p2*p1p3 - (64*mt**2*p1p3*p1p4)/s + 
     $     16*p1p4*p2p3 - (32*mt**2*p1p4*p2p3)/s + 
     $     (96*p1p2*p1p4*p2p3)/s - 
     $     (32*p1p3*p1p4*p2p3)/s - 
     $     (64*p1p4**2*p2p3**2)/s**2 + 16*p1p2*p2p4 + 
     $     16*p1p3*p2p4 - (32*mt**2*p1p3*p2p4)/s + 
     $     (96*p1p2*p1p3*p2p4)/s - (32*p1p3**2*p2p4)/s - 
     $     (64*mt**2*p2p3*p2p4)/s - (32*p1p4*p2p3*p2p4)/s + 
     $     (128*p1p3*p1p4*p2p3*p2p4)/s**2 - 
     $     (32*p1p3*p2p4**2)/s - (64*p1p3**2*p2p4**2)/s**2)

      mat(1,6)=cs1*(-4*mt**2*s**2 + 4*s**2*p1p2 - 8*mt**2*s*p1p3 + 
     $     8*s*p1p2*p1p3 + 8*mt**2*s*p1p4 - 8*s*p1p2*p1p4 - 
     $     8*s*p1p4*p2p3 - 16*p1p3*p1p4*p2p3 + 
     $     16*p1p4**2*p2p3 - 24*s*p1p3*p2p4 + 
     $     16*p1p3**2*p2p4 - 16*p1p3*p1p4*p2p4) + 
     $     cs2*(-4*mt**2*s**2 + 4*s**2*p1p2 + 8*mt**2*s*p2p3 - 
     $     8*s*p1p2*p2p3 - 8*s*p1p4*p2p3 + 
     $     16*p1p4*p2p3**2 - 8*mt**2*s*p2p4 + 
     $     8*s*p1p2*p2p4 - 24*s*p1p3*p2p4 - 
     $     16*p1p3*p2p3*p2p4 - 16*p1p4*p2p3*p2p4 + 
     $     16*p1p3*p2p4**2)
      mat(2,6)=ct1*(-8*mt**4*s + 8*mt**2*s*p1p2 - 16*mt**2*s*p1p3 + 
     $     16*s*p1p2*p1p3 + 32*mt**2*p1p3*p1p4 - 
     $     32*p1p2*p1p3*p1p4 - 16*mt**2*p1p4*p2p3 - 
     $     32*p1p3*p1p4*p2p3 + 
     $     (64*p1p3*p1p4**2*p2p3)/s - 16*mt**2*s*p2p4 + 
     $     16*mt**2*p1p3*p2p4 - 32*s*p1p3*p2p4 + 
     $     32*p1p3**2*p2p4 + 64*p1p3*p1p4*p2p4 - 
     $     (64*p1p3**2*p1p4*p2p4)/s) + 
     $     ct3*(-8*mt**4*s + 8*mt**2*s*p1p2 - 16*mt**2*s*p1p3 + 
     $     32*mt**2*p1p3*p1p4 - 16*mt**2*p1p4*p2p3 - 
     $     16*mt**2*s*p2p4 - 16*mt**2*p1p3*p2p4 - 
     $     32*s*p1p3*p2p4 - 32*p1p2*p1p3*p2p4 + 
     $     64*p1p3*p1p4*p2p4 + 32*mt**2*p2p3*p2p4 + 
     $     64*p1p3*p2p3*p2p4 - 
     $     (64*p1p3*p1p4*p2p3*p2p4)/s + 
     $     (64*p1p3**2*p2p4**2)/s) + 
     $     ct2*(-8*mt**4*s + 8*mt**2*s*p1p2 - 16*mt**2*s*p1p3 - 
     $     16*mt**2*p1p4*p2p3 - 16*mt**2*s*p2p4 + 
     $     16*s*p1p2*p2p4 + 16*mt**2*p1p3*p2p4 - 
     $     32*s*p1p3*p2p4 + 32*mt**2*p2p3*p2p4 - 
     $     32*p1p2*p2p3*p2p4 + 64*p1p3*p2p3*p2p4 - 
     $     32*p1p4*p2p3*p2p4 + 
     $     (64*p1p4*p2p3**2*p2p4)/s + 32*p1p3*p2p4**2 - 
     $     (64*p1p3*p2p3*p2p4**2)/s)
      mat(3,6)=cu1*(-8*mt**4*s + 8*mt**2*s*p1p2 - 8*s**2*p1p2 + 
     $     16*s*p1p2*p1p3 + 16*s*p1p2*p1p4 + 
     $     32*mt**2*p1p3*p1p4 - 32*p1p2*p1p3*p1p4 - 
     $     16*mt**2*p1p4*p2p3 + 16*s*p1p4*p2p3 - 
     $     32*p1p3*p1p4*p2p3 - 32*p1p4**2*p2p3 + 
     $     (64*p1p3*p1p4**2*p2p3)/s - 16*mt**2*s*p2p4 + 
     $     16*mt**2*p1p3*p2p4 + 16*s*p1p3*p2p4 - 
     $     32*p1p3**2*p2p4 + 32*p1p3*p1p4*p2p4 - 
     $     (64*p1p3**2*p1p4*p2p4)/s) + 
     $     cu3*(8*mt**4*s + 8*mt**2*s*p1p2 - 8*s**2*p1p2 - 
     $     16*s*p1p2**2 + 16*s*p1p2*p1p3 - 
     $     32*mt**2*p1p3*p1p4 - 16*mt**2*p1p4*p2p3 + 
     $     16*s*p1p4*p2p3 + 64*p1p2*p1p4*p2p3 - 
     $     32*p1p3*p1p4*p2p3 - (64*p1p4**2*p2p3**2)/s + 
     $     16*s*p1p2*p2p4 - 16*mt**2*p1p3*p2p4 + 
     $     16*s*p1p3*p2p4 + 32*p1p2*p1p3*p2p4 - 
     $     32*p1p3**2*p2p4 - 32*mt**2*p2p3*p2p4 - 
     $     32*p1p4*p2p3*p2p4 + 
     $     (64*p1p3*p1p4*p2p3*p2p4)/s - 32*p1p3*p2p4**2) +
     $     cu2*(-8*mt**4*s + 8*mt**2*s*p1p2 - 8*s**2*p1p2 - 
     $     16*mt**2*s*p1p3 + 16*s*p1p2*p2p3 - 
     $     16*mt**2*p1p4*p2p3 + 16*s*p1p4*p2p3 - 
     $     32*p1p4*p2p3**2 + 16*s*p1p2*p2p4 + 
     $     16*mt**2*p1p3*p2p4 + 16*s*p1p3*p2p4 + 
     $     32*mt**2*p2p3*p2p4 - 32*p1p2*p2p3*p2p4 + 
     $     32*p1p3*p2p3*p2p4 - 32*p1p4*p2p3*p2p4 + 
     $     (64*p1p4*p2p3**2*p2p4)/s - 32*p1p3*p2p4**2 - 
     $     (64*p1p3*p2p3*p2p4**2)/s)

      mat(1,7)=cs1*(-8*mt*s*p1p4 - 16*mt*p1p3*p1p4 + 
     $     16*mt*p1p4**2 - 8*mt*s*p2p4 + 16*mt*p1p3*p2p4 - 
     $     16*mt*p1p4*p2p4) + 
     $     cs2*(-8*mt*s*p1p4 + 16*mt*p1p4*p2p3 - 8*mt*s*p2p4 - 
     $     16*mt*p1p4*p2p4 - 16*mt*p2p3*p2p4 + 16*mt*p2p4**2)

      mat(2,7)=ct1*(8*mt*s*p1p2 - 16*mt**3*p1p4 - 
     $     32*mt*p1p3*p1p4 + (64*mt*p1p3*p1p4**2)/s - 
     $     16*mt*p1p4*p2p3 + 16*mt**3*p2p4 - 16*mt*s*p2p4 + 
     $     16*mt*p1p3*p2p4 - (64*mt*p1p3*p1p4*p2p4)/s) + 
     $     ct3*(-8*mt**3*s + 16*mt*p1p2*p1p4 - 
     $     (32*mt*p1p4**2*p2p3)/s - 16*mt*s*p2p4 - 
     $     16*mt*p1p2*p2p4 - (32*mt*p1p3*p1p4*p2p4)/s + 
     $     32*mt*p2p3*p2p4 + (32*mt*p1p4*p2p3*p2p4)/s + 
     $     (32*mt*p1p3*p2p4**2)/s) + 
     $     ct2*(-8*mt**3*s - 16*mt**3*p1p4 + 16*mt**3*p2p4 - 
     $     16*mt*s*p2p4 - 32*mt*p1p4*p2p4 + 
     $     32*mt*p2p3*p2p4 + (64*mt*p1p4*p2p3*p2p4)/s + 
     $     32*mt*p2p4**2 - (64*mt*p2p3*p2p4**2)/s)

      mat(3,7)=cu1*(8*mt*s*p1p2 - 16*mt**3*p1p4 + 
     $     (64*mt*p1p3*p1p4**2)/s - 16*mt*p1p4*p2p3 + 
     $     16*mt**3*p2p4 - 16*mt*p1p3*p2p4 - 
     $     (64*mt*p1p3*p1p4*p2p4)/s) + 
     $     cu3*(8*mt*s*p1p2 + 16*mt*p1p2*p1p4 - 
     $     16*mt*p1p4*p2p3 - (32*mt*p1p4**2*p2p3)/s - 
     $     16*mt*p1p2*p2p4 - 16*mt*p1p3*p2p4 - 
     $     (32*mt*p1p3*p1p4*p2p4)/s + 
     $     (32*mt*p1p4*p2p3*p2p4)/s + (32*mt*p1p3*p2p4**2)/s) + 
     $     cu2*(-8*mt**3*s - 16*mt**3*p1p4 + 16*mt**3*p2p4 + 
     $     32*mt*p2p3*p2p4 + (64*mt*p1p4*p2p3*p2p4)/s - 
     $     (64*mt*p2p3*p2p4**2)/s)

      mat(1,8)=cs1*(8*mt*s*p1p3 - 16*mt*p1p3**2 + 
     $     16*mt*p1p3*p1p4 + 8*mt*s*p2p3 + 
     $     16*mt*p1p3*p2p3 - 16*mt*p1p4*p2p3) + 
     $     cs2*(8*mt*s*p1p3 + 8*mt*s*p2p3 + 16*mt*p1p3*p2p3 - 
     $     16*mt*p2p3**2 - 16*mt*p1p3*p2p4 + 16*mt*p2p3*p2p4)

      mat(2,8)=ct1*(8*mt**3*s - 16*mt**3*p1p3 + 16*mt*s*p1p3 - 
     $     32*mt*p1p3**2 - 32*mt*p1p3*p1p4 + 
     $     (64*mt*p1p3**2*p1p4)/s + 16*mt**3*p2p3 + 
     $     32*mt*p1p3*p2p3 - (64*mt*p1p3*p1p4*p2p3)/s) + 
     $     ct3*(8*mt**3*s + 16*mt*s*p1p3 + 16*mt*p1p2*p1p3 - 
     $     32*mt*p1p3*p1p4 - 16*mt*p1p2*p2p3 - 
     $     (32*mt*p1p3*p1p4*p2p3)/s + 
     $     (32*mt*p1p4*p2p3**2)/s - (32*mt*p1p3**2*p2p4)/s + 
     $     (32*mt*p1p3*p2p3*p2p4)/s) + 
     $     ct2*(-8*mt*s*p1p2 - 16*mt**3*p1p3 + 16*mt*s*p1p3 + 
     $     16*mt**3*p2p3 + 16*mt*p1p4*p2p3 - 
     $     16*mt*p1p3*p2p4 + 32*mt*p2p3*p2p4 + 
     $     (64*mt*p1p3*p2p3*p2p4)/s - (64*mt*p2p3**2*p2p4)/s)

      mat(3,8)=cu1*(8*mt**3*s - 16*mt**3*p1p3 - 32*mt*p1p3*p1p4 + 
     $     (64*mt*p1p3**2*p1p4)/s + 16*mt**3*p2p3 - 
     $     (64*mt*p1p3*p1p4*p2p3)/s) + 
     $     cu3*(-8*mt*s*p1p2 + 16*mt*p1p2*p1p3 - 
     $     16*mt*p1p2*p2p3 + 16*mt*p1p4*p2p3 - 
     $     (32*mt*p1p3*p1p4*p2p3)/s + 
     $     (32*mt*p1p4*p2p3**2)/s + 16*mt*p1p3*p2p4 - 
     $     (32*mt*p1p3**2*p2p4)/s + (32*mt*p1p3*p2p3*p2p4)/s) + 
     $     cu2*(-8*mt*s*p1p2 - 16*mt**3*p1p3 + 16*mt**3*p2p3 + 
     $     16*mt*p1p4*p2p3 + 16*mt*p1p3*p2p4 + 
     $     (64*mt*p1p3*p2p3*p2p4)/s - (64*mt*p2p3**2*p2p4)/s)

       do i=1,3
          mat(i,9)=0d0
       end do

       mat(1,10)=cs1*(-2*mt**3*s**2+2*mt*s**2*p1p2-4*mt**3*s*p1p3+
     $      4*mt*s*p1p2*p1p3+4*mt**3*s*p1p4-
     $      4*mt*s*p1p2*p1p4-4*mt*s*p1p4*p2p3-
     $      8*mt*p1p3*p1p4*p2p3+8*mt*p1p4**2*p2p3-
     $      4*mt*s*p1p3*p2p4-8*mt*p1p3**2*p2p4+
     $      8*mt*p1p3*p1p4*p2p4+8*mt*s*p2p3*p2p4+
     $      16*mt*p1p3*p2p3*p2p4-16*mt*p1p4*p2p3*p2p4)+ 
     $      cs2*(-2*mt**3*s**2+2*mt*s**2*p1p2+4*mt**3*s*p2p3-
     $      4*mt*s*p1p2*p2p3-4*mt*s*p1p4*p2p3+
     $      8*mt*p1p4*p2p3**2-4*mt**3*s*p2p4+
     $      4*mt*s*p1p2*p2p4-4*mt*s*p1p3*p2p4+
     $      8*mt*s*p2p3*p2p4+8*mt*p1p3*p2p3*p2p4-
     $      8*mt*p1p4*p2p3*p2p4-16*mt*p2p3**2*p2p4-
     $      8*mt*p1p3*p2p4**2+16*mt*p2p3*p2p4**2)

       mat(2,10)=ct1*(8*mt**3*s*p1p2-8*mt*s*p1p2**2+
     $      16*mt*s*p1p2*p1p3-32*mt*p1p2*p1p3*p1p4-
     $      16*mt**3*p1p4*p2p3+32*mt*p1p2*p1p4*p2p3-
     $      32*mt*p1p3*p1p4*p2p3+
     $      (64*mt*p1p3*p1p4**2*p2p3)/s-
     $      (32*mt*p1p4**2*p2p3**2)/s+16*mt*s*p1p2*p2p4-
     $      16*mt**3*p1p3*p2p4+32*mt*p1p2*p1p3*p2p4-
     $      32*mt*p1p3**2*p2p4+
     $      (64*mt*p1p3**2*p1p4*p2p4)/s-
     $      32*mt*p1p4*p2p3*p2p4-
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      32*mt*p1p3*p2p4**2-(32*mt*p1p3**2*p2p4**2)/s)+
     $      ct3*(8*mt**3*s*p1p2-8*mt*s*p1p2**2-8*mt**3*s*p1p3-
     $      16*mt**3*p1p4*p2p3+32*mt*p1p2*p1p4*p2p3-
     $      (32*mt*p1p4**2*p2p3**2)/s+16*mt*s*p1p2*p2p4-
     $      16*mt**3*p1p3*p2p4+32*mt*p1p2*p1p3*p2p4-
     $      32*mt*p1p2*p2p3*p2p4+
     $      32*mt*p1p3*p2p3*p2p4-
     $      32*mt*p1p4*p2p3*p2p4-
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (64*mt*p1p4*p2p3**2*p2p4)/s-
     $      32*mt*p1p3*p2p4**2-(32*mt*p1p3**2*p2p4**2)/s+
     $      (64*mt*p1p3*p2p3*p2p4**2)/s)+
     $      ct2*(-8*mt**5*s+8*mt**3*s*p1p2-8*mt**3*s*p1p3-
     $      16*mt**3*p1p4*p2p3-16*mt**3*s*p2p4+
     $      8*mt*s*p1p2*p2p4-16*mt**3*p1p3*p2p4+
     $      64*mt**3*p2p3*p2p4-32*mt*p1p2*p2p3*p2p4+
     $      32*mt*p1p3*p2p3*p2p4-16*mt*p1p4*p2p3*p2p4+
     $      (64*mt*p1p4*p2p3**2*p2p4)/s-
     $      16*mt*p1p3*p2p4**2+64*mt*p2p3*p2p4**2+
     $      (64*mt*p1p3*p2p3*p2p4**2)/s-
     $      (128*mt*p2p3**2*p2p4**2)/s)

       mat(3,10)=cu1*(4*mt**3*s**2+8*mt**3*s*p1p2-8*mt*s*p1p2**2+
     $      8*mt*s*p1p2*p1p3-8*mt**3*s*p1p4-
     $      32*mt*p1p2*p1p3*p1p4-16*mt**3*p1p4*p2p3+
     $      32*mt*p1p2*p1p4*p2p3-
     $      16*mt*p1p3*p1p4*p2p3+
     $      (64*mt*p1p3*p1p4**2*p2p3)/s-
     $      (32*mt*p1p4**2*p2p3**2)/s+16*mt*s*p1p2*p2p4-
     $      16*mt**3*p1p3*p2p4+32*mt*p1p2*p1p3*p2p4-
     $      16*mt*p1p3**2*p2p4+
     $      (64*mt*p1p3**2*p1p4*p2p4)/s-
     $      16*mt*s*p2p3*p2p4-
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      32*mt*p1p3*p2p4**2-(32*mt*p1p3**2*p2p4**2)/s)+
     $      cu3*(-8*mt**5*s+4*mt**3*s**2+8*mt**3*s*p1p2+
     $      8*mt*s*p1p2*p1p3+32*mt**3*p1p3*p1p4-
     $      16*mt**3*p1p4*p2p3-16*mt*p1p3*p1p4*p2p3-
     $      16*mt**3*s*p2p4-16*mt**3*p1p3*p2p4-
     $      16*mt*p1p3**2*p2p4+32*mt**3*p2p3*p2p4-
     $      16*mt*s*p2p3*p2p4-32*mt*p1p2*p2p3*p2p4-
     $      (128*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (64*mt*p1p4*p2p3**2*p2p4)/s+
     $      64*mt*p2p3*p2p4**2+(64*mt*p1p3*p2p3*p2p4**2)/s)+ 
     $      cu2*(-8*mt**5*s+4*mt**3*s**2+8*mt**3*s*p1p2-8*mt**3*s*p1p3-
     $      8*mt**3*s*p2p3-16*mt**3*p1p4*p2p3-
     $      16*mt**3*s*p2p4-16*mt**3*p1p3*p2p4+
     $      64*mt**3*p2p3*p2p4-16*mt*s*p2p3*p2p4-
     $      32*mt*p1p2*p2p3*p2p4+32*mt*p1p3*p2p3*p2p4+
     $      32*mt*p2p3**2*p2p4+
     $      (64*mt*p1p4*p2p3**2*p2p4)/s+
     $      64*mt*p2p3*p2p4**2+
     $      (64*mt*p1p3*p2p3*p2p4**2)/s-
     $      (128*mt*p2p3**2*p2p4**2)/s)

       mat(1,11)=cs1*(2*mt**3*s**2-2*mt*s**2*p1p2+4*mt**3*s*p1p3-
     $      4*mt*s*p1p2*p1p3-4*mt**3*s*p1p4+
     $      4*mt*s*p1p2*p1p4-8*mt*s*p1p3*p1p4-
     $      16*mt*p1p3**2*p1p4+16*mt*p1p3*p1p4**2+
     $      4*mt*s*p1p4*p2p3+8*mt*p1p3*p1p4*p2p3-
     $      8*mt*p1p4**2*p2p3+4*mt*s*p1p3*p2p4+
     $      8*mt*p1p3**2*p2p4-8*mt*p1p3*p1p4*p2p4)+
     $      cs2*(2*mt**3*s**2-2*mt*s**2*p1p2-8*mt*s*p1p3*p1p4-
     $      4*mt**3*s*p2p3+4*mt*s*p1p2*p2p3+
     $      4*mt*s*p1p4*p2p3+16*mt*p1p3*p1p4*p2p3-
     $      8*mt*p1p4*p2p3**2+4*mt**3*s*p2p4-
     $      4*mt*s*p1p2*p2p4+4*mt*s*p1p3*p2p4-
     $      16*mt*p1p3*p1p4*p2p4-8*mt*p1p3*p2p3*p2p4+
     $      8*mt*p1p4*p2p3*p2p4+8*mt*p1p3*p2p4**2)

       mat(2,11)=ct1*(8*mt**5*s-8*mt**3*s*p1p2+16*mt**3*s*p1p3-
     $      64*mt**3*p1p3*p1p4+32*mt*p1p2*p1p3*p1p4-
     $      64*mt*p1p3**2*p1p4+(128*mt*p1p3**2*p1p4**2)/s+
     $      16*mt**3*p1p4*p2p3-
     $      (64*mt*p1p3*p1p4**2*p2p3)/s+16*mt**3*s*p2p4+
     $      16*mt**3*p1p3*p2p4-64*mt*p1p3*p1p4*p2p4-
     $      (64*mt*p1p3**2*p1p4*p2p4)/s)+
     $      ct3*(8*mt**5*s-8*mt**3*s*p1p2-8*mt*s*p1p2*p1p3-
     $      32*mt**3*p1p3*p1p4+32*mt*p1p2*p1p3*p1p4+
     $      16*mt**3*p1p4*p2p3+16*mt*p1p3*p1p4*p2p3-
     $      (64*mt*p1p3*p1p4**2*p2p3)/s+16*mt**3*s*p2p4+
     $      16*mt**3*p1p3*p2p4+16*mt*p1p3**2*p2p4-
     $      64*mt*p1p3*p1p4*p2p4-
     $      (64*mt*p1p3**2*p1p4*p2p4)/s-
     $      32*mt**3*p2p3*p2p4+
     $      (128*mt*p1p3*p1p4*p2p3*p2p4)/s)+
     $      ct2*(-8*mt**3*s*p1p2+8*mt*s*p1p2**2-
     $      8*mt*s*p1p2*p1p3+16*mt**3*p1p4*p2p3-
     $      32*mt*p1p2*p1p4*p2p3+16*mt*p1p3*p1p4*p2p3+
     $      (32*mt*p1p4**2*p2p3**2)/s+8*mt**3*s*p2p4-
     $      16*mt*s*p1p2*p2p4+16*mt**3*p1p3*p2p4-
     $      32*mt*p1p2*p1p3*p2p4+16*mt*p1p3**2*p2p4-
     $      32*mt*p1p3*p1p4*p2p4+32*mt*p1p2*p2p3*p2p4+
     $      32*mt*p1p4*p2p3*p2p4+
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (64*mt*p1p4*p2p3**2*p2p4)/s+
     $      32*mt*p1p3*p2p4**2+(32*mt*p1p3**2*p2p4**2)/s-
     $      (64*mt*p1p3*p2p3*p2p4**2)/s)

       mat(3,11)=cu1*(8*mt**5*s-8*mt**3*s*p1p2+4*mt*s**2*p1p2+
     $      8*mt**3*s*p1p3-8*mt*s*p1p2*p1p4-
     $      64*mt**3*p1p3*p1p4+32*mt*p1p2*p1p3*p1p4-
     $      32*mt*p1p3**2*p1p4+(128*mt*p1p3**2*p1p4**2)/s+
     $      16*mt**3*p1p4*p2p3-8*mt*s*p1p4*p2p3+
     $      16*mt*p1p4**2*p2p3-
     $      (64*mt*p1p3*p1p4**2*p2p3)/s+16*mt**3*s*p2p4+
     $      16*mt**3*p1p3*p2p4-8*mt*s*p1p3*p2p4-
     $      48*mt*p1p3*p1p4*p2p4-
     $      (64*mt*p1p3**2*p1p4*p2p4)/s)+
     $      cu3*(-8*mt**3*s*p1p2+4*mt*s**2*p1p2+8*mt*s*p1p2**2+
     $      8*mt**3*s*p1p3+32*mt*p1p2*p1p3*p1p4-
     $      32*mt*p1p3**2*p1p4+16*mt**3*p1p4*p2p3-
     $      8*mt*s*p1p4*p2p3-32*mt*p1p2*p1p4*p2p3-
     $      (64*mt*p1p3*p1p4**2*p2p3)/s+
     $      (32*mt*p1p4**2*p2p3**2)/s-16*mt*s*p1p2*p2p4+
     $      16*mt**3*p1p3*p2p4-8*mt*s*p1p3*p2p4-
     $      32*mt*p1p2*p1p3*p2p4-
     $      (64*mt*p1p3**2*p1p4*p2p4)/s+
     $      32*mt*p1p4*p2p3*p2p4+
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      32*mt*p1p3*p2p4**2+(32*mt*p1p3**2*p2p4**2)/s)+
     $      cu2*(-8*mt**3*s*p1p2+4*mt*s**2*p1p2+8*mt*s*p1p2**2-
     $      8*mt*s*p1p2*p1p3-8*mt*s*p1p2*p2p3+
     $      16*mt**3*p1p4*p2p3-8*mt*s*p1p4*p2p3-
     $      32*mt*p1p2*p1p4*p2p3+16*mt*p1p3*p1p4*p2p3+
     $      16*mt*p1p4*p2p3**2+(32*mt*p1p4**2*p2p3**2)/s-
     $      16*mt*s*p1p2*p2p4+16*mt**3*p1p3*p2p4-
     $      8*mt*s*p1p3*p2p4-32*mt*p1p2*p1p3*p2p4+
     $      16*mt*p1p3**2*p2p4+32*mt*p1p2*p2p3*p2p4+
     $      16*mt*p1p3*p2p3*p2p4+32*mt*p1p4*p2p3*p2p4+
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (64*mt*p1p4*p2p3**2*p2p4)/s+
     $      32*mt*p1p3*p2p4**2+(32*mt*p1p3**2*p2p4**2)/s-
     $      (64*mt*p1p3*p2p3*p2p4**2)/s)

       do i=1,3
          mat(i,12)=0d0
       end do

       mat(1,13)=cs1*(-2*mt**3*s**2+2*mt*s**2*p1p2-4*mt**3*s*p1p3+
     $      4*mt*s*p1p2*p1p3+4*mt**3*s*p1p4-
     $      4*mt*s*p1p2*p1p4-4*mt*s*p1p4*p2p3-
     $      8*mt*p1p3*p1p4*p2p3+8*mt*p1p4**2*p2p3-
     $      4*mt*s*p1p3*p2p4-8*mt*p1p3**2*p2p4+
     $      8*mt*p1p3*p1p4*p2p4+8*mt*s*p2p3*p2p4+
     $      16*mt*p1p3*p2p3*p2p4-16*mt*p1p4*p2p3*p2p4)+ 
     $      cs2*(-2*mt**3*s**2+2*mt*s**2*p1p2+4*mt**3*s*p2p3-
     $      4*mt*s*p1p2*p2p3-4*mt*s*p1p4*p2p3+
     $      8*mt*p1p4*p2p3**2-4*mt**3*s*p2p4+
     $      4*mt*s*p1p2*p2p4-4*mt*s*p1p3*p2p4+
     $      8*mt*s*p2p3*p2p4+8*mt*p1p3*p2p3*p2p4-
     $      8*mt*p1p4*p2p3*p2p4-16*mt*p2p3**2*p2p4-
     $      8*mt*p1p3*p2p4**2+16*mt*p2p3*p2p4**2)

       mat(2,13)=ct1*(8*mt**3*s*p1p2-8*mt*s*p1p2**2-
     $      8*mt**3*s*p1p3+16*mt*s*p1p2*p1p3-
     $      32*mt*p1p2*p1p3*p1p4-16*mt**3*p1p4*p2p3+
     $      32*mt*p1p2*p1p4*p2p3-
     $      32*mt*p1p3*p1p4*p2p3+
     $      (64*mt*p1p3*p1p4**2*p2p3)/s-
     $      (32*mt*p1p4**2*p2p3**2)/s+8*mt*s*p1p2*p2p4-
     $      16*mt**3*p1p3*p2p4+32*mt*p1p2*p1p3*p2p4-
     $      32*mt*p1p3**2*p2p4+
     $      (64*mt*p1p3**2*p1p4*p2p4)/s+
     $      32*mt*p1p3*p2p3*p2p4-
     $      16*mt*p1p4*p2p3*p2p4-
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      16*mt*p1p3*p2p4**2-(32*mt*p1p3**2*p2p4**2)/s)+
     $      ct3*(-8*mt**5*s+8*mt**3*s*p1p2-16*mt**3*s*p1p3+
     $      32*mt**3*p1p3*p1p4-16*mt**3*p1p4*p2p3+
     $      8*mt*s*p1p2*p2p4-16*mt**3*p1p3*p2p4+
     $      32*mt**3*p2p3*p2p4-32*mt*p1p2*p2p3*p2p4+
     $      64*mt*p1p3*p2p3*p2p4-
     $      16*mt*p1p4*p2p3*p2p4-
     $      (128*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (64*mt*p1p4*p2p3**2*p2p4)/s-
     $      16*mt*p1p3*p2p4**2+(64*mt*p1p3*p2p3*p2p4**2)/s)+
     $      ct2*(-8*mt**5*s+8*mt**3*s*p1p2-16*mt**3*s*p1p3-
     $      16*mt**3*p1p4*p2p3-16*mt**3*s*p2p4-
     $      16*mt**3*p1p3*p2p4+64*mt**3*p2p3*p2p4-
     $      32*mt*p1p2*p2p3*p2p4+64*mt*p1p3*p2p3*p2p4+
     $      (64*mt*p1p4*p2p3**2*p2p4)/s+
     $      64*mt*p2p3*p2p4**2+
     $      (64*mt*p1p3*p2p3*p2p4**2)/s-
     $      (128*mt*p2p3**2*p2p4**2)/s)

       mat(3,13)=cu1*(8*mt**3*s*p1p2-4*mt*s**2*p1p2-
     $      8*mt*s*p1p2**2+16*mt*s*p1p2*p1p3+
     $      8*mt*s*p1p2*p1p4-32*mt*p1p2*p1p3*p1p4-
     $      16*mt**3*p1p4*p2p3+8*mt*s*p1p4*p2p3+
     $      32*mt*p1p2*p1p4*p2p3-
     $      32*mt*p1p3*p1p4*p2p3-16*mt*p1p4**2*p2p3+
     $      (64*mt*p1p3*p1p4**2*p2p3)/s-
     $      (32*mt*p1p4**2*p2p3**2)/s+8*mt*s*p1p2*p2p4-
     $      16*mt**3*p1p3*p2p4+8*mt*s*p1p3*p2p4+
     $      32*mt*p1p2*p1p3*p2p4-32*mt*p1p3**2*p2p4-
     $      16*mt*p1p3*p1p4*p2p4+
     $      (64*mt*p1p3**2*p1p4*p2p4)/s-
     $      16*mt*p1p4*p2p3*p2p4-
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      16*mt*p1p3*p2p4**2-(32*mt*p1p3**2*p2p4**2)/s)+
     $      cu3*(8*mt**3*s*p1p2-4*mt*s**2*p1p2-8*mt*s*p1p2**2+
     $      16*mt*s*p1p2*p1p3-16*mt**3*p1p4*p2p3+
     $      8*mt*s*p1p4*p2p3+32*mt*p1p2*p1p4*p2p3-
     $      32*mt*p1p3*p1p4*p2p3-
     $      (32*mt*p1p4**2*p2p3**2)/s-8*mt**3*s*p2p4-
     $      16*mt**3*p1p3*p2p4+8*mt*s*p1p3*p2p4+
     $      32*mt*p1p2*p1p3*p2p4-32*mt*p1p3**2*p2p4-
     $      32*mt*p1p2*p2p3*p2p4-
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (64*mt*p1p4*p2p3**2*p2p4)/s-
     $      (32*mt*p1p3**2*p2p4**2)/s+32*mt*p2p3*p2p4**2+
     $      (64*mt*p1p3*p2p3*p2p4**2)/s)+
     $      cu2*(-8*mt**5*s+8*mt**3*s*p1p2-4*mt*s**2*p1p2-
     $      16*mt**3*s*p1p3+8*mt*s*p1p2*p2p3-
     $      16*mt**3*p1p4*p2p3+8*mt*s*p1p4*p2p3-
     $      16*mt*p1p4*p2p3**2-8*mt**3*s*p2p4-
     $      16*mt**3*p1p3*p2p4+8*mt*s*p1p3*p2p4+
     $      64*mt**3*p2p3*p2p4-32*mt*p1p2*p2p3*p2p4+
     $      48*mt*p1p3*p2p3*p2p4+
     $      (64*mt*p1p4*p2p3**2*p2p4)/s+
     $      32*mt*p2p3*p2p4**2+
     $      (64*mt*p1p3*p2p3*p2p4**2)/s-
     $      (128*mt*p2p3**2*p2p4**2)/s)

       mat(1,14)=cs1*(2*mt**3*s**2 - 2*mt*s**2*p1p2 + 4*mt**3*s*p1p3 - 
     $      4*mt*s*p1p2*p1p3 - 4*mt**3*s*p1p4 + 
     $      4*mt*s*p1p2*p1p4 - 8*mt*s*p1p3*p1p4 - 
     $      16*mt*p1p3**2*p1p4 + 16*mt*p1p3*p1p4**2 + 
     $      4*mt*s*p1p4*p2p3 + 8*mt*p1p3*p1p4*p2p3 - 
     $      8*mt*p1p4**2*p2p3 + 4*mt*s*p1p3*p2p4 + 
     $      8*mt*p1p3**2*p2p4 - 8*mt*p1p3*p1p4*p2p4) + 
     $      cs2*(2*mt**3*s**2 - 2*mt*s**2*p1p2 - 8*mt*s*p1p3*p1p4 - 
     $      4*mt**3*s*p2p3 + 4*mt*s*p1p2*p2p3 + 
     $      4*mt*s*p1p4*p2p3 + 16*mt*p1p3*p1p4*p2p3 - 
     $      8*mt*p1p4*p2p3**2 + 4*mt**3*s*p2p4 - 
     $      4*mt*s*p1p2*p2p4 + 4*mt*s*p1p3*p2p4 - 
     $      16*mt*p1p3*p1p4*p2p4 - 8*mt*p1p3*p2p3*p2p4 + 
     $      8*mt*p1p4*p2p3*p2p4 + 8*mt*p1p3*p2p4**2)

       mat(2,14)=ct1*(8*mt**5*s - 8*mt**3*s*p1p2 + 16*mt**3*s*p1p3 - 
     $      8*mt*s*p1p2*p1p3 - 64*mt**3*p1p3*p1p4 + 
     $      32*mt*p1p2*p1p3*p1p4 - 64*mt*p1p3**2*p1p4 + 
     $      (128*mt*p1p3**2*p1p4**2)/s + 16*mt**3*p1p4*p2p3 + 
     $      16*mt*p1p3*p1p4*p2p3 - 
     $      (64*mt*p1p3*p1p4**2*p2p3)/s + 8*mt**3*s*p2p4 + 
     $      16*mt**3*p1p3*p2p4 + 16*mt*p1p3**2*p2p4 - 
     $      32*mt*p1p3*p1p4*p2p4 - 
     $      (64*mt*p1p3**2*p1p4*p2p4)/s) + 
     $      ct3*(-8*mt**3*s*p1p2 + 8*mt*s*p1p2**2 - 
     $      16*mt*s*p1p2*p1p3 + 32*mt*p1p2*p1p3*p1p4 + 
     $      16*mt**3*p1p4*p2p3 - 32*mt*p1p2*p1p4*p2p3 + 
     $      32*mt*p1p3*p1p4*p2p3 - 
     $      (64*mt*p1p3*p1p4**2*p2p3)/s + 
     $      (32*mt*p1p4**2*p2p3**2)/s + 8*mt**3*s*p2p4 + 
     $      16*mt**3*p1p3*p2p4 - 32*mt*p1p2*p1p3*p2p4 + 
     $      32*mt*p1p3**2*p2p4 - 32*mt*p1p3*p1p4*p2p4 - 
     $      (64*mt*p1p3**2*p1p4*p2p4)/s + 
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s + 
     $      (32*mt*p1p3**2*p2p4**2)/s) + 
     $      ct2*(-8*mt**3*s*p1p2 + 8*mt*s*p1p2**2 - 
     $      16*mt*s*p1p2*p1p3 + 16*mt**3*p1p4*p2p3 - 
     $      32*mt*p1p2*p1p4*p2p3 + 32*mt*p1p3*p1p4*p2p3 + 
     $      (32*mt*p1p4**2*p2p3**2)/s - 16*mt*s*p1p2*p2p4 + 
     $      16*mt**3*p1p3*p2p4 - 32*mt*p1p2*p1p3*p2p4 + 
     $      32*mt*p1p3**2*p2p4 + 32*mt*p1p2*p2p3*p2p4 + 
     $      32*mt*p1p4*p2p3*p2p4 + 
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (64*mt*p1p4*p2p3**2*p2p4)/s + 
     $      32*mt*p1p3*p2p4**2 + (32*mt*p1p3**2*p2p4**2)/s - 
     $      (64*mt*p1p3*p2p3*p2p4**2)/s)

       mat(3,14)=cu1*(8*mt**5*s - 4*mt**3*s**2 - 8*mt**3*s*p1p2 + 
     $      16*mt**3*s*p1p3 + 8*mt**3*s*p1p4 - 
     $      64*mt**3*p1p3*p1p4 + 16*mt*s*p1p3*p1p4 + 
     $      32*mt*p1p2*p1p3*p1p4 - 64*mt*p1p3**2*p1p4 - 
     $      32*mt*p1p3*p1p4**2 + (128*mt*p1p3**2*p1p4**2)/s + 
     $      16*mt**3*p1p4*p2p3 - 
     $      (64*mt*p1p3*p1p4**2*p2p3)/s + 8*mt**3*s*p2p4 + 
     $      16*mt**3*p1p3*p2p4 - 32*mt*p1p3*p1p4*p2p4 - 
     $      (64*mt*p1p3**2*p1p4*p2p4)/s) + 
     $      cu3*(8*mt**5*s - 4*mt**3*s**2 - 8*mt**3*s*p1p2 + 
     $      16*mt**3*s*p1p3 - 32*mt**3*p1p3*p1p4 + 
     $      16*mt*s*p1p3*p1p4 + 32*mt*p1p2*p1p3*p1p4 - 
     $      64*mt*p1p3**2*p1p4 + 16*mt**3*p1p4*p2p3 - 
     $      (64*mt*p1p3*p1p4**2*p2p3)/s - 8*mt*s*p1p2*p2p4 + 
     $      16*mt**3*p1p3*p2p4 - 
     $      (64*mt*p1p3**2*p1p4*p2p4)/s - 
     $      32*mt**3*p2p3*p2p4 + 16*mt*p1p4*p2p3*p2p4 + 
     $      (128*mt*p1p3*p1p4*p2p3*p2p4)/s + 
     $      16*mt*p1p3*p2p4**2) + 
     $      cu2*(-4*mt**3*s**2 - 8*mt**3*s*p1p2 + 8*mt*s*p1p2**2 - 
     $      16*mt*s*p1p2*p1p3 + 16*mt*s*p1p3*p1p4 + 
     $      8*mt**3*s*p2p3 + 16*mt**3*p1p4*p2p3 - 
     $      32*mt*p1p2*p1p4*p2p3 + 
     $      (32*mt*p1p4**2*p2p3**2)/s - 8*mt*s*p1p2*p2p4 + 
     $      16*mt**3*p1p3*p2p4 - 32*mt*p1p2*p1p3*p2p4 + 
     $      32*mt*p1p3**2*p2p4 + 32*mt*p1p2*p2p3*p2p4 + 
     $      16*mt*p1p4*p2p3*p2p4 + 
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (64*mt*p1p4*p2p3**2*p2p4)/s + 
     $      16*mt*p1p3*p2p4**2 + (32*mt*p1p3**2*p2p4**2)/s - 
     $      (64*mt*p1p3*p2p3*p2p4**2)/s)

       do i=1,3
          mat(i,15)=0d0
       end do

       mat(1,16)=cs1*(4*mt**2*s*p1p4 + 8*mt**2*p1p3*p1p4 - 
     $      8*mt**2*p1p4**2 + 4*s*p1p2*p2p4 - 
     $      8*p1p2*p1p3*p2p4 + 8*p1p2*p1p4*p2p4 - 
     $      24*p1p4*p2p3*p2p4 - 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s + 
     $      (16*p1p4**2*p2p3*p2p4)/s - 8*p1p3*p2p4**2 + 
     $      (16*p1p3**2*p2p4**2)/s - (16*p1p3*p1p4*p2p4**2)/s) +
     $      cs2*(4*mt**2*s*p1p4 - 8*mt**2*p1p4*p2p3 + 
     $      4*s*p1p2*p2p4 + 8*mt**2*p1p4*p2p4 + 
     $      8*p1p2*p2p3*p2p4 - 24*p1p4*p2p3*p2p4 + 
     $      (16*p1p4*p2p3**2*p2p4)/s - 8*p1p2*p2p4**2 - 
     $      8*p1p3*p2p4**2 - (16*p1p3*p2p3*p2p4**2)/s - 
     $      (16*p1p4*p2p3*p2p4**2)/s + (16*p1p3*p2p4**3)/s)
       mat(2,16)=ct1*(-4*mt**4*s - 4*mt**2*s*p1p2 + 16*p1p2**2*p1p4 + 
     $      16*mt**2*p1p3*p1p4 + 8*mt**2*p1p4*p2p3 - 
     $      (64*p1p2*p1p4**2*p2p3)/s + 
     $      (64*p1p4**3*p2p3**2)/s**2 - 16*mt**2*p1p2*p2p4 + 
     $      8*s*p1p2*p2p4 + 8*mt**2*p1p3*p2p4 - 
     $      16*p1p2*p1p3*p2p4 - 32*p1p2*p1p4*p2p4 + 
     $      16*mt**2*p2p3*p2p4 - 16*p1p4*p2p3*p2p4 + 
     $      (32*mt**2*p1p4*p2p3*p2p4)/s - 
     $      (32*p1p3*p1p4*p2p3*p2p4)/s + 
     $      (64*p1p4**2*p2p3*p2p4)/s - 16*p1p3*p2p4**2 + 
     $      (32*mt**2*p1p3*p2p4**2)/s + (32*p1p3**2*p2p4**2)/s + 
     $      (64*p1p3*p1p4*p2p4**2)/s - 
     $      (64*p1p3**2*p1p4*p2p4**2)/s**2) + 
     $      ct3*(4*mt**4*s + 4*mt**2*s*p1p2 - 16*mt**2*p1p2*p1p4 - 
     $      8*mt**2*p1p4*p2p3 + (32*mt**2*p1p4**2*p2p3)/s + 
     $      8*s*p1p2*p2p4 + 16*p1p2**2*p2p4 - 
     $      8*mt**2*p1p3*p2p4 - 32*p1p2*p1p4*p2p4 + 
     $      (32*mt**2*p1p3*p1p4*p2p4)/s - 
     $      16*mt**2*p2p3*p2p4 - 16*p1p2*p2p3*p2p4 - 
     $      16*p1p4*p2p3*p2p4 + 
     $      (64*p1p4**2*p2p3*p2p4)/s + 
     $      (32*p1p4*p2p3**2*p2p4)/s - 
     $      (64*p1p4**2*p2p3**2*p2p4)/s**2 - 
     $      16*p1p3*p2p4**2 - (64*p1p2*p1p3*p2p4**2)/s + 
     $      (64*p1p3*p1p4*p2p4**2)/s + 
     $      (32*p1p3*p2p3*p2p4**2)/s + 
     $      (64*p1p3**2*p2p4**3)/s**2) + 
     $      ct2*(4*mt**4*s + 4*mt**2*s*p1p2 + 16*mt**4*p1p4 - 
     $      8*mt**2*p1p4*p2p3 - 16*mt**2*p1p2*p2p4 + 
     $      8*s*p1p2*p2p4 - 8*mt**2*p1p3*p2p4 + 
     $      32*mt**2*p1p4*p2p4 - 16*mt**2*p2p3*p2p4 - 
     $      16*p1p2*p2p3*p2p4 - 16*p1p4*p2p3*p2p4 - 
     $      (96*mt**2*p1p4*p2p3*p2p4)/s + 
     $      (32*p1p4*p2p3**2*p2p4)/s - 16*p1p2*p2p4**2 - 
     $      16*p1p3*p2p4**2 + (32*mt**2*p1p3*p2p4**2)/s + 
     $      (64*p1p2*p2p3*p2p4**2)/s + 
     $      (32*p1p3*p2p3*p2p4**2)/s - 
     $      (96*p1p4*p2p3*p2p4**2)/s + 
     $      (128*p1p4*p2p3**2*p2p4**2)/s**2 + 
     $      (32*p1p3*p2p4**3)/s - (128*p1p3*p2p3*p2p4**3)/s**2)
       mat(3,16)=cu1*(-4*mt**4*s - 4*mt**2*s*p1p2 - 8*mt**2*s*p1p4 + 
     $      16*p1p2**2*p1p4 + 16*mt**2*p1p3*p1p4 + 
     $      16*mt**2*p1p4**2 + 8*mt**2*p1p4*p2p3 - 
     $      (64*p1p2*p1p4**2*p2p3)/s + 
     $      (64*p1p4**3*p2p3**2)/s**2 - 16*mt**2*p1p2*p2p4 + 
     $      8*mt**2*p1p3*p2p4 - 32*p1p2*p1p4*p2p4 + 
     $      16*mt**2*p2p3*p2p4 + 32*p1p4*p2p3*p2p4 + 
     $      (32*mt**2*p1p4*p2p3*p2p4)/s - 
     $      (64*p1p3*p1p4*p2p3*p2p4)/s + 
     $      (32*mt**2*p1p3*p2p4**2)/s + 
     $      (64*p1p3*p1p4*p2p4**2)/s - 
     $      (64*p1p3**2*p1p4*p2p4**2)/s**2) + 
     $      cu3*(-4*mt**4*s - 4*mt**2*s*p1p2 - 8*mt**2*s*p1p4 - 
     $      16*mt**2*p1p2*p1p4 + 16*mt**2*p1p3*p1p4 + 
     $      8*mt**2*p1p4*p2p3 + (32*mt**2*p1p4**2*p2p3)/s + 
     $      16*mt**4*p2p4 + 8*mt**2*p1p3*p2p4 + 
     $      32*mt**2*p1p4*p2p4 - 
     $      (32*mt**2*p1p3*p1p4*p2p4)/s + 
     $      16*mt**2*p2p3*p2p4 + 32*p1p4*p2p3*p2p4 + 
     $      (64*p1p2*p1p4*p2p3*p2p4)/s - 
     $      (64*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (128*p1p4**2*p2p3**2*p2p4)/s**2 - 
     $      (64*mt**2*p2p3*p2p4**2)/s - 
     $      (128*p1p4*p2p3*p2p4**2)/s + 
     $      (128*p1p3*p1p4*p2p3*p2p4**2)/s**2) + 
     $      cu2*(4*mt**4*s + 4*mt**2*s*p1p2 + 16*mt**4*p1p4 - 
     $      8*mt**2*s*p1p4 + 8*mt**2*p1p4*p2p3 - 
     $      16*mt**2*p1p2*p2p4 - 8*mt**2*p1p3*p2p4 + 
     $      32*mt**2*p1p4*p2p4 - 16*mt**2*p2p3*p2p4 - 
     $      16*p1p2*p2p3*p2p4 + 32*p1p4*p2p3*p2p4 - 
     $      (96*mt**2*p1p4*p2p3*p2p4)/s - 
     $      (32*p1p4*p2p3**2*p2p4)/s + 
     $      (32*mt**2*p1p3*p2p4**2)/s + 
     $      (64*p1p2*p2p3*p2p4**2)/s + 
     $      (32*p1p3*p2p3*p2p4**2)/s - 
     $      (128*p1p4*p2p3*p2p4**2)/s + 
     $      (128*p1p4*p2p3**2*p2p4**2)/s**2 - 
     $      (128*p1p3*p2p3*p2p4**3)/s**2)

       mat(1,17)=cs1*(4*s*p1p2*p1p4 + 8*p1p2*p1p3*p1p4 - 
     $      8*p1p2*p1p4**2 - 8*p1p4**2*p2p3 - 
     $      (16*p1p3*p1p4**2*p2p3)/s + (16*p1p4**3*p2p3)/s + 
     $      4*mt**2*s*p2p4 - 8*mt**2*p1p3*p2p4 + 
     $      8*mt**2*p1p4*p2p4 - 24*p1p3*p1p4*p2p4 + 
     $      (16*p1p3**2*p1p4*p2p4)/s - 
     $      (16*p1p3*p1p4**2*p2p4)/s) + 
     $      cs2*(4*s*p1p2*p1p4 - 8*p1p2*p1p4*p2p3 - 
     $      8*p1p4**2*p2p3 + (16*p1p4**2*p2p3**2)/s + 
     $      4*mt**2*s*p2p4 + 8*p1p2*p1p4*p2p4 - 
     $      24*p1p3*p1p4*p2p4 + 8*mt**2*p2p3*p2p4 - 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (16*p1p4**2*p2p3*p2p4)/s - 8*mt**2*p2p4**2 + 
     $      (16*p1p3*p1p4*p2p4**2)/s)
       mat(2,17)=ct1*(-4*mt**4*s - 4*mt**2*s*p1p2 + 
     $      16*mt**2*p1p2*p1p4 + 16*mt**2*p1p3*p1p4 + 
     $      16*p1p2*p1p3*p1p4 - 
     $      (64*p1p2*p1p3*p1p4**2)/s + 8*mt**2*p1p4*p2p3 - 
     $      (32*mt**2*p1p4**2*p2p3)/s - 
     $      (32*p1p3*p1p4**2*p2p3)/s + 
     $      (128*p1p3*p1p4**3*p2p3)/s**2 - 16*mt**4*p2p4 + 
     $      8*mt**2*s*p2p4 - 8*mt**2*p1p3*p2p4 - 
     $      32*mt**2*p1p4*p2p4 - 32*p1p3*p1p4*p2p4 + 
     $      (96*mt**2*p1p3*p1p4*p2p4)/s + 
     $      (32*p1p3**2*p1p4*p2p4)/s + 
     $      (128*p1p3*p1p4**2*p2p4)/s - 
     $      (128*p1p3**2*p1p4**2*p2p4)/s**2) + 
     $      ct3*(4*mt**4*s + 4*mt**2*s*p1p2 - 16*mt**4*p1p4 - 
     $      16*mt**2*p1p3*p1p4 + (64*mt**2*p1p3*p1p4**2)/s - 
     $      8*mt**2*p1p4*p2p3 + 8*mt**2*s*p2p4 + 
     $      16*mt**2*p1p2*p2p4 - 8*mt**2*p1p3*p2p4 - 
     $      32*mt**2*p1p4*p2p4 - 32*p1p3*p1p4*p2p4 - 
     $      (64*p1p2*p1p3*p1p4*p2p4)/s + 
     $      (128*p1p3*p1p4**2*p2p4)/s - 16*mt**2*p2p3*p2p4 + 
     $      (32*mt**2*p1p4*p2p3*p2p4)/s + 
     $      (64*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (128*p1p3*p1p4**2*p2p3*p2p4)/s**2 - 
     $      (32*mt**2*p1p3*p2p4**2)/s + 
     $      (128*p1p3**2*p1p4*p2p4**2)/s**2) + 
     $      ct2*(4*mt**4*s + 4*mt**2*s*p1p2 + 16*mt**2*p1p2*p1p4 - 
     $      16*mt**2*p1p3*p1p4 - 8*mt**2*p1p4*p2p3 - 
     $      (32*mt**2*p1p4**2*p2p3)/s + 8*mt**2*s*p2p4 - 
     $      16*p1p2**2*p2p4 - 8*mt**2*p1p3*p2p4 + 
     $      32*p1p2*p1p4*p2p4 - 32*p1p3*p1p4*p2p4 - 
     $      (32*mt**2*p1p3*p1p4*p2p4)/s - 
     $      16*mt**2*p2p3*p2p4 + 
     $      (64*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (64*p1p4**2*p2p3*p2p4)/s + 
     $      (64*p1p4**2*p2p3**2*p2p4)/s**2 - 16*mt**2*p2p4**2 + 
     $      (64*p1p2*p1p3*p2p4**2)/s - 
     $      (64*p1p3**2*p2p4**3)/s**2)
       mat(3,17)=cu1*(-4*mt**4*s - 4*mt**2*s*p1p2 + 
     $      16*mt**2*p1p2*p1p4 - 8*s*p1p2*p1p4 + 
     $      16*mt**2*p1p3*p1p4 + 16*p1p2*p1p3*p1p4 + 
     $      16*p1p2*p1p4**2 - (64*p1p2*p1p3*p1p4**2)/s + 
     $      8*mt**2*p1p4*p2p3 + 16*p1p4**2*p2p3 - 
     $      (32*mt**2*p1p4**2*p2p3)/s - 
     $      (32*p1p3*p1p4**2*p2p3)/s - (32*p1p4**3*p2p3)/s + 
     $      (128*p1p3*p1p4**3*p2p3)/s**2 - 16*mt**4*p2p4 + 
     $      8*mt**2*p1p3*p2p4 - 32*mt**2*p1p4*p2p4 + 
     $      16*p1p3*p1p4*p2p4 + 
     $      (96*mt**2*p1p3*p1p4*p2p4)/s - 
     $      (32*p1p3**2*p1p4*p2p4)/s + 
     $      (96*p1p3*p1p4**2*p2p4)/s - 
     $      (128*p1p3**2*p1p4**2*p2p4)/s**2) + 
     $      cu3*(-4*mt**4*s - 4*mt**2*s*p1p2 - 8*s*p1p2*p1p4 - 
     $      16*p1p2**2*p1p4 + 16*mt**2*p1p3*p1p4 + 
     $      16*p1p2*p1p3*p1p4 + 8*mt**2*p1p4*p2p3 + 
     $      16*p1p4**2*p2p3 + (64*p1p2*p1p4**2*p2p3)/s - 
     $      (32*p1p3*p1p4**2*p2p3)/s - 
     $      (64*p1p4**3*p2p3**2)/s**2 + 16*mt**2*p1p2*p2p4 + 
     $      8*mt**2*p1p3*p2p4 + 32*p1p2*p1p4*p2p4 + 
     $      16*p1p3*p1p4*p2p4 - 
     $      (32*p1p3**2*p1p4*p2p4)/s - 
     $      (32*mt**2*p1p4*p2p3*p2p4)/s - 
     $      (64*p1p4**2*p2p3*p2p4)/s - 
     $      (32*mt**2*p1p3*p2p4**2)/s - 
     $      (64*p1p3*p1p4*p2p4**2)/s + 
     $      (64*p1p3**2*p1p4*p2p4**2)/s**2) + 
     $      cu2*(4*mt**4*s + 4*mt**2*s*p1p2 + 16*mt**2*p1p2*p1p4 - 
     $      8*s*p1p2*p1p4 - 16*mt**2*p1p3*p1p4 - 
     $      8*mt**2*p1p4*p2p3 + 16*p1p2*p1p4*p2p3 + 
     $      16*p1p4**2*p2p3 - (32*mt**2*p1p4**2*p2p3)/s - 
     $      (32*p1p4**2*p2p3**2)/s - 16*p1p2**2*p2p4 - 
     $      8*mt**2*p1p3*p2p4 + 32*p1p2*p1p4*p2p4 + 
     $      16*p1p3*p1p4*p2p4 - 
     $      (32*mt**2*p1p3*p1p4*p2p4)/s - 
     $      16*mt**2*p2p3*p2p4 + 
     $      (32*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (64*p1p4**2*p2p3*p2p4)/s + 
     $      (64*p1p4**2*p2p3**2*p2p4)/s**2 + 
     $      (64*p1p2*p1p3*p2p4**2)/s - 
     $      (64*p1p3*p1p4*p2p4**2)/s - 
     $      (64*p1p3**2*p2p4**3)/s**2)

       do i=1,3
          mat(i,18)=0d0
       enddo

      mat(1,19)=cs1*(-4*mt**2*s*p1p4 - 8*mt**2*p1p3*p1p4 + 
     $     8*mt**2*p1p4**2 - 4*s*p1p2*p2p4 + 
     $     8*p1p2*p1p3*p2p4 - 8*p1p2*p1p4*p2p4 + 
     $     24*p1p4*p2p3*p2p4 + 
     $     (16*p1p3*p1p4*p2p3*p2p4)/s - 
     $     (16*p1p4**2*p2p3*p2p4)/s + 8*p1p3*p2p4**2 - 
     $     (16*p1p3**2*p2p4**2)/s + (16*p1p3*p1p4*p2p4**2)/s) +
     $     cs2*(-4*mt**2*s*p1p4 + 8*mt**2*p1p4*p2p3 - 
     $     4*s*p1p2*p2p4 - 8*mt**2*p1p4*p2p4 - 
     $     8*p1p2*p2p3*p2p4 + 24*p1p4*p2p3*p2p4 - 
     $     (16*p1p4*p2p3**2*p2p4)/s + 8*p1p2*p2p4**2 + 
     $     8*p1p3*p2p4**2 + (16*p1p3*p2p3*p2p4**2)/s + 
     $     (16*p1p4*p2p3*p2p4**2)/s - (16*p1p3*p2p4**3)/s)
      mat(2,19)=ct1*(8*mt**2*s*p1p2 - 16*p1p2**2*p1p4 - 
     $     16*mt**2*p1p3*p1p4 - 16*mt**2*p1p4*p2p3 + 
     $     (64*p1p2*p1p4**2*p2p3)/s - 
     $     (64*p1p4**3*p2p3**2)/s**2 + 16*mt**2*p1p2*p2p4 - 
     $     16*mt**2*p1p3*p2p4 + 16*p1p2*p1p4*p2p4 - 
     $     (32*mt**2*p1p4*p2p3*p2p4)/s + 
     $     (64*p1p3*p1p4*p2p3*p2p4)/s - 
     $     (32*p1p4**2*p2p3*p2p4)/s - 
     $     (32*mt**2*p1p3*p2p4**2)/s - 
     $     (32*p1p3*p1p4*p2p4**2)/s + 
     $     (64*p1p3**2*p1p4*p2p4**2)/s**2) + 
     $     ct3*(-8*mt**4*s + 16*mt**2*p1p2*p1p4 - 
     $     (32*mt**2*p1p4**2*p2p3)/s - 16*mt**4*p2p4 + 
     $     16*p1p2*p1p4*p2p4 + 
     $     (32*mt**2*p1p3*p1p4*p2p4)/s + 
     $     32*mt**2*p2p3*p2p4 - 
     $     (64*p1p2*p1p4*p2p3*p2p4)/s - 
     $     (32*p1p4**2*p2p3*p2p4)/s + 
     $     (128*p1p4**2*p2p3**2*p2p4)/s**2 - 
     $     (32*p1p3*p1p4*p2p4**2)/s + 
     $     (64*mt**2*p2p3*p2p4**2)/s - 
     $     (128*p1p3*p1p4*p2p3*p2p4**2)/s**2) + 
     $     ct2*(-8*mt**4*s - 16*mt**4*p1p4 + 16*mt**2*p1p2*p2p4 - 
     $     32*mt**2*p1p4*p2p4 + 32*mt**2*p2p3*p2p4 + 
     $     (96*mt**2*p1p4*p2p3*p2p4)/s - 
     $     (32*mt**2*p1p3*p2p4**2)/s - 
     $     (64*p1p2*p2p3*p2p4**2)/s + 
     $     (128*p1p4*p2p3*p2p4**2)/s - 
     $     (128*p1p4*p2p3**2*p2p4**2)/s**2 + 
     $     (128*p1p3*p2p3*p2p4**3)/s**2)
      mat(3,19)=cu1*(8*mt**2*s*p1p2 - 16*p1p2**2*p1p4 - 
     $     16*mt**2*p1p4*p2p3 + (64*p1p2*p1p4**2*p2p3)/s - 
     $     (64*p1p4**3*p2p3**2)/s**2 + 16*mt**2*p1p2*p2p4 - 
     $     16*mt**2*p1p3*p2p4 + 32*p1p2*p1p4*p2p4 - 
     $     (32*mt**2*p1p4*p2p3*p2p4)/s - 
     $     (64*p1p4**2*p2p3*p2p4)/s - 
     $     (32*mt**2*p1p3*p2p4**2)/s - 
     $     (64*p1p3*p1p4*p2p4**2)/s + 
     $     (64*p1p3**2*p1p4*p2p4**2)/s**2) + 
     $     cu3*(8*mt**2*s*p1p2 + 16*mt**2*p1p2*p1p4 - 
     $     16*mt**2*p1p4*p2p3 - (32*mt**2*p1p4**2*p2p3)/s - 
     $     16*p1p2**2*p2p4 - 16*mt**2*p1p3*p2p4 - 
     $     16*mt**2*p1p4*p2p4 - 
     $     (32*mt**2*p1p3*p1p4*p2p4)/s + 
     $     (64*p1p4**2*p2p3**2*p2p4)/s**2 + 
     $     (64*p1p2*p1p3*p2p4**2)/s + 
     $     (64*p1p4*p2p3*p2p4**2)/s - 
     $     (64*p1p3**2*p2p4**3)/s**2) + 
     $     cu2*(-8*mt**4*s - 16*mt**4*p1p4 + 16*mt**2*p1p2*p2p4 - 
     $     16*mt**2*p1p4*p2p4 + 32*mt**2*p2p3*p2p4 + 
     $     16*p1p2*p2p3*p2p4 + 
     $     (96*mt**2*p1p4*p2p3*p2p4)/s - 
     $     (32*p1p4*p2p3**2*p2p4)/s - 
     $     (32*mt**2*p1p3*p2p4**2)/s - 
     $     (64*p1p2*p2p3*p2p4**2)/s - 
     $     (32*p1p3*p2p3*p2p4**2)/s + 
     $     (64*p1p4*p2p3*p2p4**2)/s - 
     $     (128*p1p4*p2p3**2*p2p4**2)/s**2 + 
     $     (128*p1p3*p2p3*p2p4**3)/s**2)

      mat(1,20)=cs1*(-4*s*p1p2*p1p4 - 8*p1p2*p1p3*p1p4 + 
     $      8*p1p2*p1p4**2 + 8*p1p4**2*p2p3 + 
     $      (16*p1p3*p1p4**2*p2p3)/s - (16*p1p4**3*p2p3)/s - 
     $      4*mt**2*s*p2p4 + 8*mt**2*p1p3*p2p4 - 
     $      8*mt**2*p1p4*p2p4 + 24*p1p3*p1p4*p2p4 - 
     $      (16*p1p3**2*p1p4*p2p4)/s + 
     $      (16*p1p3*p1p4**2*p2p4)/s) + 
     $      cs2*(-4*s*p1p2*p1p4 + 8*p1p2*p1p4*p2p3 + 
     $      8*p1p4**2*p2p3 - (16*p1p4**2*p2p3**2)/s - 
     $      4*mt**2*s*p2p4 - 8*p1p2*p1p4*p2p4 + 
     $      24*p1p3*p1p4*p2p4 - 8*mt**2*p2p3*p2p4 + 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s + 
     $      (16*p1p4**2*p2p3*p2p4)/s + 8*mt**2*p2p4**2 - 
     $      (16*p1p3*p1p4*p2p4**2)/s)
      mat(2,20)=ct1*(8*mt**4*s - 16*mt**2*p1p2*p1p4 - 
     $     32*mt**2*p1p3*p1p4 - 16*p1p2*p1p3*p1p4 + 
     $     (64*p1p2*p1p3*p1p4**2)/s + 
     $     (32*mt**2*p1p4**2*p2p3)/s + 
     $     (32*p1p3*p1p4**2*p2p3)/s - 
     $     (128*p1p3*p1p4**3*p2p3)/s**2 + 16*mt**4*p2p4 + 
     $     16*mt**2*p1p4*p2p4 - 
     $     (96*mt**2*p1p3*p1p4*p2p4)/s + 
     $     (32*p1p3**2*p1p4*p2p4)/s - 
     $     (64*p1p3*p1p4**2*p2p4)/s + 
     $     (128*p1p3**2*p1p4**2*p2p4)/s**2) + 
     $     ct3*(-8*mt**2*s*p1p2 + 16*p1p2**2*p1p4 + 
     $     16*mt**2*p1p4*p2p3 - (64*p1p2*p1p4**2*p2p3)/s + 
     $     (64*p1p4**3*p2p3**2)/s**2 - 16*mt**2*p1p2*p2p4 + 
     $     16*mt**2*p1p3*p2p4 + 16*mt**2*p1p4*p2p4 - 
     $     (64*p1p3*p1p4**2*p2p4)/s + 
     $     (32*mt**2*p1p4*p2p3*p2p4)/s + 
     $     (32*mt**2*p1p3*p2p4**2)/s - 
     $     (64*p1p3**2*p1p4*p2p4**2)/s**2) + 
     $     ct2*(-8*mt**2*s*p1p2 - 16*mt**2*p1p2*p1p4 + 
     $     16*mt**2*p1p4*p2p3 + (32*mt**2*p1p4**2*p2p3)/s + 
     $     16*p1p2**2*p2p4 + 16*mt**2*p1p3*p2p4 - 
     $     32*p1p2*p1p4*p2p4 + 
     $     (32*mt**2*p1p3*p1p4*p2p4)/s + 
     $     (64*p1p4**2*p2p3*p2p4)/s - 
     $     (64*p1p4**2*p2p3**2*p2p4)/s**2 - 
     $     (64*p1p2*p1p3*p2p4**2)/s + 
     $     (64*p1p3*p1p4*p2p4**2)/s + 
     $     (64*p1p3**2*p2p4**3)/s**2)
      mat(3,20)=cu1*(8*mt**4*s - 16*mt**2*p1p2*p1p4 - 
     $     32*mt**2*p1p3*p1p4 + (64*p1p2*p1p3*p1p4**2)/s + 
     $     (32*mt**2*p1p4**2*p2p3)/s - 
     $     (128*p1p3*p1p4**3*p2p3)/s**2 + 16*mt**4*p2p4 + 
     $     32*mt**2*p1p4*p2p4 - 
     $     (96*mt**2*p1p3*p1p4*p2p4)/s - 
     $     (128*p1p3*p1p4**2*p2p4)/s + 
     $     (128*p1p3**2*p1p4**2*p2p4)/s**2) + 
     $     cu3*(8*mt**4*s + 16*mt**4*p1p4 - 32*mt**2*p1p3*p1p4 - 
     $     (64*mt**2*p1p3*p1p4**2)/s - 16*mt**2*p1p2*p2p4 - 
     $     16*p1p2*p1p4*p2p4 + 
     $     (64*p1p2*p1p3*p1p4*p2p4)/s - 
     $     (32*mt**2*p1p4*p2p3*p2p4)/s + 
     $     (32*p1p4**2*p2p3*p2p4)/s + 
     $     (128*p1p3*p1p4**2*p2p3*p2p4)/s**2 + 
     $     (32*mt**2*p1p3*p2p4**2)/s + 
     $     (32*p1p3*p1p4*p2p4**2)/s - 
     $     (128*p1p3**2*p1p4*p2p4**2)/s**2) + 
     $     cu2*(-8*mt**2*s*p1p2 - 16*mt**2*p1p2*p1p4 + 
     $     16*mt**2*p1p4*p2p3 + (32*mt**2*p1p4**2*p2p3)/s + 
     $     16*p1p2**2*p2p4 + 16*mt**2*p1p3*p2p4 - 
     $     16*p1p2*p1p4*p2p4 + 
     $     (32*mt**2*p1p3*p1p4*p2p4)/s + 
     $     16*mt**2*p2p3*p2p4 - 
     $     (64*p1p3*p1p4*p2p3*p2p4)/s + 
     $     (32*p1p4**2*p2p3*p2p4)/s - 
     $     (64*p1p4**2*p2p3**2*p2p4)/s**2 - 
     $     (64*p1p2*p1p3*p2p4**2)/s + 
     $     (32*p1p3*p1p4*p2p4**2)/s + 
     $     (64*p1p3**2*p2p4**3)/s**2)

       do i=1,3
          mat(i,21)=0d0
       end do

       mat(1,22)=cs1*(-4*mt**2*s*p1p3 + 8*mt**2*p1p3**2 - 
     $      8*mt**2*p1p3*p1p4 - 4*s*p1p2*p2p3 - 
     $      8*p1p2*p1p3*p2p3 + 8*p1p2*p1p4*p2p3 + 
     $      8*p1p4*p2p3**2 + (16*p1p3*p1p4*p2p3**2)/s - 
     $      (16*p1p4**2*p2p3**2)/s + 24*p1p3*p2p3*p2p4 - 
     $      (16*p1p3**2*p2p3*p2p4)/s + 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s) + 
     $      cs2*(-4*mt**2*s*p1p3 - 4*s*p1p2*p2p3 - 
     $      8*mt**2*p1p3*p2p3 + 8*p1p2*p2p3**2 + 
     $      8*p1p4*p2p3**2 - (16*p1p4*p2p3**3)/s + 
     $      8*mt**2*p1p3*p2p4 - 8*p1p2*p2p3*p2p4 + 
     $      24*p1p3*p2p3*p2p4 + 
     $      (16*p1p3*p2p3**2*p2p4)/s + 
     $      (16*p1p4*p2p3**2*p2p4)/s - 
     $      (16*p1p3*p2p3*p2p4**2)/s)
       mat(2,22)=ct1*(-8*mt**2*s*p1p2 + 16*p1p2**2*p1p3 - 
     $      16*mt**2*p1p2*p2p3 - 32*p1p2*p1p3*p2p3 + 
     $      16*mt**2*p1p4*p2p3 + (32*mt**2*p1p4*p2p3**2)/s + 
     $      (64*p1p3*p1p4*p2p3**2)/s - 
     $      (64*p1p3*p1p4**2*p2p3**2)/s**2 + 
     $      16*mt**2*p1p3*p2p4 - (64*p1p2*p1p3**2*p2p4)/s + 
     $      (32*mt**2*p1p3*p2p3*p2p4)/s + 
     $      (64*p1p3**2*p2p3*p2p4)/s + 
     $      (64*p1p3**3*p2p4**2)/s**2) + 
     $      ct3*(-8*mt**2*s*p1p2 - 16*mt**2*p1p2*p1p3 + 
     $      16*p1p2**2*p2p3 + 16*mt**2*p1p3*p2p3 + 
     $      16*mt**2*p1p4*p2p3 + 
     $      (32*mt**2*p1p3*p1p4*p2p3)/s - 
     $      (64*p1p2*p1p4*p2p3**2)/s + 
     $      (64*p1p4**2*p2p3**3)/s**2 + 16*mt**2*p1p3*p2p4 + 
     $      (32*mt**2*p1p3**2*p2p4)/s - 
     $      (64*p1p3*p2p3**2*p2p4)/s - 
     $      (64*p1p3**2*p2p3*p2p4**2)/s**2) + 
     $      ct2*(8*mt**4*s + 16*mt**4*p1p3 - 16*mt**2*p1p2*p2p3 + 
     $      16*mt**2*p1p3*p2p3 + (32*mt**2*p1p4*p2p3**2)/s - 
     $      32*mt**2*p2p3*p2p4 - 16*p1p2*p2p3*p2p4 - 
     $      (96*mt**2*p1p3*p2p3*p2p4)/s + 
     $      (64*p1p2*p2p3**2*p2p4)/s - 
     $      (64*p1p3*p2p3**2*p2p4)/s + 
     $      (32*p1p4*p2p3**2*p2p4)/s - 
     $      (128*p1p4*p2p3**3*p2p4)/s**2 + 
     $      (32*p1p3*p2p3*p2p4**2)/s + 
     $      (128*p1p3*p2p3**2*p2p4**2)/s**2)
       mat(3,22)=cu1*(-8*mt**2*s*p1p2 + 16*p1p2**2*p1p3 + 
     $      16*mt**2*p1p3*p1p4 - 16*mt**2*p1p2*p2p3 - 
     $      16*p1p2*p1p3*p2p3 + 16*mt**2*p1p4*p2p3 + 
     $      (32*mt**2*p1p4*p2p3**2)/s + 
     $      (32*p1p3*p1p4*p2p3**2)/s - 
     $      (64*p1p3*p1p4**2*p2p3**2)/s**2 + 
     $      16*mt**2*p1p3*p2p4 - (64*p1p2*p1p3**2*p2p4)/s + 
     $      (32*mt**2*p1p3*p2p3*p2p4)/s + 
     $      (32*p1p3**2*p2p3*p2p4)/s - 
     $      (64*p1p3*p1p4*p2p3*p2p4)/s + 
     $      (64*p1p3**3*p2p4**2)/s**2) + 
     $      cu3*(8*mt**4*s - 16*mt**2*p1p2*p1p3 + 16*mt**4*p2p3 - 
     $      16*p1p2*p1p3*p2p3 - 
     $      (32*mt**2*p1p3*p1p4*p2p3)/s + 
     $      (32*p1p3*p1p4*p2p3**2)/s + 
     $      (32*mt**2*p1p3**2*p2p4)/s - 32*mt**2*p2p3*p2p4 + 
     $      (64*p1p2*p1p3*p2p3*p2p4)/s + 
     $      (32*p1p3**2*p2p3*p2p4)/s - 
     $      (64*mt**2*p2p3**2*p2p4)/s + 
     $      (128*p1p3*p1p4*p2p3**2*p2p4)/s**2 - 
     $      (128*p1p3**2*p2p3*p2p4**2)/s**2) + 
     $      cu2*(8*mt**4*s + 16*mt**4*p1p3 - 16*mt**2*p1p2*p2p3 + 
     $      32*mt**2*p1p3*p2p3 + (32*mt**2*p1p4*p2p3**2)/s - 
     $      32*mt**2*p2p3*p2p4 - 
     $      (96*mt**2*p1p3*p2p3*p2p4)/s + 
     $      (64*p1p2*p2p3**2*p2p4)/s - 
     $      (128*p1p3*p2p3**2*p2p4)/s - 
     $      (128*p1p4*p2p3**3*p2p4)/s**2 + 
     $      (128*p1p3*p2p3**2*p2p4**2)/s**2)

       mat(1,23)=cs1*(-4*s*p1p2*p1p3 + 8*p1p2*p1p3**2 - 
     $      8*p1p2*p1p3*p1p4 - 4*mt**2*s*p2p3 - 
     $      8*mt**2*p1p3*p2p3 + 8*mt**2*p1p4*p2p3 + 
     $      24*p1p3*p1p4*p2p3 + 
     $      (16*p1p3**2*p1p4*p2p3)/s - 
     $      (16*p1p3*p1p4**2*p2p3)/s + 8*p1p3**2*p2p4 - 
     $      (16*p1p3**3*p2p4)/s + (16*p1p3**2*p1p4*p2p4)/s) + 
     $      cs2*(-4*s*p1p2*p1p3 - 4*mt**2*s*p2p3 - 
     $      8*p1p2*p1p3*p2p3 + 24*p1p3*p1p4*p2p3 + 
     $      8*mt**2*p2p3**2 - (16*p1p3*p1p4*p2p3**2)/s + 
     $      8*p1p2*p1p3*p2p4 + 8*p1p3**2*p2p4 - 
     $      8*mt**2*p2p3*p2p4 + (16*p1p3**2*p2p3*p2p4)/s + 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (16*p1p3**2*p2p4**2)/s)
       mat(2,23)=ct1*(-8*mt**4*s + 16*mt**2*p1p2*p1p3 + 
     $      32*mt**2*p1p3*p1p4 - (64*p1p2*p1p3**2*p1p4)/s - 
     $      16*mt**4*p2p3 - 32*mt**2*p1p3*p2p3 + 
     $      (96*mt**2*p1p3*p1p4*p2p3)/s + 
     $      (128*p1p3**2*p1p4*p2p3)/s - 
     $      (128*p1p3**2*p1p4**2*p2p3)/s**2 - 
     $      (32*mt**2*p1p3**2*p2p4)/s + 
     $      (128*p1p3**3*p1p4*p2p4)/s**2) + 
     $      ct3*(-8*mt**4*s - 16*mt**4*p1p3 + 32*mt**2*p1p3*p1p4 + 
     $      (64*mt**2*p1p3**2*p1p4)/s + 16*mt**2*p1p2*p2p3 + 
     $      16*p1p2*p1p3*p2p3 - 
     $      (64*p1p2*p1p3*p1p4*p2p3)/s - 
     $      (32*mt**2*p1p4*p2p3**2)/s - 
     $      (32*p1p3*p1p4*p2p3**2)/s + 
     $      (128*p1p3*p1p4**2*p2p3**2)/s**2 + 
     $      (32*mt**2*p1p3*p2p3*p2p4)/s - 
     $      (32*p1p3**2*p2p3*p2p4)/s - 
     $      (128*p1p3**2*p1p4*p2p3*p2p4)/s**2) + 
     $      ct2*(8*mt**2*s*p1p2 + 16*mt**2*p1p2*p1p3 - 
     $      16*p1p2**2*p2p3 + 16*p1p2*p1p3*p2p3 - 
     $      16*mt**2*p1p4*p2p3 - 
     $      (32*mt**2*p1p3*p1p4*p2p3)/s + 
     $      (64*p1p2*p1p4*p2p3**2)/s - 
     $      (32*p1p3*p1p4*p2p3**2)/s - 
     $      (64*p1p4**2*p2p3**3)/s**2 - 16*mt**2*p1p3*p2p4 - 
     $      (32*mt**2*p1p3**2*p2p4)/s - 16*mt**2*p2p3*p2p4 - 
     $      (32*p1p3**2*p2p3*p2p4)/s + 
     $      (64*p1p3*p1p4*p2p3*p2p4)/s + 
     $      (64*p1p3**2*p2p3*p2p4**2)/s**2)
       mat(3,23)=cu1*(-8*mt**4*s + 16*mt**2*p1p2*p1p3 + 
     $      32*mt**2*p1p3*p1p4 + 16*p1p2*p1p3*p1p4 - 
     $      (64*p1p2*p1p3**2*p1p4)/s - 16*mt**4*p2p3 - 
     $      16*mt**2*p1p3*p2p3 + 
     $      (96*mt**2*p1p3*p1p4*p2p3)/s + 
     $      (64*p1p3**2*p1p4*p2p3)/s - 
     $      (32*p1p3*p1p4**2*p2p3)/s - 
     $      (128*p1p3**2*p1p4**2*p2p3)/s**2 - 
     $      (32*mt**2*p1p3**2*p2p4)/s - 
     $      (32*p1p3**2*p1p4*p2p4)/s + 
     $      (128*p1p3**3*p1p4*p2p4)/s**2) + 
     $      cu3*(8*mt**2*s*p1p2 - 16*p1p2**2*p1p3 + 
     $      16*mt**2*p1p2*p2p3 - 16*mt**2*p1p3*p2p3 - 
     $      16*mt**2*p1p4*p2p3 + (64*p1p3**2*p1p4*p2p3)/s - 
     $      (32*mt**2*p1p4*p2p3**2)/s + 
     $      (64*p1p3*p1p4**2*p2p3**2)/s**2 - 
     $      16*mt**2*p1p3*p2p4 + (64*p1p2*p1p3**2*p2p4)/s - 
     $      (32*mt**2*p1p3*p2p3*p2p4)/s - 
     $      (64*p1p3**3*p2p4**2)/s**2) + 
     $      cu2*(8*mt**2*s*p1p2 + 16*mt**2*p1p2*p1p3 - 
     $      16*p1p2**2*p2p3 + 32*p1p2*p1p3*p2p3 - 
     $      16*mt**2*p1p4*p2p3 - 
     $      (32*mt**2*p1p3*p1p4*p2p3)/s + 
     $      (64*p1p2*p1p4*p2p3**2)/s - 
     $      (64*p1p3*p1p4*p2p3**2)/s - 
     $      (64*p1p4**2*p2p3**3)/s**2 - 16*mt**2*p1p3*p2p4 - 
     $      (32*mt**2*p1p3**2*p2p4)/s - 
     $      (64*p1p3**2*p2p3*p2p4)/s + 
     $      (64*p1p3**2*p2p3*p2p4**2)/s**2)

       do i=1,3
          mat(i,24)=0d0
       end do

       mat(1,25)=cs1*(4*mt**2*s*p1p3 - 8*mt**2*p1p3**2 + 
     $      8*mt**2*p1p3*p1p4 + 4*s*p1p2*p2p3 + 
     $      8*p1p2*p1p3*p2p3 - 8*p1p2*p1p4*p2p3 - 
     $      8*p1p4*p2p3**2 - (16*p1p3*p1p4*p2p3**2)/s + 
     $      (16*p1p4**2*p2p3**2)/s - 24*p1p3*p2p3*p2p4 + 
     $      (16*p1p3**2*p2p3*p2p4)/s - 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s) + 
     $      cs2*(4*mt**2*s*p1p3 + 4*s*p1p2*p2p3 + 
     $      8*mt**2*p1p3*p2p3 - 8*p1p2*p2p3**2 - 
     $      8*p1p4*p2p3**2 + (16*p1p4*p2p3**3)/s - 
     $      8*mt**2*p1p3*p2p4 + 8*p1p2*p2p3*p2p4 - 
     $      24*p1p3*p2p3*p2p4 - 
     $      (16*p1p3*p2p3**2*p2p4)/s - 
     $      (16*p1p4*p2p3**2*p2p4)/s + 
     $      (16*p1p3*p2p3*p2p4**2)/s)
       mat(2,25)=ct1*(4*mt**4*s + 4*mt**2*s*p1p2 + 8*mt**2*s*p1p3 - 
     $      16*p1p2**2*p1p3 - 16*mt**2*p1p3**2 - 
     $      16*mt**2*p1p3*p1p4 + 16*mt**2*p1p2*p2p3 + 
     $      32*p1p2*p1p3*p2p3 - 8*mt**2*p1p4*p2p3 - 
     $      (32*mt**2*p1p4*p2p3**2)/s - 
     $      (64*p1p3*p1p4*p2p3**2)/s + 
     $      (64*p1p3*p1p4**2*p2p3**2)/s**2 - 
     $      8*mt**2*p1p3*p2p4 + (64*p1p2*p1p3**2*p2p4)/s - 
     $      16*mt**2*p2p3*p2p4 - 32*p1p3*p2p3*p2p4 - 
     $      (32*mt**2*p1p3*p2p3*p2p4)/s + 
     $      (64*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (64*p1p3**3*p2p4**2)/s**2) + 
     $      ct3*(4*mt**4*s + 4*mt**2*s*p1p2 + 8*mt**2*s*p1p3 + 
     $      16*mt**2*p1p2*p1p3 - 16*mt**2*p1p3*p1p4 - 
     $      16*mt**4*p2p3 - 32*mt**2*p1p3*p2p3 - 
     $      8*mt**2*p1p4*p2p3 + 
     $      (32*mt**2*p1p3*p1p4*p2p3)/s - 
     $      8*mt**2*p1p3*p2p4 - (32*mt**2*p1p3**2*p2p4)/s - 
     $      16*mt**2*p2p3*p2p4 - 32*p1p3*p2p3*p2p4 - 
     $      (64*p1p2*p1p3*p2p3*p2p4)/s + 
     $      (64*p1p3*p1p4*p2p3*p2p4)/s + 
     $      (64*mt**2*p2p3**2*p2p4)/s + 
     $      (128*p1p3*p2p3**2*p2p4)/s - 
     $      (128*p1p3*p1p4*p2p3**2*p2p4)/s**2 + 
     $      (128*p1p3**2*p2p3*p2p4**2)/s**2) + 
     $      ct2*(-4*mt**4*s - 4*mt**2*s*p1p2 - 16*mt**4*p1p3 + 
     $      8*mt**2*s*p1p3 + 16*mt**2*p1p2*p2p3 - 
     $      32*mt**2*p1p3*p2p3 + 8*mt**2*p1p4*p2p3 - 
     $      (32*mt**2*p1p4*p2p3**2)/s - 8*mt**2*p1p3*p2p4 + 
     $      16*mt**2*p2p3*p2p4 + 16*p1p2*p2p3*p2p4 - 
     $      32*p1p3*p2p3*p2p4 + 
     $      (96*mt**2*p1p3*p2p3*p2p4)/s - 
     $      (64*p1p2*p2p3**2*p2p4)/s + 
     $      (128*p1p3*p2p3**2*p2p4)/s - 
     $      (32*p1p4*p2p3**2*p2p4)/s + 
     $      (128*p1p4*p2p3**3*p2p4)/s**2 + 
     $      (32*p1p3*p2p3*p2p4**2)/s - 
     $      (128*p1p3*p2p3**2*p2p4**2)/s**2)
       mat(3,25)=cu1*(4*mt**4*s + 4*mt**2*s*p1p2 - 16*p1p2**2*p1p3 - 
     $      16*mt**2*p1p3*p1p4 + 16*mt**2*p1p2*p2p3 - 
     $      8*s*p1p2*p2p3 + 32*p1p2*p1p3*p2p3 - 
     $      8*mt**2*p1p4*p2p3 + 16*p1p2*p1p4*p2p3 + 
     $      16*p1p4*p2p3**2 - (32*mt**2*p1p4*p2p3**2)/s - 
     $      (64*p1p3*p1p4*p2p3**2)/s - 
     $      (32*p1p4**2*p2p3**2)/s + 
     $      (64*p1p3*p1p4**2*p2p3**2)/s**2 - 
     $      8*mt**2*p1p3*p2p4 + (64*p1p2*p1p3**2*p2p4)/s - 
     $      16*mt**2*p2p3*p2p4 + 16*p1p3*p2p3*p2p4 - 
     $      (32*mt**2*p1p3*p2p3*p2p4)/s - 
     $      (64*p1p3**2*p2p3*p2p4)/s + 
     $      (32*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (64*p1p3**3*p2p4**2)/s**2) + 
     $      cu3*(-4*mt**4*s - 4*mt**2*s*p1p2 + 16*mt**2*p1p2*p1p3 - 
     $      8*s*p1p2*p2p3 - 16*p1p2**2*p2p3 + 
     $      32*p1p2*p1p3*p2p3 + 8*mt**2*p1p4*p2p3 - 
     $      (32*mt**2*p1p3*p1p4*p2p3)/s + 16*p1p4*p2p3**2 + 
     $      (64*p1p2*p1p4*p2p3**2)/s - 
     $      (64*p1p3*p1p4*p2p3**2)/s - 
     $      (64*p1p4**2*p2p3**3)/s**2 + 8*mt**2*p1p3*p2p4 - 
     $      (32*mt**2*p1p3**2*p2p4)/s + 16*mt**2*p2p3*p2p4 + 
     $      16*p1p2*p2p3*p2p4 + 16*p1p3*p2p3*p2p4 - 
     $      (64*p1p3**2*p2p3*p2p4)/s - 
     $      (32*p1p4*p2p3**2*p2p4)/s - 
     $      (32*p1p3*p2p3*p2p4**2)/s + 
     $      (64*p1p3**2*p2p3*p2p4**2)/s**2) + 
     $      cu2*(-4*mt**4*s - 4*mt**2*s*p1p2 - 16*mt**4*p1p3 + 
     $      16*mt**2*p1p2*p2p3 - 8*s*p1p2*p2p3 - 
     $      32*mt**2*p1p3*p2p3 + 8*mt**2*p1p4*p2p3 + 
     $      16*p1p2*p2p3**2 + 16*p1p4*p2p3**2 - 
     $      (32*mt**2*p1p4*p2p3**2)/s - (32*p1p4*p2p3**3)/s + 
     $      8*mt**2*p1p3*p2p4 + 16*mt**2*p2p3*p2p4 + 
     $      16*p1p2*p2p3*p2p4 + 16*p1p3*p2p3*p2p4 + 
     $      (96*mt**2*p1p3*p2p3*p2p4)/s - 
     $      (64*p1p2*p2p3**2*p2p4)/s + 
     $      (96*p1p3*p2p3**2*p2p4)/s - 
     $      (32*p1p4*p2p3**2*p2p4)/s + 
     $      (128*p1p4*p2p3**3*p2p4)/s**2 - 
     $      (32*p1p3*p2p3*p2p4**2)/s - 
     $      (128*p1p3*p2p3**2*p2p4**2)/s**2)

       mat(1,26)=cs1*(4*s*p1p2*p1p3 - 8*p1p2*p1p3**2 + 
     $      8*p1p2*p1p3*p1p4 + 4*mt**2*s*p2p3 + 
     $      8*mt**2*p1p3*p2p3 - 8*mt**2*p1p4*p2p3 - 
     $      24*p1p3*p1p4*p2p3 - 
     $      (16*p1p3**2*p1p4*p2p3)/s + 
     $      (16*p1p3*p1p4**2*p2p3)/s - 8*p1p3**2*p2p4 + 
     $      (16*p1p3**3*p2p4)/s - (16*p1p3**2*p1p4*p2p4)/s) + 
     $      cs2*(4*s*p1p2*p1p3 + 4*mt**2*s*p2p3 + 
     $      8*p1p2*p1p3*p2p3 - 24*p1p3*p1p4*p2p3 - 
     $      8*mt**2*p2p3**2 + (16*p1p3*p1p4*p2p3**2)/s - 
     $      8*p1p2*p1p3*p2p4 - 8*p1p3**2*p2p4 + 
     $      8*mt**2*p2p3*p2p4 - (16*p1p3**2*p2p3*p2p4)/s - 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s + 
     $      (16*p1p3**2*p2p4**2)/s)
       mat(2,26)=ct1*(4*mt**4*s + 4*mt**2*s*p1p2 - 16*mt**2*p1p2*p1p3+ 
     $      8*s*p1p2*p1p3 - 16*p1p2*p1p3**2 - 
     $      16*mt**2*p1p3*p1p4 - 16*p1p2*p1p3*p1p4 + 
     $      (64*p1p2*p1p3**2*p1p4)/s + 16*mt**4*p2p3 + 
     $      32*mt**2*p1p3*p2p3 - 8*mt**2*p1p4*p2p3 - 
     $      16*p1p3*p1p4*p2p3 - 
     $      (96*mt**2*p1p3*p1p4*p2p3)/s - 
     $      (96*p1p3**2*p1p4*p2p3)/s + 
     $      (32*p1p3*p1p4**2*p2p3)/s + 
     $      (128*p1p3**2*p1p4**2*p2p3)/s**2 - 
     $      8*mt**2*p1p3*p2p4 - 16*p1p3**2*p2p4 + 
     $      (32*mt**2*p1p3**2*p2p4)/s + (32*p1p3**3*p2p4)/s + 
     $      (32*p1p3**2*p1p4*p2p4)/s - 
     $      (128*p1p3**3*p1p4*p2p4)/s**2) + 
     $      ct3*(4*mt**4*s + 4*mt**2*s*p1p2 + 8*s*p1p2*p1p3 + 
     $      16*p1p2**2*p1p3 - 16*mt**2*p1p3*p1p4 - 
     $      16*p1p2*p1p3*p1p4 - 16*mt**2*p1p2*p2p3 - 
     $      32*p1p2*p1p3*p2p3 - 8*mt**2*p1p4*p2p3 - 
     $      16*p1p3*p1p4*p2p3 + 
     $      (32*p1p3*p1p4**2*p2p3)/s + 
     $      (32*mt**2*p1p4*p2p3**2)/s + 
     $      (64*p1p3*p1p4*p2p3**2)/s - 
     $      (64*p1p3*p1p4**2*p2p3**2)/s**2 - 
     $      8*mt**2*p1p3*p2p4 - 16*p1p3**2*p2p4 - 
     $      (64*p1p2*p1p3**2*p2p4)/s + 
     $      (32*p1p3**2*p1p4*p2p4)/s + 
     $      (32*mt**2*p1p3*p2p3*p2p4)/s + 
     $      (64*p1p3**2*p2p3*p2p4)/s + 
     $      (64*p1p3**3*p2p4**2)/s**2) + 
     $      ct2*(-4*mt**4*s - 4*mt**2*s*p1p2 - 16*mt**2*p1p2*p1p3 + 
     $      8*s*p1p2*p1p3 + 16*mt**2*p1p3*p1p4 + 
     $      16*p1p2**2*p2p3 - 32*p1p2*p1p3*p2p3 + 
     $      8*mt**2*p1p4*p2p3 - 16*p1p3*p1p4*p2p3 + 
     $      (32*mt**2*p1p3*p1p4*p2p3)/s - 
     $      (64*p1p2*p1p4*p2p3**2)/s + 
     $      (64*p1p3*p1p4*p2p3**2)/s + 
     $      (64*p1p4**2*p2p3**3)/s**2 + 8*mt**2*p1p3*p2p4 - 
     $      16*p1p2*p1p3*p2p4 - 16*p1p3**2*p2p4 + 
     $      (32*mt**2*p1p3**2*p2p4)/s + 16*mt**2*p2p3*p2p4 + 
     $      (64*p1p3**2*p2p3*p2p4)/s - 
     $      (32*p1p3*p1p4*p2p3*p2p4)/s + 
     $      (32*p1p3**2*p2p4**2)/s - 
     $      (64*p1p3**2*p2p3*p2p4**2)/s**2)
       mat(3,26)=cu1*(4*mt**4*s + 4*mt**2*s*p1p2 - 16*mt**2*p1p2*p1p3- 
     $      16*mt**2*p1p3*p1p4 - 16*p1p2*p1p3*p1p4 + 
     $      (64*p1p2*p1p3**2*p1p4)/s + 16*mt**4*p2p3 - 
     $      8*mt**2*s*p2p3 + 32*mt**2*p1p3*p2p3 + 
     $      8*mt**2*p1p4*p2p3 + 32*p1p3*p1p4*p2p3 - 
     $      (96*mt**2*p1p3*p1p4*p2p3)/s - 
     $      (128*p1p3**2*p1p4*p2p3)/s - 
     $      (32*p1p3*p1p4**2*p2p3)/s + 
     $      (128*p1p3**2*p1p4**2*p2p3)/s**2 - 
     $      8*mt**2*p1p3*p2p4 + (32*mt**2*p1p3**2*p2p4)/s + 
     $      (32*p1p3**2*p1p4*p2p4)/s - 
     $      (128*p1p3**3*p1p4*p2p4)/s**2) + 
     $      cu3*(-4*mt**4*s - 4*mt**2*s*p1p2 + 16*mt**4*p1p3 + 
     $      16*mt**2*p1p3*p1p4 - (64*mt**2*p1p3**2*p1p4)/s - 
     $      8*mt**2*s*p2p3 - 16*mt**2*p1p2*p2p3 + 
     $      32*mt**2*p1p3*p2p3 + 8*mt**2*p1p4*p2p3 + 
     $      32*p1p3*p1p4*p2p3 + 
     $      (64*p1p2*p1p3*p1p4*p2p3)/s - 
     $      (128*p1p3**2*p1p4*p2p3)/s + 
     $      (32*mt**2*p1p4*p2p3**2)/s - 
     $      (128*p1p3*p1p4**2*p2p3**2)/s**2 + 
     $      8*mt**2*p1p3*p2p4 + 16*mt**2*p2p3*p2p4 - 
     $      (32*mt**2*p1p3*p2p3*p2p4)/s - 
     $      (64*p1p3*p1p4*p2p3*p2p4)/s + 
     $      (128*p1p3**2*p1p4*p2p3*p2p4)/s**2) + 
     $      cu2*(-4*mt**4*s - 4*mt**2*s*p1p2 - 16*mt**2*p1p2*p1p3 + 
     $      16*mt**2*p1p3*p1p4 - 8*mt**2*s*p2p3 + 
     $      16*p1p2**2*p2p3 - 32*p1p2*p1p3*p2p3 + 
     $      8*mt**2*p1p4*p2p3 + 32*p1p3*p1p4*p2p3 + 
     $      (32*mt**2*p1p3*p1p4*p2p3)/s + 16*mt**2*p2p3**2 - 
     $      (64*p1p2*p1p4*p2p3**2)/s + 
     $      (64*p1p4**2*p2p3**3)/s**2 + 8*mt**2*p1p3*p2p4 + 
     $      (32*mt**2*p1p3**2*p2p4)/s + 16*mt**2*p2p3*p2p4 + 
     $      (64*p1p3**2*p2p3*p2p4)/s - 
     $      (64*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (64*p1p3**2*p2p3*p2p4**2)/s**2)

       do i=1,3
          mat(i,27)=0d0
       end do

       mat(1,28)=cs1*(-8*mt**3*p1p3+8*mt*p1p2*p1p3+
     $      8*mt**3*p1p4-8*mt*p1p2*p1p4-
     $      (16*mt*p1p3*p1p4*p2p3)/s+
     $      (16*mt*p1p4**2*p2p3)/s-(16*mt*p1p3**2*p2p4)/s+
     $      (16*mt*p1p3*p1p4*p2p4)/s+
     $      (32*mt*p1p3*p2p3*p2p4)/s-
     $      (32*mt*p1p4*p2p3*p2p4)/s)+
     $      cs2*(8*mt**3*p2p3-8*mt*p1p2*p2p3+
     $      (16*mt*p1p4*p2p3**2)/s-8*mt**3*p2p4+
     $      8*mt*p1p2*p2p4+(16*mt*p1p3*p2p3*p2p4)/s-
     $      (16*mt*p1p4*p2p3*p2p4)/s-
     $      (32*mt*p2p3**2*p2p4)/s-(16*mt*p1p3*p2p4**2)/s+
     $      (32*mt*p2p3*p2p4**2)/s)

       mat(2,28)=ct1*(16*mt**3*p1p2-4*mt*s*p1p2-16*mt*p1p2**2+
     $      24*mt*p1p2*p1p3+16*mt*p1p2*p1p4-
     $      (64*mt*p1p2*p1p3*p1p4)/s+8*mt*p1p2*p2p3+
     $      8*mt*p1p4*p2p3-(32*mt**3*p1p4*p2p3)/s+
     $      (64*mt*p1p2*p1p4*p2p3)/s-
     $      (48*mt*p1p3*p1p4*p2p3)/s-
     $      (32*mt*p1p4**2*p2p3)/s+
     $      (128*mt*p1p3*p1p4**2*p2p3)/s**2-
     $      (16*mt*p1p4*p2p3**2)/s-
     $      (64*mt*p1p4**2*p2p3**2)/s**2+16*mt*p1p2*p2p4+
     $      8*mt*p1p3*p2p4-(32*mt**3*p1p3*p2p4)/s+
     $      (64*mt*p1p2*p1p3*p2p4)/s-
     $      (48*mt*p1p3**2*p2p4)/s-
     $      (32*mt*p1p3*p1p4*p2p4)/s+
     $      (128*mt*p1p3**2*p1p4*p2p4)/s**2-
     $      (16*mt*p1p3*p2p3*p2p4)/s-
     $      (32*mt*p1p4*p2p3*p2p4)/s-
     $      (128*mt*p1p3*p1p4*p2p3*p2p4)/s**2-
     $      (32*mt*p1p3*p2p4**2)/s-(64*mt*p1p3**2*p2p4**2)/s**2)+
     $      ct3*(16*mt**3*p1p2-4*mt*s*p1p2-16*mt*p1p2**2-
     $      8*mt**3*p1p3+16*mt*p1p2*p1p4-8*mt**3*p2p3+
     $      8*mt*p1p4*p2p3-(32*mt**3*p1p4*p2p3)/s+
     $      (64*mt*p1p2*p1p4*p2p3)/s-
     $      (32*mt*p1p4**2*p2p3)/s-
     $      (64*mt*p1p4**2*p2p3**2)/s**2+16*mt*p1p2*p2p4+
     $      8*mt*p1p3*p2p4-(32*mt**3*p1p3*p2p4)/s+
     $      (64*mt*p1p2*p1p3*p2p4)/s-
     $      (32*mt*p1p3*p1p4*p2p4)/s-
     $      (64*mt*p1p2*p2p3*p2p4)/s+
     $      (32*mt*p1p3*p2p3*p2p4)/s-
     $      (32*mt*p1p4*p2p3*p2p4)/s-
     $      (128*mt*p1p3*p1p4*p2p3*p2p4)/s**2+
     $      (32*mt*p2p3**2*p2p4)/s+
     $      (128*mt*p1p4*p2p3**2*p2p4)/s**2-
     $      (32*mt*p1p3*p2p4**2)/s-
     $      (64*mt*p1p3**2*p2p4**2)/s**2+
     $      (128*mt*p1p3*p2p3*p2p4**2)/s**2)+
     $      ct2*(-16*mt**5+16*mt**3*p1p2-4*mt*s*p1p2-
     $      8*mt**3*p1p3-16*mt**3*p1p4-8*mt**3*p2p3+
     $      8*mt*p1p4*p2p3-(32*mt**3*p1p4*p2p3)/s-
     $      16*mt**3*p2p4+16*mt*p1p2*p2p4+
     $      8*mt*p1p3*p2p4-(32*mt**3*p1p3*p2p4)/s+
     $      (128*mt**3*p2p3*p2p4)/s-
     $      (64*mt*p1p2*p2p3*p2p4)/s+
     $      (32*mt*p1p3*p2p3*p2p4)/s+
     $      (32*mt*p1p4*p2p3*p2p4)/s+
     $      (32*mt*p2p3**2*p2p4)/s+
     $      (128*mt*p1p4*p2p3**2*p2p4)/s**2-
     $      (32*mt*p1p3*p2p4**2)/s+(64*mt*p2p3*p2p4**2)/s+
     $      (128*mt*p1p3*p2p3*p2p4**2)/s**2-
     $      (256*mt*p2p3**2*p2p4**2)/s**2)

       mat(3,28)=cu1*(4*mt**3*s+16*mt**3*p1p2-16*mt*p1p2**2+
     $      8*mt*p1p2*p1p3-16*mt**3*p1p4+
     $      16*mt*p1p2*p1p4-(64*mt*p1p2*p1p3*p1p4)/s+
     $      8*mt*p1p2*p2p3-(32*mt**3*p1p4*p2p3)/s+
     $      (64*mt*p1p2*p1p4*p2p3)/s-
     $      (16*mt*p1p3*p1p4*p2p3)/s-
     $      (32*mt*p1p4**2*p2p3)/s+
     $      (128*mt*p1p3*p1p4**2*p2p3)/s**2-
     $      (16*mt*p1p4*p2p3**2)/s-
     $      (64*mt*p1p4**2*p2p3**2)/s**2+16*mt*p1p2*p2p4-
     $      (32*mt**3*p1p3*p2p4)/s+
     $      (64*mt*p1p2*p1p3*p2p4)/s-
     $      (16*mt*p1p3**2*p2p4)/s-
     $      (32*mt*p1p3*p1p4*p2p4)/s+
     $      (128*mt*p1p3**2*p1p4*p2p4)/s**2-
     $      16*mt*p2p3*p2p4-(16*mt*p1p3*p2p3*p2p4)/s+
     $      (32*mt*p1p4*p2p3*p2p4)/s-
     $      (128*mt*p1p3*p1p4*p2p3*p2p4)/s**2-
     $      (32*mt*p1p3*p2p4**2)/s-(64*mt*p1p3**2*p2p4**2)/s**2)+
     $      cu3*(-16*mt**5+4*mt**3*s+16*mt**3*p1p2+
     $      8*mt*p1p2*p1p3-16*mt**3*p1p4+
     $      (64*mt**3*p1p3*p1p4)/s+8*mt*p1p2*p2p3-
     $      (32*mt**3*p1p4*p2p3)/s-
     $      (16*mt*p1p3*p1p4*p2p3)/s-
     $      (16*mt*p1p4*p2p3**2)/s-16*mt**3*p2p4-
     $      (32*mt**3*p1p3*p2p4)/s-(16*mt*p1p3**2*p2p4)/s-
     $      16*mt*p2p3*p2p4+(64*mt**3*p2p3*p2p4)/s-
     $      (64*mt*p1p2*p2p3*p2p4)/s-
     $      (16*mt*p1p3*p2p3*p2p4)/s+
     $      (64*mt*p1p4*p2p3*p2p4)/s-
     $      (256*mt*p1p3*p1p4*p2p3*p2p4)/s**2+
     $      (128*mt*p1p4*p2p3**2*p2p4)/s**2+
     $      (64*mt*p2p3*p2p4**2)/s+
     $      (128*mt*p1p3*p2p3*p2p4**2)/s**2)+
     $      cu2*(-16*mt**5+4*mt**3*s+16*mt**3*p1p2-8*mt**3*p1p3-
     $      16*mt**3*p1p4-24*mt**3*p2p3-
     $      (32*mt**3*p1p4*p2p3)/s-16*mt**3*p2p4-
     $      (32*mt**3*p1p3*p2p4)/s-16*mt*p2p3*p2p4+
     $      (128*mt**3*p2p3*p2p4)/s-
     $      (64*mt*p1p2*p2p3*p2p4)/s+
     $      (32*mt*p1p3*p2p3*p2p4)/s+
     $      (64*mt*p1p4*p2p3*p2p4)/s+
     $      (96*mt*p2p3**2*p2p4)/s+
     $      (128*mt*p1p4*p2p3**2*p2p4)/s**2+
     $      (64*mt*p2p3*p2p4**2)/s+
     $      (128*mt*p1p3*p2p3*p2p4**2)/s**2-
     $      (256*mt*p2p3**2*p2p4**2)/s**2)

       mat(1,29)=cs1*(8*mt**3*p1p3 - 8*mt*p1p2*p1p3 - 
     $      8*mt**3*p1p4 + 8*mt*p1p2*p1p4 - 
     $      (32*mt*p1p3**2*p1p4)/s + (32*mt*p1p3*p1p4**2)/s + 
     $      (16*mt*p1p3*p1p4*p2p3)/s - 
     $      (16*mt*p1p4**2*p2p3)/s + (16*mt*p1p3**2*p2p4)/s - 
     $      (16*mt*p1p3*p1p4*p2p4)/s) + 
     $      cs2*(-8*mt**3*p2p3 + 8*mt*p1p2*p2p3 + 
     $      (32*mt*p1p3*p1p4*p2p3)/s - 
     $      (16*mt*p1p4*p2p3**2)/s + 8*mt**3*p2p4 - 
     $      8*mt*p1p2*p2p4 - (32*mt*p1p3*p1p4*p2p4)/s - 
     $      (16*mt*p1p3*p2p3*p2p4)/s + 
     $      (16*mt*p1p4*p2p3*p2p4)/s + (16*mt*p1p3*p2p4**2)/s)

       mat(2,29)=ct1*(16*mt**5 - 4*mt**3*s - 16*mt**3*p1p2 + 
     $      24*mt**3*p1p3 + 16*mt**3*p1p4 + 16*mt*p1p3*p1p4 - 
     $      (128*mt**3*p1p3*p1p4)/s + 
     $      (64*mt*p1p2*p1p3*p1p4)/s - 
     $      (96*mt*p1p3**2*p1p4)/s - (64*mt*p1p3*p1p4**2)/s + 
     $      (256*mt*p1p3**2*p1p4**2)/s**2 + 8*mt**3*p2p3 + 
     $      (32*mt**3*p1p4*p2p3)/s - 
     $      (32*mt*p1p3*p1p4*p2p3)/s - 
     $      (128*mt*p1p3*p1p4**2*p2p3)/s**2 + 16*mt**3*p2p4 + 
     $      (32*mt**3*p1p3*p2p4)/s - 
     $      (64*mt*p1p3*p1p4*p2p4)/s - 
     $      (128*mt*p1p3**2*p1p4*p2p4)/s**2) + 
     $      ct3*(16*mt**5 - 4*mt**3*s - 16*mt**3*p1p2 - 
     $      8*mt*p1p2*p1p3 + 16*mt**3*p1p4 + 
     $      16*mt*p1p3*p1p4 - (64*mt**3*p1p3*p1p4)/s + 
     $      (64*mt*p1p2*p1p3*p1p4)/s - 
     $      (64*mt*p1p3*p1p4**2)/s - 8*mt*p1p2*p2p3 + 
     $      (32*mt**3*p1p4*p2p3)/s + 
     $      (16*mt*p1p3*p1p4*p2p3)/s - 
     $      (128*mt*p1p3*p1p4**2*p2p3)/s**2 + 
     $      (16*mt*p1p4*p2p3**2)/s + 16*mt**3*p2p4 + 
     $      (32*mt**3*p1p3*p2p4)/s + (16*mt*p1p3**2*p2p4)/s - 
     $      (64*mt*p1p3*p1p4*p2p4)/s - 
     $      (128*mt*p1p3**2*p1p4*p2p4)/s**2 - 
     $      (64*mt**3*p2p3*p2p4)/s + 
     $      (16*mt*p1p3*p2p3*p2p4)/s + 
     $      (256*mt*p1p3*p1p4*p2p3*p2p4)/s**2) + 
     $      ct2*(-4*mt**3*s - 16*mt**3*p1p2 + 16*mt*p1p2**2 - 
     $      8*mt*p1p2*p1p3 - 16*mt*p1p2*p1p4 + 
     $      16*mt*p1p3*p1p4 - 8*mt*p1p2*p2p3 + 
     $      (32*mt**3*p1p4*p2p3)/s - 
     $      (64*mt*p1p2*p1p4*p2p3)/s + 
     $      (16*mt*p1p3*p1p4*p2p3)/s + 
     $      (32*mt*p1p4**2*p2p3)/s + (16*mt*p1p4*p2p3**2)/s + 
     $      (64*mt*p1p4**2*p2p3**2)/s**2 + 16*mt**3*p2p4 - 
     $      16*mt*p1p2*p2p4 + (32*mt**3*p1p3*p2p4)/s - 
     $      (64*mt*p1p2*p1p3*p2p4)/s + 
     $      (16*mt*p1p3**2*p2p4)/s - 
     $      (32*mt*p1p3*p1p4*p2p4)/s + 
     $      (64*mt*p1p2*p2p3*p2p4)/s + 
     $      (16*mt*p1p3*p2p3*p2p4)/s + 
     $      (32*mt*p1p4*p2p3*p2p4)/s + 
     $      (128*mt*p1p3*p1p4*p2p3*p2p4)/s**2 - 
     $      (128*mt*p1p4*p2p3**2*p2p4)/s**2 + 
     $      (32*mt*p1p3*p2p4**2)/s + 
     $      (64*mt*p1p3**2*p2p4**2)/s**2 - 
     $      (128*mt*p1p3*p2p3*p2p4**2)/s**2)

       mat(3,29)=cu1*(16*mt**5 - 16*mt**3*p1p2 + 4*mt*s*p1p2 + 
     $      8*mt**3*p1p3 + 16*mt**3*p1p4 - 16*mt*p1p2*p1p4 - 
     $      (128*mt**3*p1p3*p1p4)/s + 
     $      (64*mt*p1p2*p1p3*p1p4)/s - 
     $      (32*mt*p1p3**2*p1p4)/s - (64*mt*p1p3*p1p4**2)/s + 
     $      (256*mt*p1p3**2*p1p4**2)/s**2 + 8*mt**3*p2p3 - 
     $      8*mt*p1p4*p2p3 + (32*mt**3*p1p4*p2p3)/s - 
     $      (32*mt*p1p3*p1p4*p2p3)/s + 
     $      (32*mt*p1p4**2*p2p3)/s - 
     $      (128*mt*p1p3*p1p4**2*p2p3)/s**2 + 16*mt**3*p2p4 - 
     $      8*mt*p1p3*p2p4 + (32*mt**3*p1p3*p2p4)/s - 
     $      (32*mt*p1p3*p1p4*p2p4)/s - 
     $      (128*mt*p1p3**2*p1p4*p2p4)/s**2) + 
     $      cu3*(-16*mt**3*p1p2 + 4*mt*s*p1p2 + 16*mt*p1p2**2 + 
     $      8*mt**3*p1p3 - 16*mt*p1p2*p1p4 + 
     $      (64*mt*p1p2*p1p3*p1p4)/s - 
     $      (32*mt*p1p3**2*p1p4)/s + 8*mt**3*p2p3 - 
     $      8*mt*p1p4*p2p3 + (32*mt**3*p1p4*p2p3)/s - 
     $      (64*mt*p1p2*p1p4*p2p3)/s - 
     $      (32*mt*p1p3*p1p4*p2p3)/s + 
     $      (32*mt*p1p4**2*p2p3)/s - 
     $      (128*mt*p1p3*p1p4**2*p2p3)/s**2 + 
     $      (64*mt*p1p4**2*p2p3**2)/s**2 - 16*mt*p1p2*p2p4 - 
     $      8*mt*p1p3*p2p4 + (32*mt**3*p1p3*p2p4)/s - 
     $      (64*mt*p1p2*p1p3*p2p4)/s + 
     $      (32*mt*p1p3*p1p4*p2p4)/s - 
     $      (128*mt*p1p3**2*p1p4*p2p4)/s**2 + 
     $      (32*mt*p1p4*p2p3*p2p4)/s + 
     $      (128*mt*p1p3*p1p4*p2p3*p2p4)/s**2 + 
     $      (32*mt*p1p3*p2p4**2)/s + (64*mt*p1p3**2*p2p4**2)/s**2) + 
     $      cu2*(-16*mt**3*p1p2 + 4*mt*s*p1p2 + 16*mt*p1p2**2 - 
     $      8*mt*p1p2*p1p3 - 16*mt*p1p2*p1p4 - 
     $      24*mt*p1p2*p2p3 - 8*mt*p1p4*p2p3 + 
     $      (32*mt**3*p1p4*p2p3)/s - 
     $      (64*mt*p1p2*p1p4*p2p3)/s + 
     $      (16*mt*p1p3*p1p4*p2p3)/s + 
     $      (32*mt*p1p4**2*p2p3)/s + (48*mt*p1p4*p2p3**2)/s + 
     $      (64*mt*p1p4**2*p2p3**2)/s**2 - 16*mt*p1p2*p2p4 - 
     $      8*mt*p1p3*p2p4 + (32*mt**3*p1p3*p2p4)/s - 
     $      (64*mt*p1p2*p1p3*p2p4)/s + 
     $      (16*mt*p1p3**2*p2p4)/s + 
     $      (32*mt*p1p3*p1p4*p2p4)/s + 
     $      (64*mt*p1p2*p2p3*p2p4)/s + 
     $      (48*mt*p1p3*p2p3*p2p4)/s + 
     $      (32*mt*p1p4*p2p3*p2p4)/s + 
     $      (128*mt*p1p3*p1p4*p2p3*p2p4)/s**2 - 
     $      (128*mt*p1p4*p2p3**2*p2p4)/s**2 + 
     $      (32*mt*p1p3*p2p4**2)/s + 
     $      (64*mt*p1p3**2*p2p4**2)/s**2 - 
     $      (128*mt*p1p3*p2p3*p2p4**2)/s**2)

       do i=1,3
          mat(i,30)=0d0
       end do

       mat(1,31)=cs1*(-8*mt**3*p1p3+8*mt*p1p2*p1p3+
     $      8*mt**3*p1p4-8*mt*p1p2*p1p4-
     $      (16*mt*p1p3*p1p4*p2p3)/s+
     $      (16*mt*p1p4**2*p2p3)/s-(16*mt*p1p3**2*p2p4)/s+
     $      (16*mt*p1p3*p1p4*p2p4)/s+
     $      (32*mt*p1p3*p2p3*p2p4)/s-
     $      (32*mt*p1p4*p2p3*p2p4)/s)+
     $      cs2*(8*mt**3*p2p3-8*mt*p1p2*p2p3+
     $      (16*mt*p1p4*p2p3**2)/s-8*mt**3*p2p4+
     $      8*mt*p1p2*p2p4+(16*mt*p1p3*p2p3*p2p4)/s-
     $      (16*mt*p1p4*p2p3*p2p4)/s-
     $      (32*mt*p2p3**2*p2p4)/s-(16*mt*p1p3*p2p4**2)/s+
     $      (32*mt*p2p3*p2p4**2)/s)

       mat(2,31)=ct1*(4*mt**3*s+16*mt**3*p1p2-16*mt*p1p2**2-
     $      16*mt**3*p1p3+16*mt*p1p2*p1p3+
     $      8*mt*p1p2*p1p4-(64*mt*p1p2*p1p3*p1p4)/s+
     $      16*mt*p1p2*p2p3-(32*mt**3*p1p4*p2p3)/s+
     $      (64*mt*p1p2*p1p4*p2p3)/s-
     $      (32*mt*p1p3*p1p4*p2p3)/s-
     $      (16*mt*p1p4**2*p2p3)/s+
     $      (128*mt*p1p3*p1p4**2*p2p3)/s**2-
     $      (32*mt*p1p4*p2p3**2)/s-
     $      (64*mt*p1p4**2*p2p3**2)/s**2+8*mt*p1p2*p2p4-
     $      (32*mt**3*p1p3*p2p4)/s+
     $      (64*mt*p1p2*p1p3*p2p4)/s-
     $      (32*mt*p1p3**2*p2p4)/s-
     $      (16*mt*p1p3*p1p4*p2p4)/s+
     $      (128*mt*p1p3**2*p1p4*p2p4)/s**2-
     $      16*mt*p2p3*p2p4+(32*mt*p1p3*p2p3*p2p4)/s-
     $      (16*mt*p1p4*p2p3*p2p4)/s-
     $      (128*mt*p1p3*p1p4*p2p3*p2p4)/s**2-
     $      (16*mt*p1p3*p2p4**2)/s-(64*mt*p1p3**2*p2p4**2)/s**2)+
     $      ct3*(-16*mt**5+4*mt**3*s+16*mt**3*p1p2-16*mt**3*p1p3+
     $      8*mt*p1p2*p1p4+(64*mt**3*p1p3*p1p4)/s-
     $      16*mt**3*p2p3-(32*mt**3*p1p4*p2p3)/s-
     $      (16*mt*p1p4**2*p2p3)/s+8*mt*p1p2*p2p4-
     $      (32*mt**3*p1p3*p2p4)/s-
     $      (16*mt*p1p3*p1p4*p2p4)/s-16*mt*p2p3*p2p4+
     $      (64*mt**3*p2p3*p2p4)/s-
     $      (64*mt*p1p2*p2p3*p2p4)/s+
     $      (64*mt*p1p3*p2p3*p2p4)/s-
     $      (16*mt*p1p4*p2p3*p2p4)/s-
     $      (256*mt*p1p3*p1p4*p2p3*p2p4)/s**2+
     $      (64*mt*p2p3**2*p2p4)/s+
     $      (128*mt*p1p4*p2p3**2*p2p4)/s**2-
     $      (16*mt*p1p3*p2p4**2)/s+
     $      (128*mt*p1p3*p2p3*p2p4**2)/s**2)+
     $      ct2*(-16*mt**5+4*mt**3*s+16*mt**3*p1p2-16*mt**3*p1p3-
     $      8*mt**3*p1p4-16*mt**3*p2p3-
     $      (32*mt**3*p1p4*p2p3)/s-24*mt**3*p2p4-
     $      (32*mt**3*p1p3*p2p4)/s-16*mt*p2p3*p2p4+
     $      (128*mt**3*p2p3*p2p4)/s-
     $      (64*mt*p1p2*p2p3*p2p4)/s+
     $      (64*mt*p1p3*p2p3*p2p4)/s+
     $      (32*mt*p1p4*p2p3*p2p4)/s+
     $      (64*mt*p2p3**2*p2p4)/s+
     $      (128*mt*p1p4*p2p3**2*p2p4)/s**2+
     $      (96*mt*p2p3*p2p4**2)/s+
     $      (128*mt*p1p3*p2p3*p2p4**2)/s**2-
     $      (256*mt*p2p3**2*p2p4**2)/s**2)

       mat(3,31)=cu1*(16*mt**3*p1p2-4*mt*s*p1p2-16*mt*p1p2**2+
     $      16*mt*p1p2*p1p3+24*mt*p1p2*p1p4-
     $      (64*mt*p1p2*p1p3*p1p4)/s+16*mt*p1p2*p2p3+
     $      8*mt*p1p4*p2p3-(32*mt**3*p1p4*p2p3)/s+
     $      (64*mt*p1p2*p1p4*p2p3)/s-
     $      (32*mt*p1p3*p1p4*p2p3)/s-
     $      (48*mt*p1p4**2*p2p3)/s+
     $      (128*mt*p1p3*p1p4**2*p2p3)/s**2-
     $      (32*mt*p1p4*p2p3**2)/s-
     $      (64*mt*p1p4**2*p2p3**2)/s**2+8*mt*p1p2*p2p4+
     $      8*mt*p1p3*p2p4-(32*mt**3*p1p3*p2p4)/s+
     $      (64*mt*p1p2*p1p3*p2p4)/s-
     $      (32*mt*p1p3**2*p2p4)/s-
     $      (48*mt*p1p3*p1p4*p2p4)/s+
     $      (128*mt*p1p3**2*p1p4*p2p4)/s**2-
     $      (32*mt*p1p3*p2p3*p2p4)/s-
     $      (16*mt*p1p4*p2p3*p2p4)/s-
     $      (128*mt*p1p3*p1p4*p2p3*p2p4)/s**2-
     $      (16*mt*p1p3*p2p4**2)/s-(64*mt*p1p3**2*p2p4**2)/s**2)+
     $      cu3*(16*mt**3*p1p2-4*mt*s*p1p2-16*mt*p1p2**2+
     $      16*mt*p1p2*p1p3-8*mt**3*p1p4+
     $      16*mt*p1p2*p2p3+8*mt*p1p4*p2p3-
     $      (32*mt**3*p1p4*p2p3)/s+
     $      (64*mt*p1p2*p1p4*p2p3)/s-
     $      (32*mt*p1p3*p1p4*p2p3)/s-
     $      (32*mt*p1p4*p2p3**2)/s-
     $      (64*mt*p1p4**2*p2p3**2)/s**2-8*mt**3*p2p4+
     $      8*mt*p1p3*p2p4-(32*mt**3*p1p3*p2p4)/s+
     $      (64*mt*p1p2*p1p3*p2p4)/s-
     $      (32*mt*p1p3**2*p2p4)/s-
     $      (64*mt*p1p2*p2p3*p2p4)/s-
     $      (32*mt*p1p3*p2p3*p2p4)/s+
     $      (32*mt*p1p4*p2p3*p2p4)/s-
     $      (128*mt*p1p3*p1p4*p2p3*p2p4)/s**2+
     $      (128*mt*p1p4*p2p3**2*p2p4)/s**2-
     $      (64*mt*p1p3**2*p2p4**2)/s**2+
     $      (32*mt*p2p3*p2p4**2)/s+
     $      (128*mt*p1p3*p2p3*p2p4**2)/s**2)+
     $      cu2*(-16*mt**5+16*mt**3*p1p2-4*mt*s*p1p2-
     $      16*mt**3*p1p3-8*mt**3*p1p4-16*mt**3*p2p3+
     $      16*mt*p1p2*p2p3+8*mt*p1p4*p2p3-
     $      (32*mt**3*p1p4*p2p3)/s-(32*mt*p1p4*p2p3**2)/s-
     $      8*mt**3*p2p4+8*mt*p1p3*p2p4-
     $      (32*mt**3*p1p3*p2p4)/s+(128*mt**3*p2p3*p2p4)/s-
     $      (64*mt*p1p2*p2p3*p2p4)/s+
     $      (32*mt*p1p3*p2p3*p2p4)/s+
     $      (32*mt*p1p4*p2p3*p2p4)/s+
     $      (64*mt*p2p3**2*p2p4)/s+
     $      (128*mt*p1p4*p2p3**2*p2p4)/s**2+
     $      (32*mt*p2p3*p2p4**2)/s+
     $      (128*mt*p1p3*p2p3*p2p4**2)/s**2-
     $      (256*mt*p2p3**2*p2p4**2)/s**2)

       mat(1,32)=cs1*(8*mt**3*p1p3 - 8*mt*p1p2*p1p3 - 
     $      8*mt**3*p1p4 + 8*mt*p1p2*p1p4 - 
     $      (32*mt*p1p3**2*p1p4)/s + (32*mt*p1p3*p1p4**2)/s + 
     $      (16*mt*p1p3*p1p4*p2p3)/s - 
     $      (16*mt*p1p4**2*p2p3)/s + (16*mt*p1p3**2*p2p4)/s - 
     $      (16*mt*p1p3*p1p4*p2p4)/s) + 
     $      cs2*(-8*mt**3*p2p3 + 8*mt*p1p2*p2p3 + 
     $      (32*mt*p1p3*p1p4*p2p3)/s - 
     $      (16*mt*p1p4*p2p3**2)/s + 8*mt**3*p2p4 - 
     $      8*mt*p1p2*p2p4 - (32*mt*p1p3*p1p4*p2p4)/s - 
     $      (16*mt*p1p3*p2p3*p2p4)/s + 
     $      (16*mt*p1p4*p2p3*p2p4)/s + (16*mt*p1p3*p2p4**2)/s)

       mat(2,32)=ct1*(16*mt**5 - 16*mt**3*p1p2 + 4*mt*s*p1p2 + 
     $      16*mt**3*p1p3 - 16*mt*p1p2*p1p3 + 8*mt**3*p1p4 - 
     $      (128*mt**3*p1p3*p1p4)/s + 
     $      (64*mt*p1p2*p1p3*p1p4)/s - 
     $      (64*mt*p1p3**2*p1p4)/s - (32*mt*p1p3*p1p4**2)/s + 
     $      (256*mt*p1p3**2*p1p4**2)/s**2 + 16*mt**3*p2p3 - 
     $      8*mt*p1p4*p2p3 + (32*mt**3*p1p4*p2p3)/s - 
     $      (32*mt*p1p3*p1p4*p2p3)/s - 
     $      (128*mt*p1p3*p1p4**2*p2p3)/s**2 + 8*mt**3*p2p4 - 
     $      8*mt*p1p3*p2p4 + (32*mt**3*p1p3*p2p4)/s + 
     $      (32*mt*p1p3**2*p2p4)/s - 
     $      (32*mt*p1p3*p1p4*p2p4)/s - 
     $      (128*mt*p1p3**2*p1p4*p2p4)/s**2) + 
     $      ct3*(-16*mt**3*p1p2 + 4*mt*s*p1p2 + 16*mt*p1p2**2 - 
     $      16*mt*p1p2*p1p3 + 8*mt**3*p1p4 + 
     $      (64*mt*p1p2*p1p3*p1p4)/s - 
     $      (32*mt*p1p3*p1p4**2)/s - 16*mt*p1p2*p2p3 - 
     $      8*mt*p1p4*p2p3 + (32*mt**3*p1p4*p2p3)/s - 
     $      (64*mt*p1p2*p1p4*p2p3)/s + 
     $      (32*mt*p1p3*p1p4*p2p3)/s - 
     $      (128*mt*p1p3*p1p4**2*p2p3)/s**2 + 
     $      (32*mt*p1p4*p2p3**2)/s + 
     $      (64*mt*p1p4**2*p2p3**2)/s**2 + 8*mt**3*p2p4 - 
     $      8*mt*p1p3*p2p4 + (32*mt**3*p1p3*p2p4)/s - 
     $      (64*mt*p1p2*p1p3*p2p4)/s + 
     $      (32*mt*p1p3**2*p2p4)/s - 
     $      (32*mt*p1p3*p1p4*p2p4)/s - 
     $      (128*mt*p1p3**2*p1p4*p2p4)/s**2 + 
     $      (32*mt*p1p3*p2p3*p2p4)/s + 
     $      (128*mt*p1p3*p1p4*p2p3*p2p4)/s**2 + 
     $      (64*mt*p1p3**2*p2p4**2)/s**2) + 
     $      ct2*(-16*mt**3*p1p2 + 4*mt*s*p1p2 + 16*mt*p1p2**2 - 
     $      16*mt*p1p2*p1p3 - 8*mt*p1p2*p1p4 - 
     $      16*mt*p1p2*p2p3 - 8*mt*p1p4*p2p3 + 
     $      (32*mt**3*p1p4*p2p3)/s - 
     $      (64*mt*p1p2*p1p4*p2p3)/s + 
     $      (32*mt*p1p3*p1p4*p2p3)/s + 
     $      (16*mt*p1p4**2*p2p3)/s + (32*mt*p1p4*p2p3**2)/s + 
     $      (64*mt*p1p4**2*p2p3**2)/s**2 - 24*mt*p1p2*p2p4 - 
     $      8*mt*p1p3*p2p4 + (32*mt**3*p1p3*p2p4)/s - 
     $      (64*mt*p1p2*p1p3*p2p4)/s + 
     $      (32*mt*p1p3**2*p2p4)/s + 
     $      (16*mt*p1p3*p1p4*p2p4)/s + 
     $      (64*mt*p1p2*p2p3*p2p4)/s + 
     $      (32*mt*p1p3*p2p3*p2p4)/s + 
     $      (48*mt*p1p4*p2p3*p2p4)/s + 
     $      (128*mt*p1p3*p1p4*p2p3*p2p4)/s**2 - 
     $      (128*mt*p1p4*p2p3**2*p2p4)/s**2 + 
     $      (48*mt*p1p3*p2p4**2)/s + 
     $      (64*mt*p1p3**2*p2p4**2)/s**2 - 
     $      (128*mt*p1p3*p2p3*p2p4**2)/s**2)

       mat(3,32)=cu1*(16*mt**5 - 4*mt**3*s - 16*mt**3*p1p2 + 
     $      16*mt**3*p1p3 + 24*mt**3*p1p4 + 16*mt*p1p3*p1p4 - 
     $      (128*mt**3*p1p3*p1p4)/s + 
     $      (64*mt*p1p2*p1p3*p1p4)/s - 
     $      (64*mt*p1p3**2*p1p4)/s - (96*mt*p1p3*p1p4**2)/s + 
     $      (256*mt*p1p3**2*p1p4**2)/s**2 + 16*mt**3*p2p3 + 
     $      (32*mt**3*p1p4*p2p3)/s - 
     $      (64*mt*p1p3*p1p4*p2p3)/s - 
     $      (128*mt*p1p3*p1p4**2*p2p3)/s**2 + 8*mt**3*p2p4 + 
     $      (32*mt**3*p1p3*p2p4)/s - 
     $      (32*mt*p1p3*p1p4*p2p4)/s - 
     $      (128*mt*p1p3**2*p1p4*p2p4)/s**2) + 
     $      cu3*(16*mt**5 - 4*mt**3*s - 16*mt**3*p1p2 +16*mt**3*p1p3- 
     $      8*mt*p1p2*p1p4 + 16*mt*p1p3*p1p4 - 
     $      (64*mt**3*p1p3*p1p4)/s + 
     $      (64*mt*p1p2*p1p3*p1p4)/s - 
     $      (64*mt*p1p3**2*p1p4)/s + 16*mt**3*p2p3 + 
     $      (32*mt**3*p1p4*p2p3)/s - 
     $      (64*mt*p1p3*p1p4*p2p3)/s + 
     $      (16*mt*p1p4**2*p2p3)/s - 
     $      (128*mt*p1p3*p1p4**2*p2p3)/s**2 - 
     $      8*mt*p1p2*p2p4 + (32*mt**3*p1p3*p2p4)/s + 
     $      (16*mt*p1p3*p1p4*p2p4)/s - 
     $      (128*mt*p1p3**2*p1p4*p2p4)/s**2 - 
     $      (64*mt**3*p2p3*p2p4)/s + 
     $      (16*mt*p1p4*p2p3*p2p4)/s + 
     $      (256*mt*p1p3*p1p4*p2p3*p2p4)/s**2 + 
     $      (16*mt*p1p3*p2p4**2)/s) + 
     $      cu2*(-4*mt**3*s - 16*mt**3*p1p2 + 16*mt*p1p2**2 - 
     $      16*mt*p1p2*p1p3 - 8*mt*p1p2*p1p4 + 
     $      16*mt*p1p3*p1p4 + 16*mt**3*p2p3 - 
     $      16*mt*p1p2*p2p3 + (32*mt**3*p1p4*p2p3)/s - 
     $      (64*mt*p1p2*p1p4*p2p3)/s - 
     $      (32*mt*p1p3*p1p4*p2p3)/s + 
     $      (16*mt*p1p4**2*p2p3)/s + (32*mt*p1p4*p2p3**2)/s + 
     $      (64*mt*p1p4**2*p2p3**2)/s**2 - 8*mt*p1p2*p2p4 + 
     $      (32*mt**3*p1p3*p2p4)/s - 
     $      (64*mt*p1p2*p1p3*p2p4)/s + 
     $      (32*mt*p1p3**2*p2p4)/s + 
     $      (16*mt*p1p3*p1p4*p2p4)/s + 
     $      (64*mt*p1p2*p2p3*p2p4)/s + 
     $      (32*mt*p1p3*p2p3*p2p4)/s + 
     $      (16*mt*p1p4*p2p3*p2p4)/s + 
     $      (128*mt*p1p3*p1p4*p2p3*p2p4)/s**2 - 
     $      (128*mt*p1p4*p2p3**2*p2p4)/s**2 + 
     $      (16*mt*p1p3*p2p4**2)/s + 
     $      (64*mt*p1p3**2*p2p4**2)/s**2 - 
     $      (128*mt*p1p3*p2p3*p2p4**2)/s**2)

       do i=1,3
          do j=33,36
             mat(i,j)=0d0
          enddo
       end do

       mat(1,37)=cs1*(-2*mt**4*s**2+2*mt**2*s**2*p1p2-4*mt**4*s*p1p3+
     $      4*mt**2*s*p1p2*p1p3+4*mt**4*s*p1p4-
     $      4*mt**2*s*p1p2*p1p4-4*mt**2*s*p1p4*p2p3-
     $      8*mt**2*p1p3*p1p4*p2p3+8*mt**2*p1p4**2*p2p3-
     $      12*mt**2*s*p1p3*p2p4+8*mt**2*p1p3**2*p2p4-
     $      8*mt**2*p1p3*p1p4*p2p4+8*mt**2*s*p2p3*p2p4-
     $      8*s*p1p2*p2p3*p2p4+
     $      16*mt**2*p1p3*p2p3*p2p4-
     $      16*p1p2*p1p3*p2p3*p2p4-
     $      16*mt**2*p1p4*p2p3*p2p4+
     $      16*p1p2*p1p4*p2p3*p2p4+
     $      16*p1p4*p2p3**2*p2p4+
     $      (32*p1p3*p1p4*p2p3**2*p2p4)/s-
     $      (32*p1p4**2*p2p3**2*p2p4)/s+
     $      48*p1p3*p2p3*p2p4**2-
     $      (32*p1p3**2*p2p3*p2p4**2)/s+
     $      (32*p1p3*p1p4*p2p3*p2p4**2)/s)+
     $      cs2*(-2*mt**4*s**2+2*mt**2*s**2*p1p2+4*mt**4*s*p2p3-
     $      4*mt**2*s*p1p2*p2p3-4*mt**2*s*p1p4*p2p3+
     $      8*mt**2*p1p4*p2p3**2-4*mt**4*s*p2p4+
     $      4*mt**2*s*p1p2*p2p4-12*mt**2*s*p1p3*p2p4+
     $      8*mt**2*s*p2p3*p2p4-8*s*p1p2*p2p3*p2p4-
     $      8*mt**2*p1p3*p2p3*p2p4-
     $      8*mt**2*p1p4*p2p3*p2p4-16*mt**2*p2p3**2*p2p4+
     $      16*p1p2*p2p3**2*p2p4+16*p1p4*p2p3**2*p2p4-
     $      (32*p1p4*p2p3**3*p2p4)/s+8*mt**2*p1p3*p2p4**2+
     $      16*mt**2*p2p3*p2p4**2-16*p1p2*p2p3*p2p4**2+
     $      48*p1p3*p2p3*p2p4**2+
     $      (32*p1p3*p2p3**2*p2p4**2)/s+
     $      (32*p1p4*p2p3**2*p2p4**2)/s-
     $      (32*p1p3*p2p3*p2p4**3)/s)

       mat(2,37)=ct1*(-8*mt**2*s*p1p2**2+8*s*p1p2**3-
     $      4*mt**4*s*p1p3+12*mt**2*s*p1p2*p1p3+
     $      32*mt**2*p1p2*p1p4*p2p3-
     $      48*p1p2**2*p1p4*p2p3-
     $      24*mt**2*p1p3*p1p4*p2p3-
     $      (32*mt**2*p1p4**2*p2p3**2)/s+
     $      (96*p1p2*p1p4**2*p2p3**2)/s-
     $      (64*p1p4**3*p2p3**3)/s**2-8*s*p1p2**2*p2p4-
     $      8*mt**2*s*p1p3*p2p4+32*mt**2*p1p2*p1p3*p2p4-
     $      16*p1p2**2*p1p3*p2p4-8*mt**2*p1p3**2*p2p4+
     $      16*mt**2*p1p3*p2p3*p2p4-
     $      48*p1p2*p1p3*p2p3*p2p4+
     $      32*p1p2*p1p4*p2p3*p2p4-
     $      (64*mt**2*p1p3*p1p4*p2p3*p2p4)/s+
     $      (64*p1p2*p1p3*p1p4*p2p3*p2p4)/s+
     $      (96*p1p3*p1p4*p2p3**2*p2p4)/s-
     $      (32*p1p4**2*p2p3**2*p2p4)/s-
     $      (64*p1p3*p1p4**2*p2p3**2*p2p4)/s**2+
     $      32*p1p2*p1p3*p2p4**2-
     $      (32*mt**2*p1p3**2*p2p4**2)/s-
     $      (32*p1p2*p1p3**2*p2p4**2)/s+
     $      32*p1p3*p2p3*p2p4**2+
     $      (32*p1p3**2*p2p3*p2p4**2)/s-
     $      (64*p1p3*p1p4*p2p3*p2p4**2)/s+
     $      (64*p1p3**2*p1p4*p2p3*p2p4**2)/s**2-
     $      (32*p1p3**2*p2p4**3)/s+(64*p1p3**3*p2p4**3)/s**2)+
     $      ct3*(8*mt**4*s*p1p2-8*mt**2*s*p1p2**2-8*mt**4*s*p1p3-
     $      16*mt**4*p1p4*p2p3+32*mt**2*p1p2*p1p4*p2p3-
     $      (32*mt**2*p1p4**2*p2p3**2)/s-8*s*p1p2**2*p2p4-
     $      16*mt**4*p1p3*p2p4-8*mt**2*s*p1p3*p2p4-
     $      32*mt**2*p1p2*p2p3*p2p4+
     $      32*p1p2**2*p2p3*p2p4+
     $      64*mt**2*p1p3*p2p3*p2p4+
     $      32*p1p2*p1p4*p2p3*p2p4+
     $      (64*mt**2*p1p4*p2p3**2*p2p4)/s-
     $      (128*p1p2*p1p4*p2p3**2*p2p4)/s-
     $      (32*p1p4**2*p2p3**2*p2p4)/s+
     $      (128*p1p4**2*p2p3**3*p2p4)/s**2+
     $      32*p1p2*p1p3*p2p4**2+
     $      (32*mt**2*p1p3**2*p2p4**2)/s+
     $      32*p1p3*p2p3*p2p4**2+
     $      (64*mt**2*p1p3*p2p3*p2p4**2)/s-
     $      (64*p1p3*p1p4*p2p3*p2p4**2)/s-
     $      (128*p1p3*p2p3**2*p2p4**2)/s-
     $      (32*p1p3**2*p2p4**3)/s-
     $      (128*p1p3**2*p2p3*p2p4**3)/s**2)+
     $      ct2*(-8*mt**6*s+8*mt**4*s*p1p2-8*mt**4*s*p1p3-
     $      16*mt**4*p1p4*p2p3-4*mt**4*s*p2p4+
     $      12*mt**2*s*p1p2*p2p4+16*mt**4*p1p3*p2p4-
     $      8*mt**2*s*p1p3*p2p4+64*mt**4*p2p3*p2p4-
     $      64*mt**2*p1p2*p2p3*p2p4+
     $      64*mt**2*p1p3*p2p3*p2p4-
     $      24*mt**2*p1p4*p2p3*p2p4+
     $      (128*mt**2*p1p4*p2p3**2*p2p4)/s-
     $      8*mt**2*p1p3*p2p4**2+16*mt**2*p2p3*p2p4**2-
     $      48*p1p2*p2p3*p2p4**2+32*p1p3*p2p3*p2p4**2-
     $      (128*mt**2*p1p3*p2p3*p2p4**2)/s-
     $      (128*mt**2*p2p3**2*p2p4**2)/s+
     $      (128*p1p2*p2p3**2*p2p4**2)/s-
     $      (128*p1p3*p2p3**2*p2p4**2)/s+
     $      (96*p1p4*p2p3**2*p2p4**2)/s-
     $      (256*p1p4*p2p3**3*p2p4**2)/s**2+
     $      (32*p1p3*p2p3*p2p4**3)/s+
     $      (256*p1p3*p2p3**2*p2p4**3)/s**2)

       mat(3,37)=cu1*(2*mt**4*s**2-2*mt**2*s**2*p1p2-
     $      8*mt**2*s*p1p2**2+8*s*p1p2**3+
     $      8*mt**2*s*p1p2*p1p3-4*mt**4*s*p1p4+
     $      4*mt**2*s*p1p2*p1p4+4*mt**2*s*p1p4*p2p3+
     $      32*mt**2*p1p2*p1p4*p2p3-
     $      48*p1p2**2*p1p4*p2p3-
     $      16*mt**2*p1p3*p1p4*p2p3-8*mt**2*p1p4**2*p2p3-
     $      (32*mt**2*p1p4**2*p2p3**2)/s+
     $      (96*p1p2*p1p4**2*p2p3**2)/s-
     $      (64*p1p4**3*p2p3**3)/s**2-8*s*p1p2**2*p2p4+
     $      4*mt**2*s*p1p3*p2p4+32*mt**2*p1p2*p1p3*p2p4-
     $      16*p1p2**2*p1p3*p2p4-16*mt**2*p1p3**2*p2p4+
     $      8*mt**2*p1p3*p1p4*p2p4-8*mt**2*s*p2p3*p2p4+
     $      8*s*p1p2*p2p3*p2p4-
     $      32*p1p2*p1p3*p2p3*p2p4+
     $      16*mt**2*p1p4*p2p3*p2p4+
     $      16*p1p2*p1p4*p2p3*p2p4-
     $      (64*mt**2*p1p3*p1p4*p2p3*p2p4)/s+
     $      (64*p1p2*p1p3*p1p4*p2p3*p2p4)/s-
     $      16*p1p4*p2p3**2*p2p4+
     $      (64*p1p3*p1p4*p2p3**2*p2p4)/s-
     $      (64*p1p3*p1p4**2*p2p3**2*p2p4)/s**2+
     $      32*p1p2*p1p3*p2p4**2-
     $      (32*mt**2*p1p3**2*p2p4**2)/s-
     $      (32*p1p2*p1p3**2*p2p4**2)/s-
     $      16*p1p3*p2p3*p2p4**2+
     $      (64*p1p3**2*p2p3*p2p4**2)/s-
     $      (96*p1p3*p1p4*p2p3*p2p4**2)/s+
     $      (64*p1p3**2*p1p4*p2p3*p2p4**2)/s**2-
     $      (32*p1p3**2*p2p4**3)/s+(64*p1p3**3*p2p4**3)/s**2)+
     $      cu3*(2*mt**4*s**2+8*mt**4*s*p1p2-2*mt**2*s**2*p1p2-
     $      8*mt**2*s*p1p2**2+8*mt**2*s*p1p2*p1p3-
     $      16*mt**4*p1p4*p2p3+4*mt**2*s*p1p4*p2p3+
     $      32*mt**2*p1p2*p1p4*p2p3-
     $      16*mt**2*p1p3*p1p4*p2p3-
     $      (32*mt**2*p1p4**2*p2p3**2)/s+8*mt**2*s*p1p2*p2p4-
     $      16*mt**4*p1p3*p2p4+4*mt**2*s*p1p3*p2p4-
     $      16*mt**2*p1p3**2*p2p4-8*mt**2*s*p2p3*p2p4-
     $      32*mt**2*p1p2*p2p3*p2p4+
     $      8*s*p1p2*p2p3*p2p4+32*p1p2**2*p2p3*p2p4-
     $      32*p1p2*p1p3*p2p3*p2p4-
     $      16*mt**2*p1p4*p2p3*p2p4-
     $      16*p1p4*p2p3**2*p2p4+
     $      (64*mt**2*p1p4*p2p3**2*p2p4)/s-
     $      (128*p1p2*p1p4*p2p3**2*p2p4)/s+
     $      (64*p1p3*p1p4*p2p3**2*p2p4)/s+
     $      (128*p1p4**2*p2p3**3*p2p4)/s**2-
     $      16*mt**2*p1p3*p2p4**2+
     $      (32*mt**2*p1p3**2*p2p4**2)/s-
     $      32*p1p2*p2p3*p2p4**2-
     $      16*p1p3*p2p3*p2p4**2+
     $      (64*mt**2*p1p3*p2p3*p2p4**2)/s+
     $      (64*p1p3**2*p2p3*p2p4**2)/s+
     $      (64*p1p4*p2p3**2*p2p4**2)/s+
     $      (64*p1p3*p2p3*p2p4**3)/s-
     $      (128*p1p3**2*p2p3*p2p4**3)/s**2)+
     $      cu2*(-8*mt**6*s+2*mt**4*s**2+8*mt**4*s*p1p2-
     $      2*mt**2*s**2*p1p2-8*mt**4*s*p1p3-4*mt**4*s*p2p3+
     $      4*mt**2*s*p1p2*p2p3-16*mt**4*p1p4*p2p3+
     $      4*mt**2*s*p1p4*p2p3-8*mt**2*p1p4*p2p3**2+
     $      8*mt**2*s*p1p2*p2p4+16*mt**4*p1p3*p2p4+
     $      4*mt**2*s*p1p3*p2p4+64*mt**4*p2p3*p2p4-
     $      8*mt**2*s*p2p3*p2p4-64*mt**2*p1p2*p2p3*p2p4+
     $      8*s*p1p2*p2p3*p2p4+
     $      72*mt**2*p1p3*p2p3*p2p4-
     $      16*mt**2*p1p4*p2p3*p2p4+16*mt**2*p2p3**2*p2p4-
     $      16*p1p2*p2p3**2*p2p4-16*p1p4*p2p3**2*p2p4+
     $      (128*mt**2*p1p4*p2p3**2*p2p4)/s+
     $      (32*p1p4*p2p3**3*p2p4)/s-
     $      16*mt**2*p1p3*p2p4**2-32*p1p2*p2p3*p2p4**2-
     $      16*p1p3*p2p3*p2p4**2-
     $      (128*mt**2*p1p3*p2p3*p2p4**2)/s-
     $      (128*mt**2*p2p3**2*p2p4**2)/s+
     $      (128*p1p2*p2p3**2*p2p4**2)/s-
     $      (160*p1p3*p2p3**2*p2p4**2)/s+
     $      (64*p1p4*p2p3**2*p2p4**2)/s-
     $      (256*p1p4*p2p3**3*p2p4**2)/s**2+
     $      (64*p1p3*p2p3*p2p4**3)/s+
     $      (256*p1p3*p2p3**2*p2p4**3)/s**2)

       mat(1,38)=cs1*(-2*mt**2*s**2*p1p2+2*s**2*p1p2**2-
     $      4*mt**2*s*p1p2*p1p3+4*s*p1p2**2*p1p3+
     $      4*mt**2*s*p1p2*p1p4-4*s*p1p2**2*p1p4+
     $      4*mt**2*s*p1p4*p2p3-8*s*p1p2*p1p4*p2p3+
     $      8*mt**2*p1p3*p1p4*p2p3-
     $      16*p1p2*p1p3*p1p4*p2p3-
     $      8*mt**2*p1p4**2*p2p3+16*p1p2*p1p4**2*p2p3+
     $      8*p1p4**2*p2p3**2+
     $      (16*p1p3*p1p4**2*p2p3**2)/s-
     $      (16*p1p4**3*p2p3**2)/s+4*mt**2*s*p1p3*p2p4-
     $      16*s*p1p2*p1p3*p2p4+8*mt**2*p1p3**2*p2p4-
     $      8*mt**2*p1p3*p1p4*p2p4+
     $      32*p1p3*p1p4*p2p3*p2p4+
     $      24*p1p3**2*p2p4**2-(16*p1p3**3*p2p4**2)/s+
     $      (16*p1p3**2*p1p4*p2p4**2)/s)+
     $      cs2*(-2*mt**2*s**2*p1p2+2*s**2*p1p2**2+
     $      4*mt**2*s*p1p2*p2p3-4*s*p1p2**2*p2p3+
     $      4*mt**2*s*p1p4*p2p3-8*s*p1p2*p1p4*p2p3-
     $      8*mt**2*p1p4*p2p3**2+16*p1p2*p1p4*p2p3**2+
     $      8*p1p4**2*p2p3**2-(16*p1p4**2*p2p3**3)/s-
     $      4*mt**2*s*p1p2*p2p4+4*s*p1p2**2*p2p4+
     $      4*mt**2*s*p1p3*p2p4-16*s*p1p2*p1p3*p2p4-
     $      8*mt**2*p1p3*p2p3*p2p4+
     $      8*mt**2*p1p4*p2p3*p2p4-
     $      16*p1p2*p1p4*p2p3*p2p4+
     $      32*p1p3*p1p4*p2p3*p2p4+
     $      (16*p1p4**2*p2p3**2*p2p4)/s+
     $      8*mt**2*p1p3*p2p4**2+24*p1p3**2*p2p4**2+
     $      (16*p1p3**2*p2p3*p2p4**2)/s-
     $      (16*p1p3**2*p2p4**3)/s)

       mat(2,38)=ct1*(-8*mt**4*s*p1p2+8*mt**2*s*p1p2**2+
     $      4*mt**4*s*p1p3-4*mt**2*s*p1p2*p1p3+
     $      8*s*p1p2**2*p1p3+32*mt**2*p1p2*p1p3*p1p4-
     $      32*p1p2**2*p1p3*p1p4-16*mt**2*p1p3**2*p1p4+
     $      16*mt**4*p1p4*p2p3-32*mt**2*p1p2*p1p4*p2p3+
     $      8*mt**2*p1p3*p1p4*p2p3-
     $      32*p1p2*p1p3*p1p4*p2p3-
     $      (64*mt**2*p1p3*p1p4**2*p2p3)/s+
     $      (128*p1p2*p1p3*p1p4**2*p2p3)/s+
     $      (32*mt**2*p1p4**2*p2p3**2)/s+
     $      (32*p1p3*p1p4**2*p2p3**2)/s-
     $      (128*p1p3*p1p4**3*p2p3**2)/s**2-
     $      8*mt**2*s*p1p2*p2p4+16*mt**4*p1p3*p2p4-
     $      8*s*p1p2*p1p3*p2p4+8*mt**2*p1p3**2*p2p4-
     $      16*p1p2*p1p3**2*p2p4+
     $      32*p1p2*p1p3*p1p4*p2p4-
     $      (64*mt**2*p1p3**2*p1p4*p2p4)/s-
     $      16*mt**2*p1p3*p2p3*p2p4+
     $      16*mt**2*p1p4*p2p3*p2p4+
     $      16*p1p3*p1p4*p2p3*p2p4+
     $      (96*p1p3**2*p1p4*p2p3*p2p4)/s-
     $      (64*p1p3*p1p4**2*p2p3*p2p4)/s+
     $      16*mt**2*p1p3*p2p4**2+16*p1p3**2*p2p4**2-
     $      (32*mt**2*p1p3**2*p2p4**2)/s-
     $      (64*p1p3**2*p1p4*p2p4**2)/s+
     $      (128*p1p3**3*p1p4*p2p4**2)/s**2)+
     $      ct3*(8*mt**6*s-8*mt**4*s*p1p2-8*mt**2*s*p1p2*p1p3-
     $      32*mt**4*p1p3*p1p4+32*mt**2*p1p2*p1p3*p1p4+
     $      16*mt**4*p1p4*p2p3+16*mt**2*p1p3*p1p4*p2p3-
     $      (64*mt**2*p1p3*p1p4**2*p2p3)/s-
     $      8*mt**2*s*p1p2*p2p4-16*mt**4*p1p3*p2p4-
     $      8*s*p1p2*p1p3*p2p4+16*mt**2*p1p3**2*p2p4+
     $      32*p1p2*p1p3*p1p4*p2p4+
     $      (64*mt**2*p1p3**2*p1p4*p2p4)/s-
     $      32*mt**4*p2p3*p2p4+32*mt**2*p1p2*p2p3*p2p4+
     $      32*p1p2*p1p3*p2p3*p2p4+
     $      16*mt**2*p1p4*p2p3*p2p4+
     $      16*p1p3*p1p4*p2p3*p2p4+
     $      (128*mt**2*p1p3*p1p4*p2p3*p2p4)/s-
     $      (128*p1p2*p1p3*p1p4*p2p3*p2p4)/s-
     $      (64*p1p3*p1p4**2*p2p3*p2p4)/s-
     $      (64*mt**2*p1p4*p2p3**2*p2p4)/s-
     $      (64*p1p3*p1p4*p2p3**2*p2p4)/s+
     $      (256*p1p3*p1p4**2*p2p3**2*p2p4)/s**2+
     $      16*mt**2*p1p3*p2p4**2+16*p1p3**2*p2p4**2-
     $      (64*p1p3**2*p1p4*p2p4**2)/s+
     $      (64*mt**2*p1p3*p2p3*p2p4**2)/s-
     $      (64*p1p3**2*p2p3*p2p4**2)/s-
     $      (256*p1p3**2*p1p4*p2p3*p2p4**2)/s**2)+
     $      ct2*(-8*mt**4*s*p1p2+8*mt**2*s*p1p2**2-
     $      8*mt**2*s*p1p2*p1p3+16*mt**4*p1p4*p2p3-
     $      32*mt**2*p1p2*p1p4*p2p3+
     $      16*mt**2*p1p3*p1p4*p2p3+
     $      (32*mt**2*p1p4**2*p2p3**2)/s+4*mt**4*s*p2p4-
     $      4*mt**2*s*p1p2*p2p4+8*s*p1p2**2*p2p4+
     $      16*mt**4*p1p3*p2p4-8*s*p1p2*p1p3*p2p4+
     $      16*mt**2*p1p3**2*p2p4-16*mt**2*p1p3*p1p4*p2p4+
     $      32*mt**2*p1p2*p2p3*p2p4-
     $      32*p1p2**2*p2p3*p2p4+
     $      32*p1p2*p1p3*p2p3*p2p4+
     $      8*mt**2*p1p4*p2p3*p2p4-
     $      32*p1p2*p1p4*p2p3*p2p4+
     $      16*p1p3*p1p4*p2p3*p2p4-
     $      (64*mt**2*p1p4*p2p3**2*p2p4)/s+
     $      (128*p1p2*p1p4*p2p3**2*p2p4)/s-
     $      (64*p1p3*p1p4*p2p3**2*p2p4)/s+
     $      (32*p1p4**2*p2p3**2*p2p4)/s-
     $      (128*p1p4**2*p2p3**3*p2p4)/s**2+
     $      8*mt**2*p1p3*p2p4**2-16*p1p2*p1p3*p2p4**2+
     $      16*p1p3**2*p2p4**2-(32*mt**2*p1p3**2*p2p4**2)/s-
     $      16*mt**2*p2p3*p2p4**2-
     $      (64*mt**2*p1p3*p2p3*p2p4**2)/s-
     $      (64*p1p3**2*p2p3*p2p4**2)/s+
     $      (96*p1p3*p1p4*p2p3*p2p4**2)/s+
     $      (128*p1p3**2*p2p3*p2p4**3)/s**2)

       mat(3,38)=cu1*(2*mt**4*s**2-8*mt**4*s*p1p2+2*mt**2*s**2*p1p2+
     $      8*mt**2*s*p1p2**2-4*s**2*p1p2**2+
     $      8*s*p1p2**2*p1p3-4*mt**4*s*p1p4-
     $      4*mt**2*s*p1p2*p1p4+8*s*p1p2**2*p1p4-
     $      8*mt**2*s*p1p3*p1p4+32*mt**2*p1p2*p1p3*p1p4-
     $      32*p1p2**2*p1p3*p1p4+16*mt**2*p1p3*p1p4**2+
     $      16*mt**4*p1p4*p2p3-4*mt**2*s*p1p4*p2p3-
     $      32*mt**2*p1p2*p1p4*p2p3+
     $      16*s*p1p2*p1p4*p2p3-
     $      32*p1p2*p1p3*p1p4*p2p3+
     $      8*mt**2*p1p4**2*p2p3-32*p1p2*p1p4**2*p2p3-
     $      (64*mt**2*p1p3*p1p4**2*p2p3)/s+
     $      (128*p1p2*p1p3*p1p4**2*p2p3)/s-
     $      16*p1p4**2*p2p3**2+(32*mt**2*p1p4**2*p2p3**2)/s+
     $      (32*p1p3*p1p4**2*p2p3**2)/s+
     $      (32*p1p4**3*p2p3**2)/s-
     $      (128*p1p3*p1p4**3*p2p3**2)/s**2-
     $      8*mt**2*s*p1p2*p2p4+16*mt**4*p1p3*p2p4-
     $      4*mt**2*s*p1p3*p2p4+16*s*p1p2*p1p3*p2p4-
     $      32*p1p2*p1p3**2*p2p4+
     $      8*mt**2*p1p3*p1p4*p2p4+
     $      16*p1p2*p1p3*p1p4*p2p4-
     $      (64*mt**2*p1p3**2*p1p4*p2p4)/s-
     $      8*mt**2*s*p2p3*p2p4+32*mt**2*p1p4*p2p3*p2p4+
     $      (64*p1p3**2*p1p4*p2p3*p2p4)/s-
     $      (96*p1p3*p1p4**2*p2p3*p2p4)/s+
     $      16*mt**2*p1p3*p2p4**2-16*p1p3**2*p2p4**2-
     $      (32*mt**2*p1p3**2*p2p4**2)/s+
     $      (32*p1p3**3*p2p4**2)/s-
     $      (64*p1p3**2*p1p4*p2p4**2)/s+
     $      (128*p1p3**3*p1p4*p2p4**2)/s**2)+
     $      cu3*(2*mt**4*s**2+2*mt**2*s**2*p1p2+8*mt**2*s*p1p2**2-
     $      4*s**2*p1p2**2-8*s*p1p2**3+8*s*p1p2**2*p1p3-
     $      8*mt**2*s*p1p3*p1p4-4*mt**2*s*p1p4*p2p3-
     $      32*mt**2*p1p2*p1p4*p2p3+
     $      16*s*p1p2*p1p4*p2p3+48*p1p2**2*p1p4*p2p3-
     $      32*p1p2*p1p3*p1p4*p2p3-
     $      16*p1p4**2*p2p3**2+(32*mt**2*p1p4**2*p2p3**2)/s-
     $      (96*p1p2*p1p4**2*p2p3**2)/s+
     $      (32*p1p3*p1p4**2*p2p3**2)/s+
     $      (64*p1p4**3*p2p3**3)/s**2+8*s*p1p2**2*p2p4-
     $      4*mt**2*s*p1p3*p2p4-32*mt**2*p1p2*p1p3*p2p4+
     $      16*s*p1p2*p1p3*p2p4+16*p1p2**2*p1p3*p2p4-
     $      32*p1p2*p1p3**2*p2p4-8*mt**2*s*p2p3*p2p4-
     $      32*p1p2*p1p4*p2p3*p2p4+
     $      (64*mt**2*p1p3*p1p4*p2p3*p2p4)/s-
     $      (64*p1p2*p1p3*p1p4*p2p3*p2p4)/s+
     $      (64*p1p3**2*p1p4*p2p3*p2p4)/s+
     $      (32*p1p4**2*p2p3**2*p2p4)/s+
     $      (64*p1p3*p1p4**2*p2p3**2*p2p4)/s**2-
     $      32*p1p2*p1p3*p2p4**2-16*p1p3**2*p2p4**2+
     $      (32*mt**2*p1p3**2*p2p4**2)/s+
     $      (32*p1p2*p1p3**2*p2p4**2)/s+
     $      (32*p1p3**3*p2p4**2)/s+
     $      (64*p1p3*p1p4*p2p3*p2p4**2)/s-
     $      (64*p1p3**2*p1p4*p2p3*p2p4**2)/s**2+
     $      (32*p1p3**2*p2p4**3)/s-(64*p1p3**3*p2p4**3)/s**2)+
     $      cu2*(2*mt**4*s**2-8*mt**4*s*p1p2+2*mt**2*s**2*p1p2+
     $      8*mt**2*s*p1p2**2-4*s**2*p1p2**2-
     $      8*mt**2*s*p1p2*p1p3-8*mt**2*s*p1p3*p1p4-
     $      4*mt**4*s*p2p3-4*mt**2*s*p1p2*p2p3+
     $      8*s*p1p2**2*p2p3+16*mt**4*p1p4*p2p3-
     $      4*mt**2*s*p1p4*p2p3-32*mt**2*p1p2*p1p4*p2p3+
     $      16*s*p1p2*p1p4*p2p3+
     $      32*mt**2*p1p3*p1p4*p2p3+8*mt**2*p1p4*p2p3**2-
     $      32*p1p2*p1p4*p2p3**2-16*p1p4**2*p2p3**2+
     $      (32*mt**2*p1p4**2*p2p3**2)/s+
     $      (32*p1p4**2*p2p3**3)/s+8*s*p1p2**2*p2p4+
     $      16*mt**4*p1p3*p2p4-4*mt**2*s*p1p3*p2p4+
     $      16*s*p1p2*p1p3*p2p4+16*mt**2*p1p3**2*p2p4-
     $      8*mt**2*s*p2p3*p2p4+32*mt**2*p1p2*p2p3*p2p4-
     $      32*p1p2**2*p2p3*p2p4+
     $      8*mt**2*p1p3*p2p3*p2p4+
     $      16*p1p2*p1p3*p2p3*p2p4-
     $      32*p1p2*p1p4*p2p3*p2p4+
     $      16*mt**2*p2p3**2*p2p4-
     $      (64*mt**2*p1p4*p2p3**2*p2p4)/s+
     $      (128*p1p2*p1p4*p2p3**2*p2p4)/s-
     $      (96*p1p3*p1p4*p2p3**2*p2p4)/s+
     $      (32*p1p4**2*p2p3**2*p2p4)/s-
     $      (128*p1p4**2*p2p3**3*p2p4)/s**2-
     $      32*p1p2*p1p3*p2p4**2-16*p1p3**2*p2p4**2-
     $      (32*mt**2*p1p3**2*p2p4**2)/s-
     $      (64*mt**2*p1p3*p2p3*p2p4**2)/s-
     $      (64*p1p3**2*p2p3*p2p4**2)/s+
     $      (64*p1p3*p1p4*p2p3*p2p4**2)/s+
     $      (32*p1p3**2*p2p4**3)/s+
     $      (128*p1p3**2*p2p3*p2p4**3)/s**2)

       do i=1,3
          mat(i,39)=0d0
       end do

       mat(1,40)=cs1*(-2*mt**2*s**2*p1p2+2*s**2*p1p2**2-
     $      4*mt**2*s*p1p2*p1p3+4*s*p1p2**2*p1p3+
     $      4*mt**2*s*p1p2*p1p4-4*s*p1p2**2*p1p4+
     $      4*mt**2*s*p1p4*p2p3-8*s*p1p2*p1p4*p2p3+
     $      8*mt**2*p1p3*p1p4*p2p3-
     $      16*p1p2*p1p3*p1p4*p2p3-
     $      8*mt**2*p1p4**2*p2p3+16*p1p2*p1p4**2*p2p3+
     $      8*p1p4**2*p2p3**2+
     $      (16*p1p3*p1p4**2*p2p3**2)/s-
     $      (16*p1p4**3*p2p3**2)/s+4*mt**2*s*p1p3*p2p4-
     $      16*s*p1p2*p1p3*p2p4+8*mt**2*p1p3**2*p2p4-
     $      8*mt**2*p1p3*p1p4*p2p4+
     $      32*p1p3*p1p4*p2p3*p2p4+
     $      24*p1p3**2*p2p4**2-(16*p1p3**3*p2p4**2)/s+
     $      (16*p1p3**2*p1p4*p2p4**2)/s)+
     $      cs2*(-2*mt**2*s**2*p1p2+2*s**2*p1p2**2+
     $      4*mt**2*s*p1p2*p2p3-4*s*p1p2**2*p2p3+
     $      4*mt**2*s*p1p4*p2p3-8*s*p1p2*p1p4*p2p3-
     $      8*mt**2*p1p4*p2p3**2+16*p1p2*p1p4*p2p3**2+
     $      8*p1p4**2*p2p3**2-(16*p1p4**2*p2p3**3)/s-
     $      4*mt**2*s*p1p2*p2p4+4*s*p1p2**2*p2p4+
     $      4*mt**2*s*p1p3*p2p4-16*s*p1p2*p1p3*p2p4-
     $      8*mt**2*p1p3*p2p3*p2p4+
     $      8*mt**2*p1p4*p2p3*p2p4-
     $      16*p1p2*p1p4*p2p3*p2p4+
     $      32*p1p3*p1p4*p2p3*p2p4+
     $      (16*p1p4**2*p2p3**2*p2p4)/s+
     $      8*mt**2*p1p3*p2p4**2+24*p1p3**2*p2p4**2+
     $      (16*p1p3**2*p2p3*p2p4**2)/s-
     $      (16*p1p3**2*p2p4**3)/s)

       mat(2,40)=ct1*(-8*mt**4*s*p1p2+8*mt**2*s*p1p2**2+
     $      4*mt**4*s*p1p3-4*mt**2*s*p1p2*p1p3+
     $      8*s*p1p2**2*p1p3+32*mt**2*p1p2*p1p3*p1p4-
     $      32*p1p2**2*p1p3*p1p4-16*mt**2*p1p3**2*p1p4+
     $      16*mt**4*p1p4*p2p3-32*mt**2*p1p2*p1p4*p2p3+
     $      8*mt**2*p1p3*p1p4*p2p3-
     $      32*p1p2*p1p3*p1p4*p2p3-
     $      (64*mt**2*p1p3*p1p4**2*p2p3)/s+
     $      (128*p1p2*p1p3*p1p4**2*p2p3)/s+
     $      (32*mt**2*p1p4**2*p2p3**2)/s+
     $      (32*p1p3*p1p4**2*p2p3**2)/s-
     $      (128*p1p3*p1p4**3*p2p3**2)/s**2-
     $      8*mt**2*s*p1p2*p2p4+16*mt**4*p1p3*p2p4-
     $      8*s*p1p2*p1p3*p2p4+8*mt**2*p1p3**2*p2p4-
     $      16*p1p2*p1p3**2*p2p4+
     $      32*p1p2*p1p3*p1p4*p2p4-
     $      (64*mt**2*p1p3**2*p1p4*p2p4)/s-
     $      16*mt**2*p1p3*p2p3*p2p4+
     $      16*mt**2*p1p4*p2p3*p2p4+
     $      16*p1p3*p1p4*p2p3*p2p4+
     $      (96*p1p3**2*p1p4*p2p3*p2p4)/s-
     $      (64*p1p3*p1p4**2*p2p3*p2p4)/s+
     $      16*mt**2*p1p3*p2p4**2+16*p1p3**2*p2p4**2-
     $      (32*mt**2*p1p3**2*p2p4**2)/s-
     $      (64*p1p3**2*p1p4*p2p4**2)/s+
     $      (128*p1p3**3*p1p4*p2p4**2)/s**2)+
     $      ct3*(8*mt**2*s*p1p2**2-8*s*p1p2**3-
     $      8*mt**2*s*p1p2*p1p3-32*mt**2*p1p2*p1p4*p2p3+
     $      48*p1p2**2*p1p4*p2p3+
     $      16*mt**2*p1p3*p1p4*p2p3+
     $      (32*mt**2*p1p4**2*p2p3**2)/s-
     $      (96*p1p2*p1p4**2*p2p3**2)/s+
     $      (64*p1p4**3*p2p3**3)/s**2-8*mt**2*s*p1p2*p2p4-
     $      32*mt**2*p1p2*p1p3*p2p4-
     $      8*s*p1p2*p1p3*p2p4+16*p1p2**2*p1p3*p2p4+
     $      16*mt**2*p1p3**2*p2p4+
     $      32*p1p2*p1p3*p1p4*p2p4+
     $      32*p1p2*p1p3*p2p3*p2p4+
     $      16*mt**2*p1p4*p2p3*p2p4+
     $      16*p1p3*p1p4*p2p3*p2p4+
     $      (64*mt**2*p1p3*p1p4*p2p3*p2p4)/s-
     $      (64*p1p2*p1p3*p1p4*p2p3*p2p4)/s-
     $      (64*p1p3*p1p4**2*p2p3*p2p4)/s-
     $      (64*p1p3*p1p4*p2p3**2*p2p4)/s+
     $      (64*p1p3*p1p4**2*p2p3**2*p2p4)/s**2+
     $      16*mt**2*p1p3*p2p4**2+16*p1p3**2*p2p4**2+
     $      (32*mt**2*p1p3**2*p2p4**2)/s+
     $      (32*p1p2*p1p3**2*p2p4**2)/s-
     $      (64*p1p3**2*p1p4*p2p4**2)/s-
     $      (64*p1p3**2*p2p3*p2p4**2)/s-
     $      (64*p1p3**2*p1p4*p2p3*p2p4**2)/s**2-
     $      (64*p1p3**3*p2p4**3)/s**2)+
     $      ct2*(-8*mt**4*s*p1p2+8*mt**2*s*p1p2**2-
     $      8*mt**2*s*p1p2*p1p3+16*mt**4*p1p4*p2p3-
     $      32*mt**2*p1p2*p1p4*p2p3+
     $      16*mt**2*p1p3*p1p4*p2p3+
     $      (32*mt**2*p1p4**2*p2p3**2)/s+4*mt**4*s*p2p4-
     $      4*mt**2*s*p1p2*p2p4+8*s*p1p2**2*p2p4+
     $      16*mt**4*p1p3*p2p4-8*s*p1p2*p1p3*p2p4+
     $      16*mt**2*p1p3**2*p2p4-16*mt**2*p1p3*p1p4*p2p4+
     $      32*mt**2*p1p2*p2p3*p2p4-
     $      32*p1p2**2*p2p3*p2p4+
     $      32*p1p2*p1p3*p2p3*p2p4+
     $      8*mt**2*p1p4*p2p3*p2p4-
     $      32*p1p2*p1p4*p2p3*p2p4+
     $      16*p1p3*p1p4*p2p3*p2p4-
     $      (64*mt**2*p1p4*p2p3**2*p2p4)/s+
     $      (128*p1p2*p1p4*p2p3**2*p2p4)/s-
     $      (64*p1p3*p1p4*p2p3**2*p2p4)/s+
     $      (32*p1p4**2*p2p3**2*p2p4)/s-
     $      (128*p1p4**2*p2p3**3*p2p4)/s**2+
     $      8*mt**2*p1p3*p2p4**2-16*p1p2*p1p3*p2p4**2+
     $      16*p1p3**2*p2p4**2-(32*mt**2*p1p3**2*p2p4**2)/s-
     $      16*mt**2*p2p3*p2p4**2-
     $      (64*mt**2*p1p3*p2p3*p2p4**2)/s-
     $      (64*p1p3**2*p2p3*p2p4**2)/s+
     $      (96*p1p3*p1p4*p2p3*p2p4**2)/s+
     $      (128*p1p3**2*p2p3*p2p4**3)/s**2)

       mat(3,40)=cu1*(-2*mt**4*s**2-8*mt**4*s*p1p2+2*mt**2*s**2*p1p2+
     $      8*mt**2*s*p1p2**2+8*mt**4*s*p1p3+4*mt**4*s*p1p4-
     $      4*mt**2*s*p1p2*p1p4+8*mt**2*s*p1p3*p1p4+
     $      32*mt**2*p1p2*p1p3*p1p4-
     $      32*p1p2**2*p1p3*p1p4-32*mt**2*p1p3**2*p1p4-
     $      16*mt**2*p1p3*p1p4**2+16*mt**4*p1p4*p2p3-
     $      4*mt**2*s*p1p4*p2p3-32*mt**2*p1p2*p1p4*p2p3+
     $      8*mt**2*p1p4**2*p2p3-
     $      (64*mt**2*p1p3*p1p4**2*p2p3)/s+
     $      (128*p1p2*p1p3*p1p4**2*p2p3)/s+
     $      (32*mt**2*p1p4**2*p2p3**2)/s-
     $      (128*p1p3*p1p4**3*p2p3**2)/s**2-
     $      8*mt**2*s*p1p2*p2p4+16*mt**4*p1p3*p2p4-
     $      4*mt**2*s*p1p3*p2p4+8*mt**2*p1p3*p1p4*p2p4+
     $      48*p1p2*p1p3*p1p4*p2p4-
     $      (64*mt**2*p1p3**2*p1p4*p2p4)/s+
     $      8*mt**2*s*p2p3*p2p4-32*mt**2*p1p3*p2p3*p2p4-
     $      32*p1p3*p1p4*p2p3*p2p4+
     $      (128*p1p3**2*p1p4*p2p3*p2p4)/s-
     $      (32*p1p3*p1p4**2*p2p3*p2p4)/s+
     $      16*mt**2*p1p3*p2p4**2-
     $      (32*mt**2*p1p3**2*p2p4**2)/s-
     $      (96*p1p3**2*p1p4*p2p4**2)/s+
     $      (128*p1p3**3*p1p4*p2p4**2)/s**2)+
     $      cu3*(8*mt**6*s-2*mt**4*s**2-8*mt**4*s*p1p2+
     $      2*mt**2*s**2*p1p2+8*mt**4*s*p1p3-
     $      32*mt**4*p1p3*p1p4+8*mt**2*s*p1p3*p1p4+
     $      32*mt**2*p1p2*p1p3*p1p4-
     $      32*mt**2*p1p3**2*p1p4+16*mt**4*p1p4*p2p3-
     $      4*mt**2*s*p1p4*p2p3-
     $      (64*mt**2*p1p3*p1p4**2*p2p3)/s+8*mt**4*s*p2p4-
     $      16*mt**4*p1p3*p2p4-4*mt**2*s*p1p3*p2p4-
     $      32*mt**2*p1p3*p1p4*p2p4+
     $      (64*mt**2*p1p3**2*p1p4*p2p4)/s-
     $      32*mt**4*p2p3*p2p4+8*mt**2*s*p2p3*p2p4+
     $      32*mt**2*p1p2*p2p3*p2p4-
     $      32*mt**2*p1p3*p2p3*p2p4-
     $      32*p1p3*p1p4*p2p3*p2p4+
     $      (128*mt**2*p1p3*p1p4*p2p3*p2p4)/s-
     $      (128*p1p2*p1p3*p1p4*p2p3*p2p4)/s+
     $      (128*p1p3**2*p1p4*p2p3*p2p4)/s-
     $      (64*mt**2*p1p4*p2p3**2*p2p4)/s+
     $      (256*p1p3*p1p4**2*p2p3**2*p2p4)/s**2-
     $      32*mt**2*p2p3*p2p4**2+
     $      (64*mt**2*p1p3*p2p3*p2p4**2)/s+
     $      (128*p1p3*p1p4*p2p3*p2p4**2)/s-
     $      (256*p1p3**2*p1p4*p2p3*p2p4**2)/s**2)+
     $      cu2*(-2*mt**4*s**2-8*mt**4*s*p1p2+2*mt**2*s**2*p1p2+
     $      8*mt**2*s*p1p2**2-8*mt**2*s*p1p2*p1p3+
     $      8*mt**2*s*p1p3*p1p4+4*mt**4*s*p2p3-
     $      4*mt**2*s*p1p2*p2p3+16*mt**4*p1p4*p2p3-
     $      4*mt**2*s*p1p4*p2p3-32*mt**2*p1p2*p1p4*p2p3+
     $      8*mt**2*p1p4*p2p3**2+(32*mt**2*p1p4**2*p2p3**2)/s+
     $      8*mt**4*s*p2p4+16*mt**4*p1p3*p2p4-
     $      4*mt**2*s*p1p3*p2p4+16*mt**2*p1p3**2*p2p4-
     $      32*mt**2*p1p3*p1p4*p2p4+8*mt**2*s*p2p3*p2p4+
     $      32*mt**2*p1p2*p2p3*p2p4-
     $      32*p1p2**2*p2p3*p2p4+
     $      8*mt**2*p1p3*p2p3*p2p4+
     $      48*p1p2*p1p3*p2p3*p2p4-
     $      32*p1p3*p1p4*p2p3*p2p4-
     $      16*mt**2*p2p3**2*p2p4-
     $      (64*mt**2*p1p4*p2p3**2*p2p4)/s+
     $      (128*p1p2*p1p4*p2p3**2*p2p4)/s-
     $      (32*p1p3*p1p4*p2p3**2*p2p4)/s-
     $      (128*p1p4**2*p2p3**3*p2p4)/s**2-
     $      (32*mt**2*p1p3**2*p2p4**2)/s-
     $      32*mt**2*p2p3*p2p4**2-
     $      (64*mt**2*p1p3*p2p3*p2p4**2)/s-
     $      (96*p1p3**2*p2p3*p2p4**2)/s+
     $      (128*p1p3*p1p4*p2p3*p2p4**2)/s+
     $      (128*p1p3**2*p2p3*p2p4**3)/s**2)


       mat(1,41)=cs1*(-2*mt**4*s**2+2*mt**2*s**2*p1p2-4*mt**4*s*p1p3 + 
     $      4*mt**2*s*p1p2*p1p3 + 4*mt**4*s*p1p4 - 
     $      4*mt**2*s*p1p2*p1p4 + 8*mt**2*s*p1p3*p1p4 - 
     $      8*s*p1p2*p1p3*p1p4 + 16*mt**2*p1p3**2*p1p4 - 
     $      16*p1p2*p1p3**2*p1p4 - 16*mt**2*p1p3*p1p4**2 + 
     $      16*p1p2*p1p3*p1p4**2 - 4*mt**2*s*p1p4*p2p3 - 
     $      8*mt**2*p1p3*p1p4*p2p3 + 8*mt**2*p1p4**2*p2p3 + 
     $      16*p1p3*p1p4**2*p2p3 + 
     $      (32*p1p3**2*p1p4**2*p2p3)/s - 
     $      (32*p1p3*p1p4**3*p2p3)/s - 
     $      12*mt**2*s*p1p3*p2p4 + 8*mt**2*p1p3**2*p2p4 - 
     $      8*mt**2*p1p3*p1p4*p2p4 + 
     $      48*p1p3**2*p1p4*p2p4 - 
     $      (32*p1p3**3*p1p4*p2p4)/s + 
     $      (32*p1p3**2*p1p4**2*p2p4)/s) + 
     $      cs2*(-2*mt**4*s**2 + 2*mt**2*s**2*p1p2 + 
     $      8*mt**2*s*p1p3*p1p4 - 8*s*p1p2*p1p3*p1p4 + 
     $      4*mt**4*s*p2p3 - 4*mt**2*s*p1p2*p2p3 - 
     $      4*mt**2*s*p1p4*p2p3 - 16*mt**2*p1p3*p1p4*p2p3 + 
     $      16*p1p2*p1p3*p1p4*p2p3 + 
     $      16*p1p3*p1p4**2*p2p3 + 8*mt**2*p1p4*p2p3**2 - 
     $      (32*p1p3*p1p4**2*p2p3**2)/s - 4*mt**4*s*p2p4 + 
     $      4*mt**2*s*p1p2*p2p4 - 12*mt**2*s*p1p3*p2p4 + 
     $      16*mt**2*p1p3*p1p4*p2p4 - 
     $      16*p1p2*p1p3*p1p4*p2p4 + 
     $      48*p1p3**2*p1p4*p2p4 - 
     $      8*mt**2*p1p3*p2p3*p2p4 - 
     $      8*mt**2*p1p4*p2p3*p2p4 + 
     $      (32*p1p3**2*p1p4*p2p3*p2p4)/s + 
     $      (32*p1p3*p1p4**2*p2p3*p2p4)/s + 
     $      8*mt**2*p1p3*p2p4**2 - (32*p1p3**2*p1p4*p2p4**2)/s)

       mat(2,41)=ct1*(-8*mt**6*s + 8*mt**4*s*p1p2 - 4*mt**4*s*p1p3 + 
     $      12*mt**2*s*p1p2*p1p3 + 64*mt**4*p1p3*p1p4 - 
     $      64*mt**2*p1p2*p1p3*p1p4 + 
     $      16*mt**2*p1p3**2*p1p4 - 48*p1p2*p1p3**2*p1p4 - 
     $      (128*mt**2*p1p3**2*p1p4**2)/s + 
     $      (128*p1p2*p1p3**2*p1p4**2)/s - 
     $      16*mt**4*p1p4*p2p3 - 24*mt**2*p1p3*p1p4*p2p3 + 
     $      (128*mt**2*p1p3*p1p4**2*p2p3)/s + 
     $      (96*p1p3**2*p1p4**2*p2p3)/s - 
     $      (256*p1p3**2*p1p4**3*p2p3)/s**2 - 8*mt**4*s*p2p4 + 
     $      16*mt**4*p1p3*p2p4 - 8*mt**2*s*p1p3*p2p4 - 
     $      8*mt**2*p1p3**2*p2p4 + 64*mt**2*p1p3*p1p4*p2p4 + 
     $      32*p1p3**2*p1p4*p2p4 - 
     $      (128*mt**2*p1p3**2*p1p4*p2p4)/s + 
     $      (32*p1p3**3*p1p4*p2p4)/s - 
     $      (128*p1p3**2*p1p4**2*p2p4)/s + 
     $      (256*p1p3**3*p1p4**2*p2p4)/s**2) + 
     $      ct3*(8*mt**4*s*p1p2 - 8*mt**2*s*p1p2**2 - 
     $      8*s*p1p2**2*p1p3 - 32*mt**2*p1p2*p1p3*p1p4 + 
     $      32*p1p2**2*p1p3*p1p4 - 16*mt**4*p1p4*p2p3 + 
     $      32*mt**2*p1p2*p1p4*p2p3 + 
     $      32*p1p2*p1p3*p1p4*p2p3 + 
     $      (64*mt**2*p1p3*p1p4**2*p2p3)/s - 
     $      (128*p1p2*p1p3*p1p4**2*p2p3)/s - 
     $      (32*mt**2*p1p4**2*p2p3**2)/s - 
     $      (32*p1p3*p1p4**2*p2p3**2)/s + 
     $      (128*p1p3*p1p4**3*p2p3**2)/s**2 - 8*mt**4*s*p2p4 - 
     $      16*mt**4*p1p3*p2p4 - 8*mt**2*s*p1p3*p2p4 + 
     $      32*p1p2*p1p3**2*p2p4 + 
     $      64*mt**2*p1p3*p1p4*p2p4 + 
     $      32*p1p3**2*p1p4*p2p4 + 
     $      (64*mt**2*p1p3**2*p1p4*p2p4)/s - 
     $      (128*p1p3**2*p1p4**2*p2p4)/s - 
     $      (64*p1p3**2*p1p4*p2p3*p2p4)/s + 
     $      (32*mt**2*p1p3**2*p2p4**2)/s - 
     $      (32*p1p3**3*p2p4**2)/s - 
     $      (128*p1p3**3*p1p4*p2p4**2)/s**2) + 
     $      ct2*(-8*mt**2*s*p1p2**2 + 8*s*p1p2**3 - 
     $      8*s*p1p2**2*p1p3 + 32*mt**2*p1p2*p1p4*p2p3 - 
     $      48*p1p2**2*p1p4*p2p3 + 
     $      32*p1p2*p1p3*p1p4*p2p3 - 
     $      (32*mt**2*p1p4**2*p2p3**2)/s + 
     $      (96*p1p2*p1p4**2*p2p3**2)/s - 
     $      (32*p1p3*p1p4**2*p2p3**2)/s - 
     $      (64*p1p4**3*p2p3**3)/s**2 - 4*mt**4*s*p2p4 + 
     $      12*mt**2*s*p1p2*p2p4 - 8*mt**2*s*p1p3*p2p4 + 
     $      32*mt**2*p1p2*p1p3*p2p4 - 
     $      16*p1p2**2*p1p3*p2p4 + 32*p1p2*p1p3**2*p2p4 + 
     $      16*mt**2*p1p3*p1p4*p2p4 - 
     $      48*p1p2*p1p3*p1p4*p2p4 + 
     $      32*p1p3**2*p1p4*p2p4 - 
     $      24*mt**2*p1p4*p2p3*p2p4 - 
     $      (64*mt**2*p1p3*p1p4*p2p3*p2p4)/s + 
     $      (64*p1p2*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (64*p1p3**2*p1p4*p2p3*p2p4)/s + 
     $      (96*p1p3*p1p4**2*p2p3*p2p4)/s - 
     $      (64*p1p3*p1p4**2*p2p3**2*p2p4)/s**2 - 
     $      8*mt**2*p1p3*p2p4**2 - (32*mt**2*p1p3**2*p2p4**2)/s - 
     $      (32*p1p2*p1p3**2*p2p4**2)/s - 
     $      (32*p1p3**3*p2p4**2)/s + 
     $      (32*p1p3**2*p1p4*p2p4**2)/s + 
     $      (64*p1p3**2*p1p4*p2p3*p2p4**2)/s**2 + 
     $      (64*p1p3**3*p2p4**3)/s**2)

       mat(3,41)=cu1*(-8*mt**6*s + 2*mt**4*s**2 + 8*mt**4*s*p1p2 - 
     $      2*mt**2*s**2*p1p2 + 8*mt**2*s*p1p2*p1p3 - 
     $      4*mt**4*s*p1p4 + 4*mt**2*s*p1p2*p1p4 + 
     $      64*mt**4*p1p3*p1p4 - 8*mt**2*s*p1p3*p1p4 - 
     $      64*mt**2*p1p2*p1p3*p1p4 + 
     $      8*s*p1p2*p1p3*p1p4 - 32*p1p2*p1p3**2*p1p4 + 
     $      16*mt**2*p1p3*p1p4**2 - 16*p1p2*p1p3*p1p4**2 - 
     $      (128*mt**2*p1p3**2*p1p4**2)/s + 
     $      (128*p1p2*p1p3**2*p1p4**2)/s - 
     $      16*mt**4*p1p4*p2p3 + 4*mt**2*s*p1p4*p2p3 - 
     $      16*mt**2*p1p3*p1p4*p2p3 - 8*mt**2*p1p4**2*p2p3 - 
     $      16*p1p3*p1p4**2*p2p3 + 
     $      (128*mt**2*p1p3*p1p4**2*p2p3)/s + 
     $      (64*p1p3**2*p1p4**2*p2p3)/s + 
     $      (32*p1p3*p1p4**3*p2p3)/s - 
     $      (256*p1p3**2*p1p4**3*p2p3)/s**2 - 8*mt**4*s*p2p4 + 
     $      16*mt**4*p1p3*p2p4 + 4*mt**2*s*p1p3*p2p4 - 
     $      16*mt**2*p1p3**2*p2p4 + 
     $      72*mt**2*p1p3*p1p4*p2p4 - 
     $      16*p1p3**2*p1p4*p2p4 - 
     $      (128*mt**2*p1p3**2*p1p4*p2p4)/s + 
     $      (64*p1p3**3*p1p4*p2p4)/s - 
     $      (160*p1p3**2*p1p4**2*p2p4)/s + 
     $      (256*p1p3**3*p1p4**2*p2p4)/s**2) + 
     $      cu3*(2*mt**4*s**2 + 8*mt**4*s*p1p2 - 2*mt**2*s**2*p1p2 - 
     $      8*mt**2*s*p1p2**2 + 8*mt**2*s*p1p2*p1p3 - 
     $      8*mt**2*s*p1p3*p1p4 - 32*mt**2*p1p2*p1p3*p1p4 + 
     $      8*s*p1p2*p1p3*p1p4 + 32*p1p2**2*p1p3*p1p4 - 
     $      32*p1p2*p1p3**2*p1p4 - 16*mt**4*p1p4*p2p3 + 
     $      4*mt**2*s*p1p4*p2p3 + 32*mt**2*p1p2*p1p4*p2p3 - 
     $      16*mt**2*p1p3*p1p4*p2p3 - 
     $      16*p1p3*p1p4**2*p2p3 + 
     $      (64*mt**2*p1p3*p1p4**2*p2p3)/s - 
     $      (128*p1p2*p1p3*p1p4**2*p2p3)/s + 
     $      (64*p1p3**2*p1p4**2*p2p3)/s - 
     $      (32*mt**2*p1p4**2*p2p3**2)/s + 
     $      (128*p1p3*p1p4**3*p2p3**2)/s**2 + 
     $      8*mt**2*s*p1p2*p2p4 - 16*mt**4*p1p3*p2p4 + 
     $      4*mt**2*s*p1p3*p2p4 - 16*mt**2*p1p3**2*p2p4 - 
     $      32*p1p2*p1p3*p1p4*p2p4 - 
     $      16*p1p3**2*p1p4*p2p4 + 
     $      (64*mt**2*p1p3**2*p1p4*p2p4)/s + 
     $      (64*p1p3**3*p1p4*p2p4)/s - 
     $      16*mt**2*p1p4*p2p3*p2p4 + 
     $      (64*p1p3*p1p4**2*p2p3*p2p4)/s - 
     $      16*mt**2*p1p3*p2p4**2 + 
     $      (32*mt**2*p1p3**2*p2p4**2)/s + 
     $      (64*p1p3**2*p1p4*p2p4**2)/s - 
     $      (128*p1p3**3*p1p4*p2p4**2)/s**2) + 
     $      cu2*(2*mt**4*s**2 - 2*mt**2*s**2*p1p2 -8*mt**2*s*p1p2**2+ 
     $      8*s*p1p2**3 - 8*s*p1p2**2*p1p3 - 
     $      8*mt**2*s*p1p3*p1p4 + 8*s*p1p2*p1p3*p1p4 - 
     $      4*mt**4*s*p2p3 + 4*mt**2*s*p1p2*p2p3 + 
     $      4*mt**2*s*p1p4*p2p3 + 32*mt**2*p1p2*p1p4*p2p3 - 
     $      48*p1p2**2*p1p4*p2p3 + 
     $      16*mt**2*p1p3*p1p4*p2p3 + 
     $      16*p1p2*p1p3*p1p4*p2p3 - 
     $      16*p1p3*p1p4**2*p2p3 - 8*mt**2*p1p4*p2p3**2 - 
     $      (32*mt**2*p1p4**2*p2p3**2)/s + 
     $      (96*p1p2*p1p4**2*p2p3**2)/s - 
     $      (64*p1p4**3*p2p3**3)/s**2 + 8*mt**2*s*p1p2*p2p4 + 
     $      4*mt**2*s*p1p3*p2p4 + 32*mt**2*p1p2*p1p3*p2p4 - 
     $      16*p1p2**2*p1p3*p2p4 + 32*p1p2*p1p3**2*p2p4 - 
     $      32*p1p2*p1p3*p1p4*p2p4 - 
     $      16*p1p3**2*p1p4*p2p4 + 
     $      8*mt**2*p1p3*p2p3*p2p4 - 
     $      16*mt**2*p1p4*p2p3*p2p4 - 
     $      (64*mt**2*p1p3*p1p4*p2p3*p2p4)/s + 
     $      (64*p1p2*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (96*p1p3**2*p1p4*p2p3*p2p4)/s + 
     $      (64*p1p3*p1p4**2*p2p3*p2p4)/s - 
     $      (64*p1p3*p1p4**2*p2p3**2*p2p4)/s**2 - 
     $      16*mt**2*p1p3*p2p4**2 - 
     $      (32*mt**2*p1p3**2*p2p4**2)/s - 
     $      (32*p1p2*p1p3**2*p2p4**2)/s - 
     $      (32*p1p3**3*p2p4**2)/s + 
     $      (64*p1p3**2*p1p4*p2p4**2)/s + 
     $      (64*p1p3**2*p1p4*p2p3*p2p4**2)/s**2 + 
     $      (64*p1p3**3*p2p4**3)/s**2)
       
       do i=1,3
          do j=42,45
             mat(i,j)=0d0
          enddo
       end do

       mat(1,46)=cs1*(4*mt**3*s*p1p4+8*mt**3*p1p3*p1p4-
     $      8*mt**3*p1p4**2+4*mt**3*s*p2p4-
     $      8*mt**3*p1p3*p2p4+8*mt**3*p1p4*p2p4-
     $      16*mt*p1p4*p2p3*p2p4-
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (32*mt*p1p4**2*p2p3*p2p4)/s-
     $      16*mt*p2p3*p2p4**2+
     $      (32*mt*p1p3*p2p3*p2p4**2)/s-
     $      (32*mt*p1p4*p2p3*p2p4**2)/s)+
     $      cs2*(4*mt**3*s*p1p4-8*mt**3*p1p4*p2p3+
     $      4*mt**3*s*p2p4+8*mt**3*p1p4*p2p4+
     $      8*mt**3*p2p3*p2p4-16*mt*p1p4*p2p3*p2p4+
     $      (32*mt*p1p4*p2p3**2*p2p4)/s-8*mt**3*p2p4**2-
     $      16*mt*p2p3*p2p4**2-
     $      (32*mt*p1p4*p2p3*p2p4**2)/s-
     $      (32*mt*p2p3**2*p2p4**2)/s+(32*mt*p2p3*p2p4**3)/s)

       mat(2,46)=ct1*(-4*mt**3*s*p1p2-4*mt*s*p1p2**2+
     $      16*mt*p1p2**2*p1p4+8*mt**3*p1p3*p1p4+
     $      8*mt**3*p1p4*p2p3+16*mt*p1p2*p1p4*p2p3-
     $      (64*mt*p1p2*p1p4**2*p2p3)/s-
     $      (16*mt*p1p4**2*p2p3**2)/s+
     $      (64*mt*p1p4**3*p2p3**2)/s**2+4*mt**3*s*p2p4-
     $      16*mt*p1p2**2*p2p4+16*mt*p1p2*p1p3*p2p4-
     $      (64*mt*p1p2*p1p3*p1p4*p2p4)/s+
     $      16*mt*p1p2*p2p3*p2p4+
     $      (64*mt*p1p2*p1p4*p2p3*p2p4)/s-
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (128*mt*p1p3*p1p4**2*p2p3*p2p4)/s**2-
     $      (32*mt*p1p4*p2p3**2*p2p4)/s-
     $      (64*mt*p1p4**2*p2p3**2*p2p4)/s**2+
     $      (64*mt*p1p2*p1p3*p2p4**2)/s-
     $      (16*mt*p1p3**2*p2p4**2)/s+
     $      (64*mt*p1p3**2*p1p4*p2p4**2)/s**2-
     $      16*mt*p2p3*p2p4**2-
     $      (128*mt*p1p3*p1p4*p2p3*p2p4**2)/s**2-
     $      (64*mt*p1p3**2*p2p4**3)/s**2)+
     $      ct3*(4*mt**5*s+4*mt**3*s*p1p2-16*mt**3*p1p2*p1p4-
     $      8*mt**3*p1p4*p2p3+(32*mt**3*p1p4**2*p2p3)/s+
     $      4*mt**3*s*p2p4+16*mt**3*p1p2*p2p4-
     $      8*mt**3*p1p3*p2p4+
     $      (32*mt**3*p1p3*p1p4*p2p4)/s-
     $      32*mt**3*p2p3*p2p4-16*mt*p1p2*p2p3*p2p4-
     $      (32*mt**3*p1p4*p2p3*p2p4)/s+
     $      (64*mt*p1p2*p1p4*p2p3*p2p4)/s+
     $      (32*mt*p1p4*p2p3**2*p2p4)/s-
     $      (128*mt*p1p4**2*p2p3**2*p2p4)/s**2-
     $      (32*mt**3*p1p3*p2p4**2)/s-16*mt*p2p3*p2p4**2-
     $      (64*mt*p1p2*p2p3*p2p4**2)/s+
     $      (32*mt*p1p3*p2p3*p2p4**2)/s-
     $      (128*mt*p1p3*p1p4*p2p3*p2p4**2)/s**2+
     $      (64*mt*p2p3**2*p2p4**2)/s+
     $      (128*mt*p1p4*p2p3**2*p2p4**2)/s**2+
     $      (128*mt*p1p3*p2p3*p2p4**3)/s**2)+
     $      ct2*(4*mt**5*s+4*mt**3*s*p1p2+16*mt**5*p1p4-
     $      8*mt**3*p1p4*p2p3-16*mt**5*p2p4+4*mt**3*s*p2p4-
     $      8*mt**3*p1p3*p2p4+8*mt**3*p1p4*p2p4-
     $      32*mt**3*p2p3*p2p4-16*mt*p1p2*p2p3*p2p4-
     $      (128*mt**3*p1p4*p2p3*p2p4)/s+
     $      (32*mt*p1p4*p2p3**2*p2p4)/s-8*mt**3*p2p4**2-
     $      16*mt*p2p3*p2p4**2+(128*mt**3*p2p3*p2p4**2)/s+
     $      (32*mt*p1p3*p2p3*p2p4**2)/s-
     $      (32*mt*p1p4*p2p3*p2p4**2)/s+
     $      (64*mt*p2p3**2*p2p4**2)/s+
     $      (256*mt*p1p4*p2p3**2*p2p4**2)/s**2+
     $      (32*mt*p2p3*p2p4**3)/s-(256*mt*p2p3**2*p2p4**3)/s**2)

       mat(3,46)=cu1*(-4*mt**3*s*p1p2-4*mt*s*p1p2**2-
     $      4*mt**3*s*p1p4+16*mt*p1p2**2*p1p4+
     $      8*mt**3*p1p4**2+8*mt**3*p1p4*p2p3+
     $      16*mt*p1p2*p1p4*p2p3-
     $      (64*mt*p1p2*p1p4**2*p2p3)/s-
     $      (16*mt*p1p4**2*p2p3**2)/s+
     $      (64*mt*p1p4**3*p2p3**2)/s**2-16*mt*p1p2**2*p2p4+
     $      8*mt**3*p1p3*p2p4+16*mt*p1p2*p1p3*p2p4-
     $      8*mt**3*p1p4*p2p4-
     $      (64*mt*p1p2*p1p3*p1p4*p2p4)/s+
     $      16*mt*p1p2*p2p3*p2p4+
     $      16*mt*p1p4*p2p3*p2p4+
     $      (64*mt*p1p2*p1p4*p2p3*p2p4)/s-
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (32*mt*p1p4**2*p2p3*p2p4)/s+
     $      (128*mt*p1p3*p1p4**2*p2p3*p2p4)/s**2-
     $      (32*mt*p1p4*p2p3**2*p2p4)/s-
     $      (64*mt*p1p4**2*p2p3**2*p2p4)/s**2+
     $      (64*mt*p1p2*p1p3*p2p4**2)/s-
     $      (16*mt*p1p3**2*p2p4**2)/s+
     $      (64*mt*p1p3**2*p1p4*p2p4**2)/s**2-
     $      (32*mt*p1p3*p2p3*p2p4**2)/s+
     $      (32*mt*p1p4*p2p3*p2p4**2)/s-
     $      (128*mt*p1p3*p1p4*p2p3*p2p4**2)/s**2-
     $      (64*mt*p1p3**2*p2p4**3)/s**2)+
     $      cu3*(-4*mt**3*s*p1p2-4*mt*s*p1p2**2-4*mt**3*s*p1p4-
     $      16*mt**3*p1p2*p1p4+8*mt**3*p1p4*p2p3+
     $      16*mt*p1p2*p1p4*p2p3+
     $      (32*mt**3*p1p4**2*p2p3)/s-
     $      (16*mt*p1p4**2*p2p3**2)/s+16*mt**3*p1p2*p2p4+
     $      8*mt**3*p1p3*p2p4+16*mt*p1p2*p1p3*p2p4+
     $      (32*mt**3*p1p3*p1p4*p2p4)/s+
     $      16*mt*p1p2*p2p3*p2p4+
     $      16*mt*p1p4*p2p3*p2p4-
     $      (32*mt**3*p1p4*p2p3*p2p4)/s+
     $      (64*mt*p1p2*p1p4*p2p3*p2p4)/s-
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (32*mt*p1p4*p2p3**2*p2p4)/s-
     $      (128*mt*p1p4**2*p2p3**2*p2p4)/s**2-
     $      (32*mt**3*p1p3*p2p4**2)/s-
     $      (16*mt*p1p3**2*p2p4**2)/s-
     $      (64*mt*p1p2*p2p3*p2p4**2)/s-
     $      (32*mt*p1p3*p2p3*p2p4**2)/s-
     $      (128*mt*p1p3*p1p4*p2p3*p2p4**2)/s**2+
     $      (128*mt*p1p4*p2p3**2*p2p4**2)/s**2+
     $      (128*mt*p1p3*p2p3*p2p4**3)/s**2)+
     $      cu2*(4*mt**5*s+4*mt**3*s*p1p2+16*mt**5*p1p4-
     $      4*mt**3*s*p1p4-16*mt**5*p2p4-8*mt**3*p1p3*p2p4-
     $      40*mt**3*p2p3*p2p4-16*mt*p1p2*p2p3*p2p4+
     $      16*mt*p1p4*p2p3*p2p4-
     $      (128*mt**3*p1p4*p2p3*p2p4)/s+
     $      (128*mt**3*p2p3*p2p4**2)/s+
     $      (32*mt*p1p3*p2p3*p2p4**2)/s+
     $      (96*mt*p2p3**2*p2p4**2)/s+
     $      (256*mt*p1p4*p2p3**2*p2p4**2)/s**2-
     $      (256*mt*p2p3**2*p2p4**3)/s**2)

       mat(1,47)=cs1*(4*mt*s*p1p2*p1p4+8*mt*p1p2*p1p3*p1p4-
     $      8*mt*p1p2*p1p4**2-8*mt*p1p4**2*p2p3-
     $      (16*mt*p1p3*p1p4**2*p2p3)/s+
     $      (16*mt*p1p4**3*p2p3)/s+4*mt*s*p1p2*p2p4-
     $      8*mt*p1p2*p1p3*p2p4+8*mt*p1p2*p1p4*p2p4-
     $      8*mt*p1p3*p1p4*p2p4-
     $      (16*mt*p1p3**2*p1p4*p2p4)/s+
     $      (16*mt*p1p3*p1p4**2*p2p4)/s-
     $      8*mt*p1p4*p2p3*p2p4+
     $      (16*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (16*mt*p1p4**2*p2p3*p2p4)/s-
     $      8*mt*p1p3*p2p4**2+(16*mt*p1p3**2*p2p4**2)/s-
     $      (16*mt*p1p3*p1p4*p2p4**2)/s)+
     $      cs2*(4*mt*s*p1p2*p1p4-8*mt*p1p2*p1p4*p2p3-
     $      8*mt*p1p4**2*p2p3+(16*mt*p1p4**2*p2p3**2)/s+
     $      4*mt*s*p1p2*p2p4+8*mt*p1p2*p1p4*p2p4-
     $      8*mt*p1p3*p1p4*p2p4+8*mt*p1p2*p2p3*p2p4-
     $      8*mt*p1p4*p2p3*p2p4+
     $      (16*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (16*mt*p1p4**2*p2p3*p2p4)/s-
     $      (16*mt*p1p4*p2p3**2*p2p4)/s-8*mt*p1p2*p2p4**2-
     $      8*mt*p1p3*p2p4**2-(16*mt*p1p3*p1p4*p2p4**2)/s-
     $      (16*mt*p1p3*p2p3*p2p4**2)/s+
     $      (16*mt*p1p4*p2p3*p2p4**2)/s+
     $      (16*mt*p1p3*p2p4**3)/s)

       mat(2,47)=ct1*(-4*mt**3*s*p1p2-4*mt*s*p1p2**2+
     $      16*mt**3*p1p2*p1p4+24*mt*p1p2*p1p3*p1p4-
     $      (64*mt*p1p2*p1p3*p1p4**2)/s+
     $      8*mt**3*p1p4*p2p3+16*mt*p1p2*p1p4*p2p3-
     $      (32*mt**3*p1p4**2*p2p3)/s-
     $      (48*mt*p1p3*p1p4**2*p2p3)/s+
     $      (128*mt*p1p3*p1p4**3*p2p3)/s**2-
     $      (16*mt*p1p4**2*p2p3**2)/s-16*mt**3*p1p2*p2p4+
     $      4*mt*s*p1p2*p2p4+8*mt**3*p1p3*p2p4+
     $      8*mt*p1p2*p1p3*p2p4-
     $      (32*mt**3*p1p3*p1p4*p2p4)/s+
     $      (64*mt*p1p2*p1p3*p1p4*p2p4)/s-
     $      (48*mt*p1p3**2*p1p4*p2p4)/s+
     $      (128*mt*p1p3**2*p1p4**2*p2p4)/s**2-
     $      8*mt*p1p4*p2p3*p2p4+
     $      (32*mt**3*p1p4*p2p3*p2p4)/s-
     $      (16*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (128*mt*p1p3*p1p4**2*p2p3*p2p4)/s**2-
     $      8*mt*p1p3*p2p4**2+(32*mt**3*p1p3*p2p4**2)/s-
     $      (128*mt*p1p3**2*p1p4*p2p4**2)/s**2)+
     $      ct3*(4*mt**5*s+4*mt**3*s*p1p2-16*mt**5*p1p4-
     $      16*mt**3*p1p3*p1p4+(64*mt**3*p1p3*p1p4**2)/s-
     $      8*mt**3*p1p4*p2p3+16*mt**5*p2p4+
     $      4*mt*s*p1p2*p2p4-8*mt**3*p1p3*p2p4-
     $      (64*mt**3*p1p3*p1p4*p2p4)/s-
     $      16*mt**3*p2p3*p2p4-16*mt*p1p2*p2p3*p2p4-
     $      8*mt*p1p4*p2p3*p2p4+
     $      (64*mt**3*p1p4*p2p3*p2p4)/s+
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (256*mt*p1p3*p1p4**2*p2p3*p2p4)/s**2+
     $      (32*mt*p1p4*p2p3**2*p2p4)/s-
     $      8*mt*p1p3*p2p4**2-(64*mt**3*p2p3*p2p4**2)/s+
     $      (32*mt*p1p3*p2p3*p2p4**2)/s+
     $      (256*mt*p1p3*p1p4*p2p3*p2p4**2)/s**2)+
     $      ct2*(4*mt**5*s+4*mt**3*s*p1p2+16*mt**3*p1p2*p1p4-
     $      16*mt**3*p1p3*p1p4-8*mt**3*p1p4*p2p3-
     $      (32*mt**3*p1p4**2*p2p3)/s-16*mt**3*p1p2*p2p4+
     $      4*mt*s*p1p2*p2p4-8*mt**3*p1p3*p2p4+
     $      8*mt*p1p2*p1p4*p2p4-
     $      (32*mt**3*p1p3*p1p4*p2p4)/s-
     $      16*mt**3*p2p3*p2p4-16*mt*p1p2*p2p3*p2p4-
     $      8*mt*p1p4*p2p3*p2p4+
     $      (32*mt**3*p1p4*p2p3*p2p4)/s-
     $      (64*mt*p1p2*p1p4*p2p3*p2p4)/s+
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (16*mt*p1p4**2*p2p3*p2p4)/s+
     $      (32*mt*p1p4*p2p3**2*p2p4)/s+
     $      (128*mt*p1p4**2*p2p3**2*p2p4)/s**2-
     $      8*mt*p1p2*p2p4**2-8*mt*p1p3*p2p4**2+
     $      (32*mt**3*p1p3*p2p4**2)/s-
     $      (16*mt*p1p3*p1p4*p2p4**2)/s+
     $      (64*mt*p1p2*p2p3*p2p4**2)/s+
     $      (32*mt*p1p3*p2p3*p2p4**2)/s+
     $      (16*mt*p1p4*p2p3*p2p4**2)/s+
     $      (128*mt*p1p3*p1p4*p2p3*p2p4**2)/s**2-
     $      (128*mt*p1p4*p2p3**2*p2p4**2)/s**2+
     $      (16*mt*p1p3*p2p4**3)/s-
     $      (128*mt*p1p3*p2p3*p2p4**3)/s**2)

       mat(3,47)=cu1*(-4*mt**3*s*p1p2-4*mt*s*p1p2**2+
     $      16*mt**3*p1p2*p1p4-4*mt*s*p1p2*p1p4+
     $      16*mt*p1p2*p1p3*p1p4+8*mt*p1p2*p1p4**2-
     $      (64*mt*p1p2*p1p3*p1p4**2)/s+
     $      8*mt**3*p1p4*p2p3+16*mt*p1p2*p1p4*p2p3+
     $      8*mt*p1p4**2*p2p3-(32*mt**3*p1p4**2*p2p3)/s-
     $      (32*mt*p1p3*p1p4**2*p2p3)/s-
     $      (16*mt*p1p4**3*p2p3)/s+
     $      (128*mt*p1p3*p1p4**3*p2p3)/s**2-
     $      (16*mt*p1p4**2*p2p3**2)/s-16*mt**3*p1p2*p2p4+
     $      8*mt**3*p1p3*p2p4+16*mt*p1p2*p1p3*p2p4-
     $      8*mt*p1p2*p1p4*p2p4+8*mt*p1p3*p1p4*p2p4-
     $      (32*mt**3*p1p3*p1p4*p2p4)/s+
     $      (64*mt*p1p2*p1p3*p1p4*p2p4)/s-
     $      (32*mt*p1p3**2*p1p4*p2p4)/s-
     $      (16*mt*p1p3*p1p4**2*p2p4)/s+
     $      (128*mt*p1p3**2*p1p4**2*p2p4)/s**2+
     $      (32*mt**3*p1p4*p2p3*p2p4)/s-
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (16*mt*p1p4**2*p2p3*p2p4)/s-
     $      (128*mt*p1p3*p1p4**2*p2p3*p2p4)/s**2+
     $      (32*mt**3*p1p3*p2p4**2)/s-
     $      (16*mt*p1p3**2*p2p4**2)/s+
     $      (16*mt*p1p3*p1p4*p2p4**2)/s-
     $      (128*mt*p1p3**2*p1p4*p2p4**2)/s**2)+
     $      cu3*(-4*mt**3*s*p1p2-4*mt*s*p1p2**2-
     $      4*mt*s*p1p2*p1p4-16*mt*p1p2**2*p1p4+
     $      16*mt*p1p2*p1p3*p1p4+8*mt**3*p1p4*p2p3+
     $      16*mt*p1p2*p1p4*p2p3+8*mt*p1p4**2*p2p3+
     $      (64*mt*p1p2*p1p4**2*p2p3)/s-
     $      (32*mt*p1p3*p1p4**2*p2p3)/s-
     $      (16*mt*p1p4**2*p2p3**2)/s-
     $      (64*mt*p1p4**3*p2p3**2)/s**2+16*mt*p1p2**2*p2p4+
     $      8*mt**3*p1p3*p2p4+16*mt*p1p2*p1p3*p2p4+
     $      8*mt*p1p3*p1p4*p2p4+
     $      (64*mt*p1p2*p1p3*p1p4*p2p4)/s-
     $      (32*mt*p1p3**2*p1p4*p2p4)/s-
     $      (64*mt*p1p2*p1p4*p2p3*p2p4)/s-
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (128*mt*p1p3*p1p4**2*p2p3*p2p4)/s**2+
     $      (64*mt*p1p4**2*p2p3**2*p2p4)/s**2-
     $      (64*mt*p1p2*p1p3*p2p4**2)/s-
     $      (16*mt*p1p3**2*p2p4**2)/s-
     $      (64*mt*p1p3**2*p1p4*p2p4**2)/s**2+
     $      (128*mt*p1p3*p1p4*p2p3*p2p4**2)/s**2+
     $      (64*mt*p1p3**2*p2p4**3)/s**2)+
     $      cu2*(4*mt**5*s+4*mt**3*s*p1p2+16*mt**3*p1p2*p1p4-
     $      4*mt*s*p1p2*p1p4-16*mt**3*p1p3*p1p4-
     $      8*mt**3*p1p4*p2p3+8*mt*p1p2*p1p4*p2p3+
     $      8*mt*p1p4**2*p2p3-(32*mt**3*p1p4**2*p2p3)/s-
     $      (16*mt*p1p4**2*p2p3**2)/s-16*mt**3*p1p2*p2p4-
     $      8*mt**3*p1p3*p2p4+8*mt*p1p3*p1p4*p2p4-
     $      (32*mt**3*p1p3*p1p4*p2p4)/s-
     $      16*mt**3*p2p3*p2p4-24*mt*p1p2*p2p3*p2p4+
     $      (32*mt**3*p1p4*p2p3*p2p4)/s-
     $      (64*mt*p1p2*p1p4*p2p3*p2p4)/s+
     $      (48*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (48*mt*p1p4*p2p3**2*p2p4)/s+
     $      (128*mt*p1p4**2*p2p3**2*p2p4)/s**2+
     $      (32*mt**3*p1p3*p2p4**2)/s+
     $      (64*mt*p1p2*p2p3*p2p4**2)/s+
     $      (48*mt*p1p3*p2p3*p2p4**2)/s+
     $      (128*mt*p1p3*p1p4*p2p3*p2p4**2)/s**2-
     $      (128*mt*p1p4*p2p3**2*p2p4**2)/s**2-
     $      (128*mt*p1p3*p2p3*p2p4**3)/s**2)

       do i=1,3
          mat(i,48)=0d0
       end do

       mat(1,49)=cs1*(4*mt*s*p1p2*p1p4+8*mt*p1p2*p1p3*p1p4-
     $      8*mt*p1p2*p1p4**2-8*mt*p1p4**2*p2p3-
     $      (16*mt*p1p3*p1p4**2*p2p3)/s+
     $      (16*mt*p1p4**3*p2p3)/s+4*mt*s*p1p2*p2p4-
     $      8*mt*p1p2*p1p3*p2p4+8*mt*p1p2*p1p4*p2p4-
     $      8*mt*p1p3*p1p4*p2p4-
     $      (16*mt*p1p3**2*p1p4*p2p4)/s+
     $      (16*mt*p1p3*p1p4**2*p2p4)/s-
     $      8*mt*p1p4*p2p3*p2p4+
     $      (16*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (16*mt*p1p4**2*p2p3*p2p4)/s-
     $      8*mt*p1p3*p2p4**2+(16*mt*p1p3**2*p2p4**2)/s-
     $      (16*mt*p1p3*p1p4*p2p4**2)/s)+
     $      cs2*(4*mt*s*p1p2*p1p4-8*mt*p1p2*p1p4*p2p3-
     $      8*mt*p1p4**2*p2p3+(16*mt*p1p4**2*p2p3**2)/s+
     $      4*mt*s*p1p2*p2p4+8*mt*p1p2*p1p4*p2p4-
     $      8*mt*p1p3*p1p4*p2p4+8*mt*p1p2*p2p3*p2p4-
     $      8*mt*p1p4*p2p3*p2p4+
     $      (16*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (16*mt*p1p4**2*p2p3*p2p4)/s-
     $      (16*mt*p1p4*p2p3**2*p2p4)/s-8*mt*p1p2*p2p4**2-
     $      8*mt*p1p3*p2p4**2-(16*mt*p1p3*p1p4*p2p4**2)/s-
     $      (16*mt*p1p3*p2p3*p2p4**2)/s+
     $      (16*mt*p1p4*p2p3*p2p4**2)/s+
     $      (16*mt*p1p3*p2p4**3)/s)

       mat(2,49)=ct1*(-4*mt**5*s-4*mt**3*s*p1p2+
     $      16*mt**3*p1p2*p1p4+16*mt**3*p1p3*p1p4+
     $      24*mt*p1p2*p1p3*p1p4-
     $      (64*mt*p1p2*p1p3*p1p4**2)/s+
     $      8*mt**3*p1p4*p2p3-(32*mt**3*p1p4**2*p2p3)/s-
     $      (48*mt*p1p3*p1p4**2*p2p3)/s+
     $      (128*mt*p1p3*p1p4**3*p2p3)/s**2-
     $      16*mt**3*p1p2*p2p4+4*mt*s*p1p2*p2p4+
     $      8*mt**3*p1p3*p2p4-8*mt*p1p2*p1p3*p2p4-
     $      (32*mt**3*p1p3*p1p4*p2p4)/s+
     $      (64*mt*p1p2*p1p3*p1p4*p2p4)/s-
     $      (48*mt*p1p3**2*p1p4*p2p4)/s+
     $      (128*mt*p1p3**2*p1p4**2*p2p4)/s**2+
     $      16*mt**3*p2p3*p2p4-8*mt*p1p4*p2p3*p2p4+
     $      (32*mt**3*p1p4*p2p3*p2p4)/s-
     $      (48*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (128*mt*p1p3*p1p4**2*p2p3*p2p4)/s**2-
     $      8*mt*p1p3*p2p4**2+(32*mt**3*p1p3*p2p4**2)/s+
     $      (16*mt*p1p3**2*p2p4**2)/s-
     $      (128*mt*p1p3**2*p1p4*p2p4**2)/s**2)+
     $      ct3*(4*mt**3*s*p1p2+4*mt*s*p1p2**2-
     $      16*mt*p1p2**2*p1p4-8*mt**3*p1p4*p2p3-
     $      16*mt*p1p2*p1p4*p2p3+
     $      (64*mt*p1p2*p1p4**2*p2p3)/s+
     $      (16*mt*p1p4**2*p2p3**2)/s-
     $      (64*mt*p1p4**3*p2p3**2)/s**2+4*mt*s*p1p2*p2p4+
     $      16*mt*p1p2**2*p2p4-8*mt**3*p1p3*p2p4-
     $      16*mt*p1p2*p1p3*p2p4+
     $      (64*mt*p1p2*p1p3*p1p4*p2p4)/s-
     $      16*mt*p1p2*p2p3*p2p4-8*mt*p1p4*p2p3*p2p4-
     $      (64*mt*p1p2*p1p4*p2p3*p2p4)/s+
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (128*mt*p1p3*p1p4**2*p2p3*p2p4)/s**2+
     $      (32*mt*p1p4*p2p3**2*p2p4)/s+
     $      (64*mt*p1p4**2*p2p3**2*p2p4)/s**2-
     $      8*mt*p1p3*p2p4**2-
     $      (64*mt*p1p2*p1p3*p2p4**2)/s+
     $      (16*mt*p1p3**2*p2p4**2)/s-
     $      (64*mt*p1p3**2*p1p4*p2p4**2)/s**2+
     $      (32*mt*p1p3*p2p3*p2p4**2)/s+
     $      (128*mt*p1p3*p1p4*p2p3*p2p4**2)/s**2+
     $      (64*mt*p1p3**2*p2p4**3)/s**2)+
     $      ct2*(4*mt**3*s*p1p2+4*mt*s*p1p2**2+
     $      16*mt**3*p1p2*p1p4-8*mt**3*p1p4*p2p3-
     $      16*mt*p1p2*p1p4*p2p3-
     $      (32*mt**3*p1p4**2*p2p3)/s+
     $      (16*mt*p1p4**2*p2p3**2)/s-16*mt**3*p1p2*p2p4+
     $      4*mt*s*p1p2*p2p4-8*mt**3*p1p3*p2p4-
     $      16*mt*p1p2*p1p3*p2p4+8*mt*p1p2*p1p4*p2p4-
     $      (32*mt**3*p1p3*p1p4*p2p4)/s-
     $      16*mt*p1p2*p2p3*p2p4-8*mt*p1p4*p2p3*p2p4+
     $      (32*mt**3*p1p4*p2p3*p2p4)/s-
     $      (64*mt*p1p2*p1p4*p2p3*p2p4)/s+
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (16*mt*p1p4**2*p2p3*p2p4)/s+
     $      (32*mt*p1p4*p2p3**2*p2p4)/s+
     $      (128*mt*p1p4**2*p2p3**2*p2p4)/s**2-
     $      8*mt*p1p2*p2p4**2-8*mt*p1p3*p2p4**2+
     $      (32*mt**3*p1p3*p2p4**2)/s+
     $      (16*mt*p1p3**2*p2p4**2)/s-
     $      (16*mt*p1p3*p1p4*p2p4**2)/s+
     $      (64*mt*p1p2*p2p3*p2p4**2)/s+
     $      (32*mt*p1p3*p2p3*p2p4**2)/s+
     $      (16*mt*p1p4*p2p3*p2p4**2)/s+
     $      (128*mt*p1p3*p1p4*p2p3*p2p4**2)/s**2-
     $      (128*mt*p1p4*p2p3**2*p2p4**2)/s**2+
     $      (16*mt*p1p3*p2p4**3)/s-
     $      (128*mt*p1p3*p2p3*p2p4**3)/s**2)

       mat(3,49)=cu1*(-4*mt**5*s-4*mt**3*s*p1p2+
     $      16*mt**3*p1p2*p1p4-4*mt*s*p1p2*p1p4+
     $      16*mt**3*p1p3*p1p4+16*mt*p1p2*p1p3*p1p4+
     $      8*mt*p1p2*p1p4**2-
     $      (64*mt*p1p2*p1p3*p1p4**2)/s+
     $      8*mt**3*p1p4*p2p3+8*mt*p1p4**2*p2p3-
     $      (32*mt**3*p1p4**2*p2p3)/s-
     $      (32*mt*p1p3*p1p4**2*p2p3)/s-
     $      (16*mt*p1p4**3*p2p3)/s+
     $      (128*mt*p1p3*p1p4**3*p2p3)/s**2-
     $      16*mt**3*p1p2*p2p4+8*mt**3*p1p3*p2p4-
     $      8*mt*p1p2*p1p4*p2p4+8*mt*p1p3*p1p4*p2p4-
     $      (32*mt**3*p1p3*p1p4*p2p4)/s+
     $      (64*mt*p1p2*p1p3*p1p4*p2p4)/s-
     $      (32*mt*p1p3**2*p1p4*p2p4)/s-
     $      (16*mt*p1p3*p1p4**2*p2p4)/s+
     $      (128*mt*p1p3**2*p1p4**2*p2p4)/s**2+
     $      16*mt**3*p2p3*p2p4+
     $      (32*mt**3*p1p4*p2p3*p2p4)/s-
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (16*mt*p1p4**2*p2p3*p2p4)/s-
     $      (128*mt*p1p3*p1p4**2*p2p3*p2p4)/s**2+
     $      (32*mt**3*p1p3*p2p4**2)/s+
     $      (16*mt*p1p3*p1p4*p2p4**2)/s-
     $      (128*mt*p1p3**2*p1p4*p2p4**2)/s**2)+
     $      cu3*(-4*mt**5*s-4*mt**3*s*p1p2-16*mt**5*p1p4-
     $      4*mt*s*p1p2*p1p4+16*mt**3*p1p3*p1p4+
     $      16*mt*p1p2*p1p3*p1p4+
     $      (64*mt**3*p1p3*p1p4**2)/s+8*mt**3*p1p4*p2p3+
     $      8*mt*p1p4**2*p2p3-
     $      (32*mt*p1p3*p1p4**2*p2p3)/s+16*mt**5*p2p4+
     $      8*mt**3*p1p3*p2p4+8*mt*p1p3*p1p4*p2p4-
     $      (64*mt**3*p1p3*p1p4*p2p4)/s-
     $      (32*mt*p1p3**2*p1p4*p2p4)/s+
     $      16*mt**3*p2p3*p2p4+
     $      (64*mt**3*p1p4*p2p3*p2p4)/s-
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (256*mt*p1p3*p1p4**2*p2p3*p2p4)/s**2-
     $      (64*mt**3*p2p3*p2p4**2)/s+
     $      (256*mt*p1p3*p1p4*p2p3*p2p4**2)/s**2)+
     $      cu2*(4*mt**3*s*p1p2+4*mt*s*p1p2**2+
     $      16*mt**3*p1p2*p1p4-4*mt*s*p1p2*p1p4-
     $      8*mt**3*p1p4*p2p3-8*mt*p1p2*p1p4*p2p3+
     $      8*mt*p1p4**2*p2p3-(32*mt**3*p1p4**2*p2p3)/s-
     $      16*mt**3*p1p2*p2p4-8*mt**3*p1p3*p2p4-
     $      16*mt*p1p2*p1p3*p2p4+8*mt*p1p3*p1p4*p2p4-
     $      (32*mt**3*p1p3*p1p4*p2p4)/s-
     $      24*mt*p1p2*p2p3*p2p4+
     $      (32*mt**3*p1p4*p2p3*p2p4)/s-
     $      (64*mt*p1p2*p1p4*p2p3*p2p4)/s+
     $      (16*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (48*mt*p1p4*p2p3**2*p2p4)/s+
     $      (128*mt*p1p4**2*p2p3**2*p2p4)/s**2+
     $      (32*mt**3*p1p3*p2p4**2)/s+
     $      (16*mt*p1p3**2*p2p4**2)/s+
     $      (64*mt*p1p2*p2p3*p2p4**2)/s+
     $      (48*mt*p1p3*p2p3*p2p4**2)/s+
     $      (128*mt*p1p3*p1p4*p2p3*p2p4**2)/s**2-
     $      (128*mt*p1p4*p2p3**2*p2p4**2)/s**2-
     $      (128*mt*p1p3*p2p3*p2p4**3)/s**2)

       mat(1,50)=cs1*(4*mt**3*s*p1p4+8*mt**3*p1p3*p1p4-
     $      8*mt**3*p1p4**2-16*mt*p1p3*p1p4**2-
     $      (32*mt*p1p3**2*p1p4**2)/s+(32*mt*p1p3*p1p4**3)/s+
     $      4*mt**3*s*p2p4-8*mt**3*p1p3*p2p4+
     $      8*mt**3*p1p4*p2p4-16*mt*p1p3*p1p4*p2p4+
     $      (32*mt*p1p3**2*p1p4*p2p4)/s-
     $      (32*mt*p1p3*p1p4**2*p2p4)/s)+
     $      cs2*(4*mt**3*s*p1p4-16*mt*p1p3*p1p4**2-
     $      8*mt**3*p1p4*p2p3+(32*mt*p1p3*p1p4**2*p2p3)/s+
     $      4*mt**3*s*p2p4+8*mt**3*p1p4*p2p4-
     $      16*mt*p1p3*p1p4*p2p4-
     $      (32*mt*p1p3*p1p4**2*p2p4)/s+8*mt**3*p2p3*p2p4-
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s-8*mt**3*p2p4**2+
     $      (32*mt*p1p3*p1p4*p2p4**2)/s)

       mat(2,50)=ct1*(-4*mt**5*s-4*mt**3*s*p1p2+16*mt**5*p1p4+
     $      40*mt**3*p1p3*p1p4+16*mt*p1p2*p1p3*p1p4-
     $      (128*mt**3*p1p3*p1p4**2)/s-
     $      (96*mt*p1p3**2*p1p4**2)/s+
     $      (256*mt*p1p3**2*p1p4**3)/s**2+8*mt**3*p1p4*p2p3-
     $      (32*mt*p1p3*p1p4**2*p2p3)/s-16*mt**5*p2p4+
     $      4*mt**3*s*p2p4-16*mt*p1p3*p1p4*p2p4+
     $      (128*mt**3*p1p3*p1p4*p2p4)/s-
     $      (256*mt*p1p3**2*p1p4**2*p2p4)/s**2)+
     $      ct3*(4*mt**3*s*p1p2+4*mt*s*p1p2**2-
     $      16*mt**3*p1p2*p1p4-16*mt*p1p2*p1p3*p1p4+
     $      (64*mt*p1p2*p1p3*p1p4**2)/s-
     $      8*mt**3*p1p4*p2p3-16*mt*p1p2*p1p4*p2p3+
     $      (32*mt**3*p1p4**2*p2p3)/s+
     $      (32*mt*p1p3*p1p4**2*p2p3)/s-
     $      (128*mt*p1p3*p1p4**3*p2p3)/s**2+
     $      (16*mt*p1p4**2*p2p3**2)/s+4*mt**3*s*p2p4+
     $      16*mt**3*p1p2*p2p4-8*mt**3*p1p3*p2p4-
     $      16*mt*p1p2*p1p3*p2p4-
     $      16*mt*p1p3*p1p4*p2p4+
     $      (32*mt**3*p1p3*p1p4*p2p4)/s-
     $      (64*mt*p1p2*p1p3*p1p4*p2p4)/s+
     $      (32*mt*p1p3**2*p1p4*p2p4)/s-
     $      (128*mt*p1p3**2*p1p4**2*p2p4)/s**2-
     $      (32*mt**3*p1p4*p2p3*p2p4)/s+
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (128*mt*p1p3*p1p4**2*p2p3*p2p4)/s**2-
     $      (32*mt**3*p1p3*p2p4**2)/s+
     $      (16*mt*p1p3**2*p2p4**2)/s+
     $      (128*mt*p1p3**2*p1p4*p2p4**2)/s**2)+
     $      ct2*(4*mt**3*s*p1p2+4*mt*s*p1p2**2+
     $      16*mt*p1p2**2*p1p4-16*mt*p1p2*p1p3*p1p4-
     $      8*mt**3*p1p4*p2p3-16*mt*p1p2*p1p4*p2p3-
     $      (64*mt*p1p2*p1p4**2*p2p3)/s+
     $      (32*mt*p1p3*p1p4**2*p2p3)/s+
     $      (16*mt*p1p4**2*p2p3**2)/s+
     $      (64*mt*p1p4**3*p2p3**2)/s**2+4*mt**3*s*p2p4-
     $      16*mt*p1p2**2*p2p4-8*mt**3*p1p3*p2p4-
     $      16*mt*p1p2*p1p3*p2p4+8*mt**3*p1p4*p2p4-
     $      16*mt*p1p3*p1p4*p2p4-
     $      (64*mt*p1p2*p1p3*p1p4*p2p4)/s+
     $      (32*mt*p1p3**2*p1p4*p2p4)/s-
     $      (32*mt*p1p3*p1p4**2*p2p4)/s+
     $      (64*mt*p1p2*p1p4*p2p3*p2p4)/s+
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (128*mt*p1p3*p1p4**2*p2p3*p2p4)/s**2-
     $      (64*mt*p1p4**2*p2p3**2*p2p4)/s**2-8*mt**3*p2p4**2+
     $      (64*mt*p1p2*p1p3*p2p4**2)/s+
     $      (16*mt*p1p3**2*p2p4**2)/s+
     $      (32*mt*p1p3*p1p4*p2p4**2)/s+
     $      (64*mt*p1p3**2*p1p4*p2p4**2)/s**2-
     $      (128*mt*p1p3*p1p4*p2p3*p2p4**2)/s**2-
     $      (64*mt*p1p3**2*p2p4**3)/s**2)

       mat(3,50)=cu1*(-4*mt**5*s-4*mt**3*s*p1p2+16*mt**5*p1p4-
     $      4*mt**3*s*p1p4+32*mt**3*p1p3*p1p4+
     $      16*mt*p1p2*p1p3*p1p4+8*mt**3*p1p4**2+
     $      16*mt*p1p3*p1p4**2-(128*mt**3*p1p3*p1p4**2)/s-
     $      (64*mt*p1p3**2*p1p4**2)/s-(32*mt*p1p3*p1p4**3)/s+
     $      (256*mt*p1p3**2*p1p4**3)/s**2+8*mt**3*p1p4*p2p3-
     $      (32*mt*p1p3*p1p4**2*p2p3)/s-16*mt**5*p2p4+
     $      8*mt**3*p1p3*p2p4-8*mt**3*p1p4*p2p4+
     $      (128*mt**3*p1p3*p1p4*p2p4)/s-
     $      (32*mt*p1p3**2*p1p4*p2p4)/s+
     $      (32*mt*p1p3*p1p4**2*p2p4)/s-
     $      (256*mt*p1p3**2*p1p4**2*p2p4)/s**2)+
     $      cu3*(-4*mt**5*s-4*mt**3*s*p1p2-4*mt**3*s*p1p4-
     $      16*mt**3*p1p2*p1p4+32*mt**3*p1p3*p1p4+
     $      16*mt*p1p2*p1p3*p1p4+16*mt*p1p3*p1p4**2+
     $      (64*mt*p1p2*p1p3*p1p4**2)/s-
     $      (64*mt*p1p3**2*p1p4**2)/s+8*mt**3*p1p4*p2p3+
     $      (32*mt**3*p1p4**2*p2p3)/s-
     $      (32*mt*p1p3*p1p4**2*p2p3)/s-
     $      (128*mt*p1p3*p1p4**3*p2p3)/s**2+
     $      16*mt**3*p1p2*p2p4+8*mt**3*p1p3*p2p4+
     $      (32*mt**3*p1p3*p1p4*p2p4)/s-
     $      (64*mt*p1p2*p1p3*p1p4*p2p4)/s-
     $      (32*mt*p1p3**2*p1p4*p2p4)/s-
     $      (128*mt*p1p3**2*p1p4**2*p2p4)/s**2-
     $      (32*mt**3*p1p4*p2p3*p2p4)/s+
     $      (128*mt*p1p3*p1p4**2*p2p3*p2p4)/s**2-
     $      (32*mt**3*p1p3*p2p4**2)/s+
     $      (128*mt*p1p3**2*p1p4*p2p4**2)/s**2)+
     $      cu2*(4*mt**3*s*p1p2+4*mt*s*p1p2**2-4*mt**3*s*p1p4+
     $      16*mt*p1p2**2*p1p4-16*mt*p1p2*p1p3*p1p4+
     $      16*mt*p1p3*p1p4**2-16*mt*p1p2*p1p4*p2p3-
     $      (64*mt*p1p2*p1p4**2*p2p3)/s+
     $      (16*mt*p1p4**2*p2p3**2)/s+
     $      (64*mt*p1p4**3*p2p3**2)/s**2-16*mt*p1p2**2*p2p4-
     $      8*mt**3*p1p3*p2p4-16*mt*p1p2*p1p3*p2p4-
     $      (64*mt*p1p2*p1p3*p1p4*p2p4)/s+
     $      (32*mt*p1p3**2*p1p4*p2p4)/s-8*mt**3*p2p3*p2p4+
     $      (64*mt*p1p2*p1p4*p2p3*p2p4)/s+
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (128*mt*p1p3*p1p4**2*p2p3*p2p4)/s**2-
     $      (64*mt*p1p4**2*p2p3**2*p2p4)/s**2+
     $      (64*mt*p1p2*p1p3*p2p4**2)/s+
     $      (16*mt*p1p3**2*p2p4**2)/s+
     $      (64*mt*p1p3**2*p1p4*p2p4**2)/s**2-
     $      (128*mt*p1p3*p1p4*p2p3*p2p4**2)/s**2-
     $      (64*mt*p1p3**2*p2p4**3)/s**2)

       do i=1,3
          do j=51,54
             mat(i,j)=0d0
          enddo
       end do

       mat(1,55)=cs1*(-4*mt**3*s*p1p3+8*mt**3*p1p3**2-
     $      8*mt**3*p1p3*p1p4-4*mt**3*s*p2p3-
     $      8*mt**3*p1p3*p2p3+8*mt**3*p1p4*p2p3+
     $      16*mt*p1p3*p2p3*p2p4-
     $      (32*mt*p1p3**2*p2p3*p2p4)/s+
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      16*mt*p2p3**2*p2p4+
     $      (32*mt*p1p3*p2p3**2*p2p4)/s-
     $      (32*mt*p1p4*p2p3**2*p2p4)/s)+
     $      cs2*(-4*mt**3*s*p1p3-4*mt**3*s*p2p3-
     $      8*mt**3*p1p3*p2p3+8*mt**3*p2p3**2+
     $      8*mt**3*p1p3*p2p4-8*mt**3*p2p3*p2p4+
     $      16*mt*p1p3*p2p3*p2p4+16*mt*p2p3**2*p2p4+
     $      (32*mt*p1p3*p2p3**2*p2p4)/s-
     $      (32*mt*p2p3**3*p2p4)/s-
     $      (32*mt*p1p3*p2p3*p2p4**2)/s+
     $      (32*mt*p2p3**2*p2p4**2)/s)

       mat(2,55)=ct1*(-4*mt**3*s*p1p2-4*mt*s*p1p2**2-
     $      4*mt**3*s*p1p3+16*mt*p1p2**2*p1p3+
     $      8*mt**3*p1p3**2-16*mt*p1p2**2*p2p3-
     $      8*mt**3*p1p3*p2p3+8*mt**3*p1p4*p2p3+
     $      16*mt*p1p2*p1p4*p2p3-
     $      (64*mt*p1p2*p1p3*p1p4*p2p3)/s+
     $      (64*mt*p1p2*p1p4*p2p3**2)/s-
     $      (16*mt*p1p4**2*p2p3**2)/s+
     $      (64*mt*p1p3*p1p4**2*p2p3**2)/s**2-
     $      (64*mt*p1p4**2*p2p3**3)/s**2+8*mt**3*p1p3*p2p4+
     $      16*mt*p1p2*p1p3*p2p4-
     $      (64*mt*p1p2*p1p3**2*p2p4)/s+
     $      16*mt*p1p2*p2p3*p2p4+
     $      16*mt*p1p3*p2p3*p2p4+
     $      (64*mt*p1p2*p1p3*p2p3*p2p4)/s-
     $      (32*mt*p1p3**2*p2p3*p2p4)/s-
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (128*mt*p1p3**2*p1p4*p2p3*p2p4)/s**2+
     $      (32*mt*p1p3*p2p3**2*p2p4)/s-
     $      (32*mt*p1p4*p2p3**2*p2p4)/s-
     $      (128*mt*p1p3*p1p4*p2p3**2*p2p4)/s**2-
     $      (16*mt*p1p3**2*p2p4**2)/s+
     $      (64*mt*p1p3**3*p2p4**2)/s**2-
     $      (32*mt*p1p3*p2p3*p2p4**2)/s-
     $      (64*mt*p1p3**2*p2p3*p2p4**2)/s**2)+
     $      ct3*(-4*mt**3*s*p1p2-4*mt*s*p1p2**2-4*mt**3*s*p1p3-
     $      16*mt**3*p1p2*p1p3+16*mt**3*p1p2*p2p3+
     $      8*mt**3*p1p4*p2p3+16*mt*p1p2*p1p4*p2p3+
     $      (32*mt**3*p1p3*p1p4*p2p3)/s-
     $      (32*mt**3*p1p4*p2p3**2)/s-
     $      (16*mt*p1p4**2*p2p3**2)/s+8*mt**3*p1p3*p2p4+
     $      16*mt*p1p2*p1p3*p2p4+
     $      (32*mt**3*p1p3**2*p2p4)/s+
     $      16*mt*p1p2*p2p3*p2p4+
     $      16*mt*p1p3*p2p3*p2p4-
     $      (32*mt**3*p1p3*p2p3*p2p4)/s+
     $      (64*mt*p1p2*p1p3*p2p3*p2p4)/s-
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (64*mt*p1p2*p2p3**2*p2p4)/s-
     $      (32*mt*p1p4*p2p3**2*p2p4)/s-
     $      (128*mt*p1p3*p1p4*p2p3**2*p2p4)/s**2+
     $      (128*mt*p1p4*p2p3**3*p2p4)/s**2-
     $      (16*mt*p1p3**2*p2p4**2)/s-
     $      (32*mt*p1p3*p2p3*p2p4**2)/s-
     $      (128*mt*p1p3**2*p2p3*p2p4**2)/s**2+
     $      (128*mt*p1p3*p2p3**2*p2p4**2)/s**2)+
     $      ct2*(4*mt**5*s+4*mt**3*s*p1p2+16*mt**5*p1p3-
     $      4*mt**3*s*p1p3-16*mt**5*p2p3-8*mt**3*p1p4*p2p3-
     $      40*mt**3*p2p3*p2p4-16*mt*p1p2*p2p3*p2p4+
     $      16*mt*p1p3*p2p3*p2p4-
     $      (128*mt**3*p1p3*p2p3*p2p4)/s+
     $      (128*mt**3*p2p3**2*p2p4)/s+
     $      (32*mt*p1p4*p2p3**2*p2p4)/s+
     $      (96*mt*p2p3**2*p2p4**2)/s+
     $      (256*mt*p1p3*p2p3**2*p2p4**2)/s**2-
     $      (256*mt*p2p3**3*p2p4**2)/s**2)

       mat(3,55)=cu1*(-4*mt**3*s*p1p2-4*mt*s*p1p2**2+
     $      16*mt*p1p2**2*p1p3+8*mt**3*p1p3*p1p4+
     $      4*mt**3*s*p2p3-16*mt*p1p2**2*p2p3+
     $      16*mt*p1p2*p1p4*p2p3-
     $      (64*mt*p1p2*p1p3*p1p4*p2p3)/s+
     $      (64*mt*p1p2*p1p4*p2p3**2)/s-
     $      (16*mt*p1p4**2*p2p3**2)/s+
     $      (64*mt*p1p3*p1p4**2*p2p3**2)/s**2-
     $      (64*mt*p1p4**2*p2p3**3)/s**2+8*mt**3*p1p3*p2p4+
     $      16*mt*p1p2*p1p3*p2p4-
     $      (64*mt*p1p2*p1p3**2*p2p4)/s+
     $      16*mt*p1p2*p2p3*p2p4+
     $      (64*mt*p1p2*p1p3*p2p3*p2p4)/s-
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (128*mt*p1p3**2*p1p4*p2p3*p2p4)/s**2-
     $      16*mt*p2p3**2*p2p4-
     $      (128*mt*p1p3*p1p4*p2p3**2*p2p4)/s**2-
     $      (16*mt*p1p3**2*p2p4**2)/s+
     $      (64*mt*p1p3**3*p2p4**2)/s**2-
     $      (32*mt*p1p3*p2p3*p2p4**2)/s-
     $      (64*mt*p1p3**2*p2p3*p2p4**2)/s**2)+
     $      cu3*(4*mt**5*s+4*mt**3*s*p1p2-16*mt**3*p1p2*p1p3+
     $      4*mt**3*s*p2p3+16*mt**3*p1p2*p2p3-
     $      8*mt**3*p1p4*p2p3+
     $      (32*mt**3*p1p3*p1p4*p2p3)/s-
     $      (32*mt**3*p1p4*p2p3**2)/s-8*mt**3*p1p3*p2p4+
     $      (32*mt**3*p1p3**2*p2p4)/s-32*mt**3*p2p3*p2p4-
     $      16*mt*p1p2*p2p3*p2p4-
     $      (32*mt**3*p1p3*p2p3*p2p4)/s+
     $      (64*mt*p1p2*p1p3*p2p3*p2p4)/s-
     $      16*mt*p2p3**2*p2p4-
     $      (64*mt*p1p2*p2p3**2*p2p4)/s+
     $      (32*mt*p1p4*p2p3**2*p2p4)/s-
     $      (128*mt*p1p3*p1p4*p2p3**2*p2p4)/s**2+
     $      (128*mt*p1p4*p2p3**3*p2p4)/s**2+
     $      (32*mt*p1p3*p2p3*p2p4**2)/s-
     $      (128*mt*p1p3**2*p2p3*p2p4**2)/s**2+
     $      (64*mt*p2p3**2*p2p4**2)/s+
     $      (128*mt*p1p3*p2p3**2*p2p4**2)/s**2)+
     $      cu2*(4*mt**5*s+4*mt**3*s*p1p2+16*mt**5*p1p3-
     $      16*mt**5*p2p3+4*mt**3*s*p2p3+8*mt**3*p1p3*p2p3-
     $      8*mt**3*p1p4*p2p3-8*mt**3*p2p3**2-
     $      8*mt**3*p1p3*p2p4-32*mt**3*p2p3*p2p4-
     $      16*mt*p1p2*p2p3*p2p4-
     $      (128*mt**3*p1p3*p2p3*p2p4)/s-
     $      16*mt*p2p3**2*p2p4+(128*mt**3*p2p3**2*p2p4)/s-
     $      (32*mt*p1p3*p2p3**2*p2p4)/s+
     $      (32*mt*p1p4*p2p3**2*p2p4)/s+
     $      (32*mt*p2p3**3*p2p4)/s+
     $      (32*mt*p1p3*p2p3*p2p4**2)/s+
     $      (64*mt*p2p3**2*p2p4**2)/s+
     $      (256*mt*p1p3*p2p3**2*p2p4**2)/s**2-
     $      (256*mt*p2p3**3*p2p4**2)/s**2)

       mat(1,56)=cs1*(-4*mt*s*p1p2*p1p3+8*mt*p1p2*p1p3**2-
     $      8*mt*p1p2*p1p3*p1p4-4*mt*s*p1p2*p2p3-
     $      8*mt*p1p2*p1p3*p2p3+8*mt*p1p2*p1p4*p2p3+
     $      8*mt*p1p3*p1p4*p2p3-
     $      (16*mt*p1p3**2*p1p4*p2p3)/s+
     $      (16*mt*p1p3*p1p4**2*p2p3)/s+
     $      8*mt*p1p4*p2p3**2+
     $      (16*mt*p1p3*p1p4*p2p3**2)/s-
     $      (16*mt*p1p4**2*p2p3**2)/s+8*mt*p1p3**2*p2p4-
     $      (16*mt*p1p3**3*p2p4)/s+
     $      (16*mt*p1p3**2*p1p4*p2p4)/s+
     $      8*mt*p1p3*p2p3*p2p4+
     $      (16*mt*p1p3**2*p2p3*p2p4)/s-
     $      (16*mt*p1p3*p1p4*p2p3*p2p4)/s)+
     $      cs2*(-4*mt*s*p1p2*p1p3-4*mt*s*p1p2*p2p3-
     $      8*mt*p1p2*p1p3*p2p3+8*mt*p1p3*p1p4*p2p3+
     $      8*mt*p1p2*p2p3**2+8*mt*p1p4*p2p3**2+
     $      (16*mt*p1p3*p1p4*p2p3**2)/s-
     $      (16*mt*p1p4*p2p3**3)/s+8*mt*p1p2*p1p3*p2p4+
     $      8*mt*p1p3**2*p2p4-8*mt*p1p2*p2p3*p2p4+
     $      8*mt*p1p3*p2p3*p2p4+
     $      (16*mt*p1p3**2*p2p3*p2p4)/s-
     $      (16*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (16*mt*p1p3*p2p3**2*p2p4)/s+
     $      (16*mt*p1p4*p2p3**2*p2p4)/s-
     $      (16*mt*p1p3**2*p2p4**2)/s+
     $      (16*mt*p1p3*p2p3*p2p4**2)/s)

       mat(2,56)=ct1*(-4*mt**5*s-4*mt**3*s*p1p2+
     $      16*mt**3*p1p2*p1p3-4*mt*s*p1p2*p1p3+
     $      8*mt*p1p2*p1p3**2+16*mt**3*p1p3*p1p4+
     $      16*mt*p1p2*p1p3*p1p4-
     $      (64*mt*p1p2*p1p3**2*p1p4)/s-
     $      16*mt**3*p1p2*p2p3-8*mt*p1p2*p1p3*p2p3+
     $      8*mt**3*p1p4*p2p3+8*mt*p1p3*p1p4*p2p3-
     $      (32*mt**3*p1p3*p1p4*p2p3)/s+
     $      (64*mt*p1p2*p1p3*p1p4*p2p3)/s-
     $      (16*mt*p1p3**2*p1p4*p2p3)/s-
     $      (32*mt*p1p3*p1p4**2*p2p3)/s+
     $      (128*mt*p1p3**2*p1p4**2*p2p3)/s**2+
     $      (32*mt**3*p1p4*p2p3**2)/s+
     $      (16*mt*p1p3*p1p4*p2p3**2)/s-
     $      (128*mt*p1p3*p1p4**2*p2p3**2)/s**2+
     $      8*mt**3*p1p3*p2p4+8*mt*p1p3**2*p2p4-
     $      (32*mt**3*p1p3**2*p2p4)/s-(16*mt*p1p3**3*p2p4)/s-
     $      (32*mt*p1p3**2*p1p4*p2p4)/s+
     $      (128*mt*p1p3**3*p1p4*p2p4)/s**2+
     $      16*mt**3*p2p3*p2p4+
     $      (32*mt**3*p1p3*p2p3*p2p4)/s+
     $      (16*mt*p1p3**2*p2p3*p2p4)/s-
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (128*mt*p1p3**2*p1p4*p2p3*p2p4)/s**2)+
     $      ct3*(-4*mt**5*s-4*mt**3*s*p1p2-16*mt**5*p1p3-
     $      4*mt*s*p1p2*p1p3+16*mt**3*p1p3*p1p4+
     $      16*mt*p1p2*p1p3*p1p4+
     $      (64*mt**3*p1p3**2*p1p4)/s+16*mt**5*p2p3+
     $      8*mt**3*p1p4*p2p3+8*mt*p1p3*p1p4*p2p3-
     $      (64*mt**3*p1p3*p1p4*p2p3)/s-
     $      (32*mt*p1p3*p1p4**2*p2p3)/s+
     $      8*mt**3*p1p3*p2p4+8*mt*p1p3**2*p2p4-
     $      (32*mt*p1p3**2*p1p4*p2p4)/s+
     $      16*mt**3*p2p3*p2p4+
     $      (64*mt**3*p1p3*p2p3*p2p4)/s-
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (256*mt*p1p3**2*p1p4*p2p3*p2p4)/s**2-
     $      (64*mt**3*p2p3**2*p2p4)/s+
     $      (256*mt*p1p3*p1p4*p2p3**2*p2p4)/s**2)+
     $      ct2*(4*mt**3*s*p1p2+4*mt*s*p1p2**2+
     $      16*mt**3*p1p2*p1p3-4*mt*s*p1p2*p1p3-
     $      16*mt**3*p1p2*p2p3-8*mt**3*p1p4*p2p3-
     $      16*mt*p1p2*p1p4*p2p3+8*mt*p1p3*p1p4*p2p3-
     $      (32*mt**3*p1p3*p1p4*p2p3)/s+
     $      (32*mt**3*p1p4*p2p3**2)/s+
     $      (16*mt*p1p4**2*p2p3**2)/s-8*mt**3*p1p3*p2p4-
     $      8*mt*p1p2*p1p3*p2p4+8*mt*p1p3**2*p2p4-
     $      (32*mt**3*p1p3**2*p2p4)/s-
     $      24*mt*p1p2*p2p3*p2p4+
     $      (32*mt**3*p1p3*p2p3*p2p4)/s-
     $      (64*mt*p1p2*p1p3*p2p3*p2p4)/s+
     $      (16*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (64*mt*p1p2*p2p3**2*p2p4)/s+
     $      (48*mt*p1p4*p2p3**2*p2p4)/s+
     $      (128*mt*p1p3*p1p4*p2p3**2*p2p4)/s**2-
     $      (128*mt*p1p4*p2p3**3*p2p4)/s**2+
     $      (48*mt*p1p3*p2p3*p2p4**2)/s+
     $      (128*mt*p1p3**2*p2p3*p2p4**2)/s**2-
     $      (128*mt*p1p3*p2p3**2*p2p4**2)/s**2)

       mat(3,56)=cu1*(-4*mt**5*s-4*mt**3*s*p1p2+
     $      16*mt**3*p1p2*p1p3+16*mt**3*p1p3*p1p4+
     $      24*mt*p1p2*p1p3*p1p4-
     $      (64*mt*p1p2*p1p3**2*p1p4)/s-
     $      16*mt**3*p1p2*p2p3+4*mt*s*p1p2*p2p3+
     $      8*mt**3*p1p4*p2p3-8*mt*p1p2*p1p4*p2p3-
     $      (32*mt**3*p1p3*p1p4*p2p3)/s+
     $      (64*mt*p1p2*p1p3*p1p4*p2p3)/s-
     $      (48*mt*p1p3*p1p4**2*p2p3)/s+
     $      (128*mt*p1p3**2*p1p4**2*p2p3)/s**2-
     $      8*mt*p1p4*p2p3**2+(32*mt**3*p1p4*p2p3**2)/s+
     $      (16*mt*p1p4**2*p2p3**2)/s-
     $      (128*mt*p1p3*p1p4**2*p2p3**2)/s**2+
     $      8*mt**3*p1p3*p2p4-(32*mt**3*p1p3**2*p2p4)/s-
     $      (48*mt*p1p3**2*p1p4*p2p4)/s+
     $      (128*mt*p1p3**3*p1p4*p2p4)/s**2+
     $      16*mt**3*p2p3*p2p4-8*mt*p1p3*p2p3*p2p4+
     $      (32*mt**3*p1p3*p2p3*p2p4)/s-
     $      (48*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (128*mt*p1p3**2*p1p4*p2p3*p2p4)/s**2)+
     $      cu3*(4*mt**3*s*p1p2+4*mt*s*p1p2**2-
     $      16*mt*p1p2**2*p1p3+4*mt*s*p1p2*p2p3+
     $      16*mt*p1p2**2*p2p3-8*mt**3*p1p4*p2p3-
     $      16*mt*p1p2*p1p4*p2p3+
     $      (64*mt*p1p2*p1p3*p1p4*p2p3)/s-
     $      8*mt*p1p4*p2p3**2-
     $      (64*mt*p1p2*p1p4*p2p3**2)/s+
     $      (16*mt*p1p4**2*p2p3**2)/s-
     $      (64*mt*p1p3*p1p4**2*p2p3**2)/s**2+
     $      (64*mt*p1p4**2*p2p3**3)/s**2-8*mt**3*p1p3*p2p4-
     $      16*mt*p1p2*p1p3*p2p4+
     $      (64*mt*p1p2*p1p3**2*p2p4)/s-
     $      16*mt*p1p2*p2p3*p2p4-8*mt*p1p3*p2p3*p2p4-
     $      (64*mt*p1p2*p1p3*p2p3*p2p4)/s+
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (128*mt*p1p3**2*p1p4*p2p3*p2p4)/s**2+
     $      (32*mt*p1p4*p2p3**2*p2p4)/s+
     $      (128*mt*p1p3*p1p4*p2p3**2*p2p4)/s**2+
     $      (16*mt*p1p3**2*p2p4**2)/s-
     $      (64*mt*p1p3**3*p2p4**2)/s**2+
     $      (32*mt*p1p3*p2p3*p2p4**2)/s+
     $      (64*mt*p1p3**2*p2p3*p2p4**2)/s**2)+
     $      cu2*(4*mt**3*s*p1p2+4*mt*s*p1p2**2+
     $      16*mt**3*p1p2*p1p3-16*mt**3*p1p2*p2p3+
     $      4*mt*s*p1p2*p2p3+8*mt*p1p2*p1p3*p2p3-
     $      8*mt**3*p1p4*p2p3-16*mt*p1p2*p1p4*p2p3-
     $      (32*mt**3*p1p3*p1p4*p2p3)/s-8*mt*p1p2*p2p3**2-
     $      8*mt*p1p4*p2p3**2+(32*mt**3*p1p4*p2p3**2)/s-
     $      (16*mt*p1p3*p1p4*p2p3**2)/s+
     $      (16*mt*p1p4**2*p2p3**2)/s+(16*mt*p1p4*p2p3**3)/s-
     $      8*mt**3*p1p3*p2p4-16*mt*p1p2*p1p3*p2p4-
     $      (32*mt**3*p1p3**2*p2p4)/s-
     $      16*mt*p1p2*p2p3*p2p4-8*mt*p1p3*p2p3*p2p4+
     $      (32*mt**3*p1p3*p2p3*p2p4)/s-
     $      (64*mt*p1p2*p1p3*p2p3*p2p4)/s-
     $      (16*mt*p1p3**2*p2p3*p2p4)/s+
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (64*mt*p1p2*p2p3**2*p2p4)/s+
     $      (16*mt*p1p3*p2p3**2*p2p4)/s+
     $      (32*mt*p1p4*p2p3**2*p2p4)/s+
     $      (128*mt*p1p3*p1p4*p2p3**2*p2p4)/s**2-
     $      (128*mt*p1p4*p2p3**3*p2p4)/s**2+
     $      (16*mt*p1p3**2*p2p4**2)/s+
     $      (32*mt*p1p3*p2p3*p2p4**2)/s+
     $      (128*mt*p1p3**2*p2p3*p2p4**2)/s**2-
     $      (128*mt*p1p3*p2p3**2*p2p4**2)/s**2)

       do i=1,3
          mat(i,57)=0d0
       end do

       mat(1,58)=cs1*(-4*mt*s*p1p2*p1p3+8*mt*p1p2*p1p3**2-
     $      8*mt*p1p2*p1p3*p1p4-4*mt*s*p1p2*p2p3-
     $      8*mt*p1p2*p1p3*p2p3+8*mt*p1p2*p1p4*p2p3+
     $      8*mt*p1p3*p1p4*p2p3-
     $      (16*mt*p1p3**2*p1p4*p2p3)/s+
     $      (16*mt*p1p3*p1p4**2*p2p3)/s+
     $      8*mt*p1p4*p2p3**2+
     $      (16*mt*p1p3*p1p4*p2p3**2)/s-
     $      (16*mt*p1p4**2*p2p3**2)/s+8*mt*p1p3**2*p2p4-
     $      (16*mt*p1p3**3*p2p4)/s+
     $      (16*mt*p1p3**2*p1p4*p2p4)/s+
     $      8*mt*p1p3*p2p3*p2p4+
     $      (16*mt*p1p3**2*p2p3*p2p4)/s-
     $      (16*mt*p1p3*p1p4*p2p3*p2p4)/s)+
     $      cs2*(-4*mt*s*p1p2*p1p3-4*mt*s*p1p2*p2p3-
     $      8*mt*p1p2*p1p3*p2p3+8*mt*p1p3*p1p4*p2p3+
     $      8*mt*p1p2*p2p3**2+8*mt*p1p4*p2p3**2+
     $      (16*mt*p1p3*p1p4*p2p3**2)/s-
     $      (16*mt*p1p4*p2p3**3)/s+8*mt*p1p2*p1p3*p2p4+
     $      8*mt*p1p3**2*p2p4-8*mt*p1p2*p2p3*p2p4+
     $      8*mt*p1p3*p2p3*p2p4+
     $      (16*mt*p1p3**2*p2p3*p2p4)/s-
     $      (16*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (16*mt*p1p3*p2p3**2*p2p4)/s+
     $      (16*mt*p1p4*p2p3**2*p2p4)/s-
     $      (16*mt*p1p3**2*p2p4**2)/s+
     $      (16*mt*p1p3*p2p3*p2p4**2)/s)

       mat(2,58)=ct1*(-4*mt**3*s*p1p2-4*mt*s*p1p2**2+
     $      16*mt**3*p1p2*p1p3-4*mt*s*p1p2*p1p3+
     $      8*mt*p1p2*p1p3**2+16*mt*p1p2*p1p3*p1p4-
     $      (64*mt*p1p2*p1p3**2*p1p4)/s-
     $      16*mt**3*p1p2*p2p3-8*mt*p1p2*p1p3*p2p3+
     $      8*mt**3*p1p4*p2p3+16*mt*p1p2*p1p4*p2p3+
     $      8*mt*p1p3*p1p4*p2p3-
     $      (32*mt**3*p1p3*p1p4*p2p3)/s+
     $      (64*mt*p1p2*p1p3*p1p4*p2p3)/s-
     $      (16*mt*p1p3**2*p1p4*p2p3)/s-
     $      (32*mt*p1p3*p1p4**2*p2p3)/s+
     $      (128*mt*p1p3**2*p1p4**2*p2p3)/s**2+
     $      (32*mt**3*p1p4*p2p3**2)/s+
     $      (16*mt*p1p3*p1p4*p2p3**2)/s-
     $      (16*mt*p1p4**2*p2p3**2)/s-
     $      (128*mt*p1p3*p1p4**2*p2p3**2)/s**2+
     $      8*mt**3*p1p3*p2p4+16*mt*p1p2*p1p3*p2p4+
     $      8*mt*p1p3**2*p2p4-(32*mt**3*p1p3**2*p2p4)/s-
     $      (16*mt*p1p3**3*p2p4)/s-
     $      (32*mt*p1p3**2*p1p4*p2p4)/s+
     $      (128*mt*p1p3**3*p1p4*p2p4)/s**2+
     $      (32*mt**3*p1p3*p2p3*p2p4)/s+
     $      (16*mt*p1p3**2*p2p3*p2p4)/s-
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (128*mt*p1p3**2*p1p4*p2p3*p2p4)/s**2-
     $      (16*mt*p1p3**2*p2p4**2)/s)+
     $      ct3*(-4*mt**3*s*p1p2-4*mt*s*p1p2**2-
     $      4*mt*s*p1p2*p1p3-16*mt*p1p2**2*p1p3+
     $      16*mt*p1p2*p1p3*p1p4+16*mt*p1p2**2*p2p3+
     $      8*mt**3*p1p4*p2p3+16*mt*p1p2*p1p4*p2p3+
     $      8*mt*p1p3*p1p4*p2p3+
     $      (64*mt*p1p2*p1p3*p1p4*p2p3)/s-
     $      (32*mt*p1p3*p1p4**2*p2p3)/s-
     $      (64*mt*p1p2*p1p4*p2p3**2)/s-
     $      (16*mt*p1p4**2*p2p3**2)/s-
     $      (64*mt*p1p3*p1p4**2*p2p3**2)/s**2+
     $      (64*mt*p1p4**2*p2p3**3)/s**2+8*mt**3*p1p3*p2p4+
     $      16*mt*p1p2*p1p3*p2p4+8*mt*p1p3**2*p2p4+
     $      (64*mt*p1p2*p1p3**2*p2p4)/s-
     $      (32*mt*p1p3**2*p1p4*p2p4)/s-
     $      (64*mt*p1p2*p1p3*p2p3*p2p4)/s-
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (128*mt*p1p3**2*p1p4*p2p3*p2p4)/s**2+
     $      (128*mt*p1p3*p1p4*p2p3**2*p2p4)/s**2-
     $      (16*mt*p1p3**2*p2p4**2)/s-
     $      (64*mt*p1p3**3*p2p4**2)/s**2+
     $      (64*mt*p1p3**2*p2p3*p2p4**2)/s**2)+
     $      ct2*(4*mt**5*s+4*mt**3*s*p1p2+16*mt**3*p1p2*p1p3-
     $      4*mt*s*p1p2*p1p3-16*mt**3*p1p3*p1p4-
     $      16*mt**3*p1p2*p2p3-8*mt**3*p1p4*p2p3+
     $      8*mt*p1p3*p1p4*p2p3-
     $      (32*mt**3*p1p3*p1p4*p2p3)/s+
     $      (32*mt**3*p1p4*p2p3**2)/s-8*mt**3*p1p3*p2p4+
     $      8*mt*p1p2*p1p3*p2p4+8*mt*p1p3**2*p2p4-
     $      (32*mt**3*p1p3**2*p2p4)/s-16*mt**3*p2p3*p2p4-
     $      24*mt*p1p2*p2p3*p2p4+
     $      (32*mt**3*p1p3*p2p3*p2p4)/s-
     $      (64*mt*p1p2*p1p3*p2p3*p2p4)/s+
     $      (48*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (64*mt*p1p2*p2p3**2*p2p4)/s+
     $      (48*mt*p1p4*p2p3**2*p2p4)/s+
     $      (128*mt*p1p3*p1p4*p2p3**2*p2p4)/s**2-
     $      (128*mt*p1p4*p2p3**3*p2p4)/s**2-
     $      (16*mt*p1p3**2*p2p4**2)/s+
     $      (48*mt*p1p3*p2p3*p2p4**2)/s+
     $      (128*mt*p1p3**2*p2p3*p2p4**2)/s**2-
     $      (128*mt*p1p3*p2p3**2*p2p4**2)/s**2)

       mat(3,58)=cu1*(-4*mt**3*s*p1p2-4*mt*s*p1p2**2+
     $      16*mt**3*p1p2*p1p3+24*mt*p1p2*p1p3*p1p4-
     $      (64*mt*p1p2*p1p3**2*p1p4)/s-
     $      16*mt**3*p1p2*p2p3+4*mt*s*p1p2*p2p3+
     $      8*mt**3*p1p4*p2p3+8*mt*p1p2*p1p4*p2p3-
     $      (32*mt**3*p1p3*p1p4*p2p3)/s+
     $      (64*mt*p1p2*p1p3*p1p4*p2p3)/s-
     $      (48*mt*p1p3*p1p4**2*p2p3)/s+
     $      (128*mt*p1p3**2*p1p4**2*p2p3)/s**2-
     $      8*mt*p1p4*p2p3**2+(32*mt**3*p1p4*p2p3**2)/s-
     $      (128*mt*p1p3*p1p4**2*p2p3**2)/s**2+
     $      8*mt**3*p1p3*p2p4+16*mt*p1p2*p1p3*p2p4-
     $      (32*mt**3*p1p3**2*p2p4)/s-
     $      (48*mt*p1p3**2*p1p4*p2p4)/s+
     $      (128*mt*p1p3**3*p1p4*p2p4)/s**2-
     $      8*mt*p1p3*p2p3*p2p4+
     $      (32*mt**3*p1p3*p2p3*p2p4)/s-
     $      (16*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (128*mt*p1p3**2*p1p4*p2p3*p2p4)/s**2-
     $      (16*mt*p1p3**2*p2p4**2)/s)+
     $      cu3*(4*mt**5*s+4*mt**3*s*p1p2-16*mt**5*p1p3-
     $      16*mt**3*p1p3*p1p4+(64*mt**3*p1p3**2*p1p4)/s+
     $      16*mt**5*p2p3+4*mt*s*p1p2*p2p3-
     $      8*mt**3*p1p4*p2p3-
     $      (64*mt**3*p1p3*p1p4*p2p3)/s-
     $      8*mt*p1p4*p2p3**2-8*mt**3*p1p3*p2p4-
     $      16*mt**3*p2p3*p2p4-16*mt*p1p2*p2p3*p2p4-
     $      8*mt*p1p3*p2p3*p2p4+
     $      (64*mt**3*p1p3*p2p3*p2p4)/s+
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s-
     $      (256*mt*p1p3**2*p1p4*p2p3*p2p4)/s**2-
     $      (64*mt**3*p2p3**2*p2p4)/s+
     $      (32*mt*p1p4*p2p3**2*p2p4)/s+
     $      (256*mt*p1p3*p1p4*p2p3**2*p2p4)/s**2+
     $      (32*mt*p1p3*p2p3*p2p4**2)/s)+
     $      cu2*(4*mt**5*s+4*mt**3*s*p1p2+16*mt**3*p1p2*p1p3-
     $      16*mt**3*p1p3*p1p4-16*mt**3*p1p2*p2p3+
     $      4*mt*s*p1p2*p2p3+8*mt*p1p2*p1p3*p2p3-
     $      8*mt**3*p1p4*p2p3-(32*mt**3*p1p3*p1p4*p2p3)/s-
     $      8*mt*p1p2*p2p3**2-8*mt*p1p4*p2p3**2+
     $      (32*mt**3*p1p4*p2p3**2)/s-
     $      (16*mt*p1p3*p1p4*p2p3**2)/s+
     $      (16*mt*p1p4*p2p3**3)/s-8*mt**3*p1p3*p2p4-
     $      (32*mt**3*p1p3**2*p2p4)/s-16*mt**3*p2p3*p2p4-
     $      16*mt*p1p2*p2p3*p2p4-8*mt*p1p3*p2p3*p2p4+
     $      (32*mt**3*p1p3*p2p3*p2p4)/s-
     $      (64*mt*p1p2*p1p3*p2p3*p2p4)/s-
     $      (16*mt*p1p3**2*p2p3*p2p4)/s+
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (64*mt*p1p2*p2p3**2*p2p4)/s+
     $      (16*mt*p1p3*p2p3**2*p2p4)/s+
     $      (32*mt*p1p4*p2p3**2*p2p4)/s+
     $      (128*mt*p1p3*p1p4*p2p3**2*p2p4)/s**2-
     $      (128*mt*p1p4*p2p3**3*p2p4)/s**2+
     $      (32*mt*p1p3*p2p3*p2p4**2)/s+
     $      (128*mt*p1p3**2*p2p3*p2p4**2)/s**2-
     $      (128*mt*p1p3*p2p3**2*p2p4**2)/s**2)

       mat(1,59)=cs1*(-4*mt**3*s*p1p3+8*mt**3*p1p3**2-
     $      8*mt**3*p1p3*p1p4+16*mt*p1p3**2*p1p4-
     $      (32*mt*p1p3**3*p1p4)/s+(32*mt*p1p3**2*p1p4**2)/s-
     $      4*mt**3*s*p2p3-8*mt**3*p1p3*p2p3+
     $      8*mt**3*p1p4*p2p3+16*mt*p1p3*p1p4*p2p3+
     $      (32*mt*p1p3**2*p1p4*p2p3)/s-
     $      (32*mt*p1p3*p1p4**2*p2p3)/s)+
     $      cs2*(-4*mt**3*s*p1p3+16*mt*p1p3**2*p1p4-
     $      4*mt**3*s*p2p3-8*mt**3*p1p3*p2p3+
     $      16*mt*p1p3*p1p4*p2p3+
     $      (32*mt*p1p3**2*p1p4*p2p3)/s+8*mt**3*p2p3**2-
     $      (32*mt*p1p3*p1p4*p2p3**2)/s+8*mt**3*p1p3*p2p4-
     $      (32*mt*p1p3**2*p1p4*p2p4)/s-8*mt**3*p2p3*p2p4+
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s)

       mat(2,59)=ct1*(-4*mt**5*s-4*mt**3*s*p1p2+16*mt**5*p1p3-
     $      4*mt**3*s*p1p3+8*mt**3*p1p3**2+
     $      32*mt**3*p1p3*p1p4+16*mt*p1p2*p1p3*p1p4+
     $      16*mt*p1p3**2*p1p4-(128*mt**3*p1p3**2*p1p4)/s-
     $      (32*mt*p1p3**3*p1p4)/s-(64*mt*p1p3**2*p1p4**2)/s+
     $      (256*mt*p1p3**3*p1p4**2)/s**2-16*mt**5*p2p3-
     $      8*mt**3*p1p3*p2p3+8*mt**3*p1p4*p2p3+
     $      (128*mt**3*p1p3*p1p4*p2p3)/s+
     $      (32*mt*p1p3**2*p1p4*p2p3)/s-
     $      (32*mt*p1p3*p1p4**2*p2p3)/s-
     $      (256*mt*p1p3**2*p1p4**2*p2p3)/s**2+
     $      8*mt**3*p1p3*p2p4-(32*mt*p1p3**2*p1p4*p2p4)/s)+
     $      ct3*(-4*mt**5*s-4*mt**3*s*p1p2-4*mt**3*s*p1p3-
     $      16*mt**3*p1p2*p1p3+32*mt**3*p1p3*p1p4+
     $      16*mt*p1p2*p1p3*p1p4+16*mt*p1p3**2*p1p4+
     $      (64*mt*p1p2*p1p3**2*p1p4)/s-
     $      (64*mt*p1p3**2*p1p4**2)/s+16*mt**3*p1p2*p2p3+
     $      8*mt**3*p1p4*p2p3+
     $      (32*mt**3*p1p3*p1p4*p2p3)/s-
     $      (64*mt*p1p2*p1p3*p1p4*p2p3)/s-
     $      (32*mt*p1p3*p1p4**2*p2p3)/s-
     $      (128*mt*p1p3**2*p1p4**2*p2p3)/s**2-
     $      (32*mt**3*p1p4*p2p3**2)/s+
     $      (128*mt*p1p3*p1p4**2*p2p3**2)/s**2+
     $      8*mt**3*p1p3*p2p4+(32*mt**3*p1p3**2*p2p4)/s-
     $      (32*mt*p1p3**2*p1p4*p2p4)/s-
     $      (128*mt*p1p3**3*p1p4*p2p4)/s**2-
     $      (32*mt**3*p1p3*p2p3*p2p4)/s+
     $      (128*mt*p1p3**2*p1p4*p2p3*p2p4)/s**2)+
     $      ct2*(4*mt**3*s*p1p2+4*mt*s*p1p2**2-4*mt**3*s*p1p3+
     $      16*mt*p1p2**2*p1p3-16*mt*p1p2*p1p3*p1p4+
     $      16*mt*p1p3**2*p1p4-16*mt*p1p2**2*p2p3-
     $      8*mt**3*p1p4*p2p3-16*mt*p1p2*p1p4*p2p3-
     $      (64*mt*p1p2*p1p3*p1p4*p2p3)/s+
     $      (32*mt*p1p3*p1p4**2*p2p3)/s+
     $      (64*mt*p1p2*p1p4*p2p3**2)/s+
     $      (16*mt*p1p4**2*p2p3**2)/s+
     $      (64*mt*p1p3*p1p4**2*p2p3**2)/s**2-
     $      (64*mt*p1p4**2*p2p3**3)/s**2-
     $      16*mt*p1p2*p1p3*p2p4-
     $      (64*mt*p1p2*p1p3**2*p2p4)/s-8*mt**3*p2p3*p2p4+
     $      (64*mt*p1p2*p1p3*p2p3*p2p4)/s+
     $      (64*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (128*mt*p1p3**2*p1p4*p2p3*p2p4)/s**2-
     $      (128*mt*p1p3*p1p4*p2p3**2*p2p4)/s**2+
     $      (16*mt*p1p3**2*p2p4**2)/s+
     $      (64*mt*p1p3**3*p2p4**2)/s**2-
     $      (64*mt*p1p3**2*p2p3*p2p4**2)/s**2)

       mat(3,59)=cu1*(-4*mt**5*s-4*mt**3*s*p1p2+16*mt**5*p1p3+
     $      40*mt**3*p1p3*p1p4+16*mt*p1p2*p1p3*p1p4-
     $      (128*mt**3*p1p3**2*p1p4)/s-
     $      (96*mt*p1p3**2*p1p4**2)/s+
     $      (256*mt*p1p3**3*p1p4**2)/s**2-16*mt**5*p2p3+
     $      4*mt**3*s*p2p3-16*mt*p1p3*p1p4*p2p3+
     $      (128*mt**3*p1p3*p1p4*p2p3)/s-
     $      (256*mt*p1p3**2*p1p4**2*p2p3)/s**2+
     $      8*mt**3*p1p3*p2p4-(32*mt*p1p3**2*p1p4*p2p4)/s)+
     $      cu3*(4*mt**3*s*p1p2+4*mt*s*p1p2**2-
     $      16*mt**3*p1p2*p1p3-16*mt*p1p2*p1p3*p1p4+
     $      (64*mt*p1p2*p1p3**2*p1p4)/s+4*mt**3*s*p2p3+
     $      16*mt**3*p1p2*p2p3-8*mt**3*p1p4*p2p3-
     $      16*mt*p1p2*p1p4*p2p3-
     $      16*mt*p1p3*p1p4*p2p3+
     $      (32*mt**3*p1p3*p1p4*p2p3)/s-
     $      (64*mt*p1p2*p1p3*p1p4*p2p3)/s+
     $      (32*mt*p1p3*p1p4**2*p2p3)/s-
     $      (128*mt*p1p3**2*p1p4**2*p2p3)/s**2-
     $      (32*mt**3*p1p4*p2p3**2)/s+
     $      (16*mt*p1p4**2*p2p3**2)/s+
     $      (128*mt*p1p3*p1p4**2*p2p3**2)/s**2-
     $      8*mt**3*p1p3*p2p4-16*mt*p1p2*p1p3*p2p4+
     $      (32*mt**3*p1p3**2*p2p4)/s+
     $      (32*mt*p1p3**2*p1p4*p2p4)/s-
     $      (128*mt*p1p3**3*p1p4*p2p4)/s**2-
     $      (32*mt**3*p1p3*p2p3*p2p4)/s+
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (128*mt*p1p3**2*p1p4*p2p3*p2p4)/s**2+
     $      (16*mt*p1p3**2*p2p4**2)/s)+
     $      cu2*(4*mt**3*s*p1p2+4*mt*s*p1p2**2+
     $      16*mt*p1p2**2*p1p3-16*mt*p1p2*p1p3*p1p4+
     $      4*mt**3*s*p2p3-16*mt*p1p2**2*p2p3+
     $      8*mt**3*p1p3*p2p3-8*mt**3*p1p4*p2p3-
     $      16*mt*p1p2*p1p4*p2p3-16*mt*p1p3*p1p4*p2p3-
     $      (64*mt*p1p2*p1p3*p1p4*p2p3)/s-
     $      (32*mt*p1p3**2*p1p4*p2p3)/s+
     $      (32*mt*p1p3*p1p4**2*p2p3)/s-8*mt**3*p2p3**2+
     $      (64*mt*p1p2*p1p4*p2p3**2)/s+
     $      (32*mt*p1p3*p1p4*p2p3**2)/s+
     $      (16*mt*p1p4**2*p2p3**2)/s+
     $      (64*mt*p1p3*p1p4**2*p2p3**2)/s**2-
     $      (64*mt*p1p4**2*p2p3**3)/s**2-8*mt**3*p1p3*p2p4-
     $      16*mt*p1p2*p1p3*p2p4-
     $      (64*mt*p1p2*p1p3**2*p2p4)/s+
     $      (32*mt*p1p3**2*p1p4*p2p4)/s+
     $      (64*mt*p1p2*p1p3*p2p3*p2p4)/s+
     $      (32*mt*p1p3*p1p4*p2p3*p2p4)/s+
     $      (128*mt*p1p3**2*p1p4*p2p3*p2p4)/s**2-
     $      (128*mt*p1p3*p1p4*p2p3**2*p2p4)/s**2+
     $      (16*mt*p1p3**2*p2p4**2)/s+
     $      (64*mt*p1p3**3*p2p4**2)/s**2-
     $      (64*mt*p1p3**2*p2p3*p2p4**2)/s**2)

       do i=1,3
          do j=60,63
             mat(i,j)=0d0
          end do
       end do

       mat(1,64)=cs1*(-8*mt**4*p1p3 + 8*mt**2*p1p2*p1p3 + 
     $      8*mt**4*p1p4 - 8*mt**2*p1p2*p1p4 + 
     $      8*mt**2*p1p4*p2p3 - 8*mt**2*p1p3*p2p4 + 
     $      (32*mt**2*p1p3*p2p3*p2p4)/s - 
     $      (32*p1p2*p1p3*p2p3*p2p4)/s - 
     $      (32*mt**2*p1p4*p2p3*p2p4)/s + 
     $      (32*p1p2*p1p4*p2p3*p2p4)/s - 
     $      (32*p1p4*p2p3**2*p2p4)/s + 
     $      (32*p1p3*p2p3*p2p4**2)/s) + 
     $      cs2*(8*mt**4*p2p3 - 8*mt**2*p1p2*p2p3 + 
     $      8*mt**2*p1p4*p2p3 - 8*mt**4*p2p4 + 
     $      8*mt**2*p1p2*p2p4 - 8*mt**2*p1p3*p2p4 - 
     $      (32*mt**2*p2p3**2*p2p4)/s + 
     $      (32*p1p2*p2p3**2*p2p4)/s - 
     $      (32*p1p4*p2p3**2*p2p4)/s + 
     $      (32*mt**2*p2p3*p2p4**2)/s - 
     $      (32*p1p2*p2p3*p2p4**2)/s + 
     $      (32*p1p3*p2p3*p2p4**2)/s)
       mat(2,64)=ct1*(2*mt**4*s - 2*mt**2*s*p1p2 - 16*mt**2*p1p2**2 + 
     $      16*p1p2**3 - 8*mt**4*p1p3 + 16*mt**2*p1p2*p1p3 + 
     $      8*mt**2*p1p2*p1p4 - 8*p1p2**2*p2p3 + 
     $      4*mt**2*p1p4*p2p3 + 
     $      (64*mt**2*p1p2*p1p4*p2p3)/s - 
     $      (64*p1p2**2*p1p4*p2p3)/s - 
     $      (16*mt**2*p1p3*p1p4*p2p3)/s - 
     $      (16*mt**2*p1p4**2*p2p3)/s + 
     $      (32*p1p2*p1p4*p2p3**2)/s - 
     $      (64*mt**2*p1p4**2*p2p3**2)/s**2 + 
     $      (64*p1p2*p1p4**2*p2p3**2)/s**2 - 
     $      (32*p1p4**2*p2p3**3)/s**2 - 8*p1p2**2*p2p4 - 
     $      4*mt**2*p1p3*p2p4 + 
     $      (64*mt**2*p1p2*p1p3*p2p4)/s - 
     $      (64*p1p2**2*p1p3*p2p4)/s - 
     $      (16*mt**2*p1p3**2*p2p4)/s - 
     $      (16*mt**2*p1p3*p1p4*p2p4)/s - 
     $      8*mt**2*p2p3*p2p4 + 8*p1p2*p2p3*p2p4 + 
     $      (32*mt**2*p1p3*p2p3*p2p4)/s - 
     $      (32*p1p2*p1p3*p2p3*p2p4)/s - 
     $      (128*mt**2*p1p3*p1p4*p2p3*p2p4)/s**2 + 
     $      (128*p1p2*p1p3*p1p4*p2p3*p2p4)/s**2 - 
     $      (16*p1p4*p2p3**2*p2p4)/s + 
     $      (32*p1p4**2*p2p3**2*p2p4)/s**2 + 
     $      (32*p1p2*p1p3*p2p4**2)/s - 
     $      (64*mt**2*p1p3**2*p2p4**2)/s**2 + 
     $      (64*p1p2*p1p3**2*p2p4**2)/s**2 + 
     $      (16*p1p3*p2p3*p2p4**2)/s + 
     $      (32*p1p3**2*p2p3*p2p4**2)/s**2 - 
     $      (32*p1p3**2*p2p4**3)/s**2) + 
     $      ct3*(2*mt**4*s + 16*mt**4*p1p2 - 2*mt**2*s*p1p2 - 
     $      16*mt**2*p1p2**2 - 8*mt**4*p1p3 + 
     $      8*mt**2*p1p2*p1p4 + 8*mt**2*p1p2*p2p3 + 
     $      4*mt**2*p1p4*p2p3 - (32*mt**4*p1p4*p2p3)/s + 
     $      (32*mt**2*p1p2*p1p4*p2p3)/s - 
     $      (16*mt**2*p1p4**2*p2p3)/s - 
     $      (16*mt**2*p1p4*p2p3**2)/s - 8*p1p2**2*p2p4 - 
     $      4*mt**2*p1p3*p2p4 - (32*mt**4*p1p3*p2p4)/s + 
     $      (32*mt**2*p1p2*p1p3*p2p4)/s - 
     $      (16*mt**2*p1p3*p1p4*p2p4)/s - 
     $      8*mt**2*p2p3*p2p4 + 8*p1p2*p2p3*p2p4 - 
     $      (64*mt**2*p1p2*p2p3*p2p4)/s + 
     $      (64*p1p2**2*p2p3*p2p4)/s + 
     $      (48*mt**2*p1p3*p2p3*p2p4)/s - 
     $      (32*p1p2*p2p3**2*p2p4)/s + 
     $      (128*mt**2*p1p4*p2p3**2*p2p4)/s**2 - 
     $      (16*p1p4*p2p3**2*p2p4)/s - 
     $      (128*p1p2*p1p4*p2p3**2*p2p4)/s**2 + 
     $      (32*p1p4**2*p2p3**2*p2p4)/s**2 + 
     $      (64*p1p4*p2p3**3*p2p4)/s**2 + 
     $      (32*p1p2*p1p3*p2p4**2)/s + 
     $      (128*mt**2*p1p3*p2p3*p2p4**2)/s**2 + 
     $      (16*p1p3*p2p3*p2p4**2)/s - 
     $      (128*p1p2*p1p3*p2p3*p2p4**2)/s**2 - 
     $      (64*p1p3*p2p3**2*p2p4**2)/s**2 - 
     $      (32*p1p3**2*p2p4**3)/s**2) + 
     $      ct2*(-16*mt**6 + 2*mt**4*s + 16*mt**4*p1p2-2*mt**2*s*p1p2- 
     $      8*mt**4*p1p3 - 8*mt**4*p1p4 + 8*mt**2*p1p2*p2p3 + 
     $      4*mt**2*p1p4*p2p3 - (16*mt**2*p1p4*p2p3**2)/s - 
     $      8*mt**4*p2p4 + 16*mt**2*p1p2*p2p4 - 
     $      4*mt**2*p1p3*p2p4 - 8*mt**2*p2p3*p2p4 + 
     $      (128*mt**4*p2p3*p2p4)/s + 8*p1p2*p2p3*p2p4 - 
     $      (128*mt**2*p1p2*p2p3*p2p4)/s + 
     $      (48*mt**2*p1p3*p2p3*p2p4)/s + 
     $      (48*mt**2*p1p4*p2p3*p2p4)/s - 
     $      (32*p1p2*p2p3**2*p2p4)/s - 
     $      (16*p1p4*p2p3**2*p2p4)/s + 
     $      (64*p1p4*p2p3**3*p2p4)/s**2 - 
     $      (16*mt**2*p1p3*p2p4**2)/s + 
     $      (32*mt**2*p2p3*p2p4**2)/s - 
     $      (64*p1p2*p2p3*p2p4**2)/s + 
     $      (16*p1p3*p2p3*p2p4**2)/s - 
     $      (256*mt**2*p2p3**2*p2p4**2)/s**2 + 
     $      (256*p1p2*p2p3**2*p2p4**2)/s**2 - 
     $      (64*p1p3*p2p3**2*p2p4**2)/s**2 - 
     $      (64*p1p4*p2p3**2*p2p4**2)/s**2 + 
     $      (64*p1p3*p2p3*p2p4**3)/s**2)
       mat(3,64)=cu1*(2*mt**4*s - 2*mt**2*s*p1p2 - 16*mt**2*p1p2**2 + 
     $      16*p1p2**3 + 8*mt**2*p1p2*p1p3 - 8*mt**4*p1p4 + 
     $      16*mt**2*p1p2*p1p4 - 8*p1p2**2*p2p3 - 
     $      4*mt**2*p1p4*p2p3 + 
     $      (64*mt**2*p1p2*p1p4*p2p3)/s - 
     $      (64*p1p2**2*p1p4*p2p3)/s - 
     $      (16*mt**2*p1p3*p1p4*p2p3)/s - 
     $      (16*mt**2*p1p4**2*p2p3)/s + 
     $      (32*p1p2*p1p4*p2p3**2)/s - 
     $      (64*mt**2*p1p4**2*p2p3**2)/s**2 + 
     $      (64*p1p2*p1p4**2*p2p3**2)/s**2 - 
     $      (32*p1p4**2*p2p3**3)/s**2 - 8*p1p2**2*p2p4 + 
     $      4*mt**2*p1p3*p2p4 + 
     $      (64*mt**2*p1p2*p1p3*p2p4)/s - 
     $      (64*p1p2**2*p1p3*p2p4)/s - 
     $      (16*mt**2*p1p3**2*p2p4)/s - 
     $      (16*mt**2*p1p3*p1p4*p2p4)/s - 
     $      8*mt**2*p2p3*p2p4 + 8*p1p2*p2p3*p2p4 + 
     $      (32*mt**2*p1p4*p2p3*p2p4)/s - 
     $      (32*p1p2*p1p4*p2p3*p2p4)/s - 
     $      (128*mt**2*p1p3*p1p4*p2p3*p2p4)/s**2 + 
     $      (128*p1p2*p1p3*p1p4*p2p3*p2p4)/s**2 + 
     $      (16*p1p4*p2p3**2*p2p4)/s + 
     $      (32*p1p4**2*p2p3**2*p2p4)/s**2 + 
     $      (32*p1p2*p1p3*p2p4**2)/s - 
     $      (64*mt**2*p1p3**2*p2p4**2)/s**2 + 
     $      (64*p1p2*p1p3**2*p2p4**2)/s**2 - 
     $      (16*p1p3*p2p3*p2p4**2)/s + 
     $      (32*p1p3**2*p2p3*p2p4**2)/s**2 - 
     $      (32*p1p3**2*p2p4**3)/s**2) + 
     $      cu3*(2*mt**4*s + 16*mt**4*p1p2 - 2*mt**2*s*p1p2 - 
     $      16*mt**2*p1p2**2 + 8*mt**2*p1p2*p1p3 - 
     $      8*mt**4*p1p4 - 8*p1p2**2*p2p3 - 
     $      4*mt**2*p1p4*p2p3 - (32*mt**4*p1p4*p2p3)/s + 
     $      (32*mt**2*p1p2*p1p4*p2p3)/s - 
     $      (16*mt**2*p1p3*p1p4*p2p3)/s + 
     $      (32*p1p2*p1p4*p2p3**2)/s - 
     $      (32*p1p4**2*p2p3**3)/s**2 + 8*mt**2*p1p2*p2p4 + 
     $      4*mt**2*p1p3*p2p4 - (32*mt**4*p1p3*p2p4)/s + 
     $      (32*mt**2*p1p2*p1p3*p2p4)/s - 
     $      (16*mt**2*p1p3**2*p2p4)/s - 8*mt**2*p2p3*p2p4 + 
     $      8*p1p2*p2p3*p2p4 - 
     $      (64*mt**2*p1p2*p2p3*p2p4)/s + 
     $      (64*p1p2**2*p2p3*p2p4)/s + 
     $      (48*mt**2*p1p4*p2p3*p2p4)/s + 
     $      (128*mt**2*p1p4*p2p3**2*p2p4)/s**2 + 
     $      (16*p1p4*p2p3**2*p2p4)/s - 
     $      (128*p1p2*p1p4*p2p3**2*p2p4)/s**2 - 
     $      (16*mt**2*p1p3*p2p4**2)/s - 
     $      (32*p1p2*p2p3*p2p4**2)/s + 
     $      (128*mt**2*p1p3*p2p3*p2p4**2)/s**2 - 
     $      (16*p1p3*p2p3*p2p4**2)/s - 
     $      (128*p1p2*p1p3*p2p3*p2p4**2)/s**2 + 
     $      (32*p1p3**2*p2p3*p2p4**2)/s**2 - 
     $      (64*p1p4*p2p3**2*p2p4**2)/s**2 + 
     $      (64*p1p3*p2p3*p2p4**3)/s**2) + 
     $      cu2*(-16*mt**6 + 2*mt**4*s + 16*mt**4*p1p2-2*mt**2*s*p1p2- 
     $      8*mt**4*p1p3 - 8*mt**4*p1p4 - 8*mt**4*p2p3 + 
     $      16*mt**2*p1p2*p2p3 - 4*mt**2*p1p4*p2p3 - 
     $      (16*mt**2*p1p4*p2p3**2)/s + 8*mt**2*p1p2*p2p4 + 
     $      4*mt**2*p1p3*p2p4 - 8*mt**2*p2p3*p2p4 + 
     $      (128*mt**4*p2p3*p2p4)/s + 8*p1p2*p2p3*p2p4 - 
     $      (128*mt**2*p1p2*p2p3*p2p4)/s + 
     $      (48*mt**2*p1p3*p2p3*p2p4)/s + 
     $      (48*mt**2*p1p4*p2p3*p2p4)/s + 
     $      (32*mt**2*p2p3**2*p2p4)/s - 
     $      (64*p1p2*p2p3**2*p2p4)/s + 
     $      (16*p1p4*p2p3**2*p2p4)/s + 
     $      (64*p1p4*p2p3**3*p2p4)/s**2 - 
     $      (16*mt**2*p1p3*p2p4**2)/s - 
     $      (32*p1p2*p2p3*p2p4**2)/s - 
     $      (16*p1p3*p2p3*p2p4**2)/s - 
     $      (256*mt**2*p2p3**2*p2p4**2)/s**2 + 
     $      (256*p1p2*p2p3**2*p2p4**2)/s**2 - 
     $      (64*p1p3*p2p3**2*p2p4**2)/s**2 - 
     $      (64*p1p4*p2p3**2*p2p4**2)/s**2 + 
     $      (64*p1p3*p2p3*p2p4**3)/s**2)

       mat(1,65)=cs1*(-8*mt**2*p1p2*p1p3 + 8*p1p2**2*p1p3 + 
     $      8*mt**2*p1p2*p1p4 - 8*p1p2**2*p1p4 + 
     $      8*p1p2*p1p4*p2p3 + 
     $      (16*mt**2*p1p3*p1p4*p2p3)/s - 
     $      (16*p1p2*p1p3*p1p4*p2p3)/s - 
     $      (16*mt**2*p1p4**2*p2p3)/s + 
     $      (16*p1p2*p1p4**2*p2p3)/s - 
     $      (16*p1p4**2*p2p3**2)/s - 8*p1p2*p1p3*p2p4 + 
     $      (16*mt**2*p1p3**2*p2p4)/s - 
     $      (16*p1p2*p1p3**2*p2p4)/s - 
     $      (16*mt**2*p1p3*p1p4*p2p4)/s + 
     $      (16*p1p2*p1p3*p1p4*p2p4)/s + 
     $      (16*p1p3**2*p2p4**2)/s) + 
     $      cs2*(8*mt**2*p1p2*p2p3 - 8*p1p2**2*p2p3 + 
     $      8*p1p2*p1p4*p2p3 - (16*mt**2*p1p4*p2p3**2)/s + 
     $      (16*p1p2*p1p4*p2p3**2)/s - 
     $      (16*p1p4**2*p2p3**2)/s - 8*mt**2*p1p2*p2p4 + 
     $      8*p1p2**2*p2p4 - 8*p1p2*p1p3*p2p4 - 
     $      (16*mt**2*p1p3*p2p3*p2p4)/s + 
     $      (16*p1p2*p1p3*p2p3*p2p4)/s + 
     $      (16*mt**2*p1p4*p2p3*p2p4)/s - 
     $      (16*p1p2*p1p4*p2p3*p2p4)/s + 
     $      (16*mt**2*p1p3*p2p4**2)/s - 
     $      (16*p1p2*p1p3*p2p4**2)/s + (16*p1p3**2*p2p4**2)/s)
       mat(2,65)=ct1*(-2*mt**4*s - 16*mt**4*p1p2 + 2*mt**2*s*p1p2 + 
     $      16*mt**2*p1p2**2 + 8*mt**4*p1p3 - 
     $      8*mt**2*p1p2*p1p3 + 8*p1p2**2*p1p3 + 
     $      8*mt**4*p1p4 + 8*mt**2*p1p3*p1p4 + 
     $      (64*mt**2*p1p2*p1p3*p1p4)/s - 
     $      (64*p1p2**2*p1p3*p1p4)/s - 
     $      (32*mt**2*p1p3**2*p1p4)/s - 
     $      (32*mt**2*p1p3*p1p4**2)/s - 8*mt**2*p1p2*p2p3 - 
     $      4*mt**2*p1p4*p2p3 + (32*mt**4*p1p4*p2p3)/s - 
     $      (32*mt**2*p1p2*p1p4*p2p3)/s + 
     $      (16*mt**2*p1p3*p1p4*p2p3)/s + 
     $      (16*p1p2*p1p3*p1p4*p2p3)/s - 
     $      (128*mt**2*p1p3*p1p4**2*p2p3)/s**2 + 
     $      (128*p1p2*p1p3*p1p4**2*p2p3)/s**2 + 
     $      (16*mt**2*p1p4*p2p3**2)/s - 
     $      (64*p1p3*p1p4**2*p2p3**2)/s**2 - 
     $      8*mt**2*p1p2*p2p4 - 4*mt**2*p1p3*p2p4 + 
     $      (32*mt**4*p1p3*p2p4)/s - 8*p1p2*p1p3*p2p4 - 
     $      (32*mt**2*p1p2*p1p3*p2p4)/s + 
     $      (16*mt**2*p1p3**2*p2p4)/s - 
     $      (16*p1p2*p1p3**2*p2p4)/s + 
     $      (32*p1p2*p1p3*p1p4*p2p4)/s - 
     $      (128*mt**2*p1p3**2*p1p4*p2p4)/s**2 + 
     $      (128*p1p2*p1p3**2*p1p4*p2p4)/s**2 + 
     $      8*mt**2*p2p3*p2p4 - 
     $      (16*mt**2*p1p3*p2p3*p2p4)/s - 
     $      (16*mt**2*p1p4*p2p3*p2p4)/s - 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s + 
     $      (64*p1p3**2*p1p4*p2p3*p2p4)/s**2 + 
     $      (64*p1p3*p1p4**2*p2p3*p2p4)/s**2 + 
     $      (16*mt**2*p1p3*p2p4**2)/s + (16*p1p3**2*p2p4**2)/s - 
     $      (64*p1p3**2*p1p4*p2p4**2)/s**2) + 
     $      ct3*(16*mt**6 - 2*mt**4*s - 16*mt**4*p1p2 +2*mt**2*s*p1p2- 
     $      8*mt**2*p1p2*p1p3 + 8*mt**4*p1p4 + 
     $      8*mt**2*p1p3*p1p4 - (64*mt**4*p1p3*p1p4)/s + 
     $      (64*mt**2*p1p2*p1p3*p1p4)/s - 
     $      (32*mt**2*p1p3*p1p4**2)/s + 8*mt**4*p2p3 - 
     $      4*mt**2*p1p4*p2p3 - 
     $      (16*mt**2*p1p3*p1p4*p2p3)/s - 
     $      8*mt**2*p1p2*p2p4 - 4*mt**2*p1p3*p2p4 - 
     $      8*p1p2*p1p3*p2p4 + (16*mt**2*p1p3**2*p2p4)/s + 
     $      (32*p1p2*p1p3*p1p4*p2p4)/s + 
     $      8*mt**2*p2p3*p2p4 - (64*mt**4*p2p3*p2p4)/s + 
     $      (64*mt**2*p1p2*p2p3*p2p4)/s + 
     $      (32*p1p2*p1p3*p2p3*p2p4)/s - 
     $      (16*mt**2*p1p4*p2p3*p2p4)/s + 
     $      (256*mt**2*p1p3*p1p4*p2p3*p2p4)/s**2 - 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (256*p1p2*p1p3*p1p4*p2p3*p2p4)/s**2 + 
     $      (64*p1p3*p1p4**2*p2p3*p2p4)/s**2 - 
     $      (32*mt**2*p2p3**2*p2p4)/s + 
     $      (64*p1p3*p1p4*p2p3**2*p2p4)/s**2 + 
     $      (16*mt**2*p1p3*p2p4**2)/s + (16*p1p3**2*p2p4**2)/s - 
     $      (64*p1p3**2*p1p4*p2p4**2)/s**2 - 
     $      (64*p1p3**2*p2p3*p2p4**2)/s**2) + 
     $      ct2*(-2*mt**4*s - 16*mt**4*p1p2 + 2*mt**2*s*p1p2 + 
     $      16*mt**2*p1p2**2 - 8*mt**2*p1p2*p1p3 - 
     $      8*mt**2*p1p2*p1p4 + 8*mt**2*p1p3*p1p4 + 
     $      8*mt**4*p2p3 - 4*mt**2*p1p4*p2p3 + 
     $      (32*mt**4*p1p4*p2p3)/s - 
     $      (32*mt**2*p1p2*p1p4*p2p3)/s - 
     $      (16*mt**2*p1p3*p1p4*p2p3)/s + 
     $      (16*mt**2*p1p4**2*p2p3)/s + 8*mt**4*p2p4 - 
     $      8*mt**2*p1p2*p2p4 + 8*p1p2**2*p2p4 - 
     $      4*mt**2*p1p3*p2p4 + (32*mt**4*p1p3*p2p4)/s - 
     $      8*p1p2*p1p3*p2p4 - 
     $      (32*mt**2*p1p2*p1p3*p2p4)/s + 
     $      (16*mt**2*p1p3**2*p2p4)/s - 
     $      (16*mt**2*p1p3*p1p4*p2p4)/s + 8*mt**2*p2p3*p2p4 + 
     $      (64*mt**2*p1p2*p2p3*p2p4)/s - 
     $      (64*p1p2**2*p2p3*p2p4)/s + 
     $      (32*p1p2*p1p3*p2p3*p2p4)/s + 
     $      (16*mt**2*p1p4*p2p3*p2p4)/s + 
     $      (16*p1p2*p1p4*p2p3*p2p4)/s - 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (32*mt**2*p2p3**2*p2p4)/s - 
     $      (128*mt**2*p1p4*p2p3**2*p2p4)/s**2 + 
     $      (128*p1p2*p1p4*p2p3**2*p2p4)/s**2 + 
     $      (64*p1p3*p1p4*p2p3**2*p2p4)/s**2 - 
     $      (64*p1p4**2*p2p3**2*p2p4)/s**2 + 
     $      (16*mt**2*p1p3*p2p4**2)/s - 
     $      (16*p1p2*p1p3*p2p4**2)/s + 
     $      (16*p1p3**2*p2p4**2)/s - (32*mt**2*p2p3*p2p4**2)/s - 
     $      (128*mt**2*p1p3*p2p3*p2p4**2)/s**2 + 
     $      (128*p1p2*p1p3*p2p3*p2p4**2)/s**2 - 
     $      (64*p1p3**2*p2p3*p2p4**2)/s**2 + 
     $      (64*p1p3*p1p4*p2p3*p2p4**2)/s**2)
       mat(3,65)=cu1*(2*mt**4*s - 16*mt**4*p1p2 + 2*mt**2*s*p1p2 + 
     $      16*mt**2*p1p2**2 - 4*s*p1p2**2 + 8*p1p2**2*p1p3 - 
     $      8*mt**2*p1p2*p1p4 + 16*p1p2**2*p1p4 - 
     $      8*mt**2*p1p3*p1p4 + 
     $      (64*mt**2*p1p2*p1p3*p1p4)/s - 
     $      (64*p1p2**2*p1p3*p1p4)/s - 8*mt**2*p1p2*p2p3 - 
     $      4*mt**2*p1p4*p2p3 + (32*mt**4*p1p4*p2p3)/s + 
     $      8*p1p2*p1p4*p2p3 - 
     $      (32*mt**2*p1p2*p1p4*p2p3)/s + 
     $      (16*mt**2*p1p4**2*p2p3)/s - 
     $      (48*p1p2*p1p4**2*p2p3)/s - 
     $      (128*mt**2*p1p3*p1p4**2*p2p3)/s**2 + 
     $      (128*p1p2*p1p3*p1p4**2*p2p3)/s**2 + 
     $      (16*mt**2*p1p4*p2p3**2)/s - 
     $      (32*p1p3*p1p4**2*p2p3**2)/s**2 + 
     $      (32*p1p4**3*p2p3**2)/s**2 - 8*mt**2*p1p2*p2p4 - 
     $      4*mt**2*p1p3*p2p4 + (32*mt**4*p1p3*p2p4)/s + 
     $      16*p1p2*p1p3*p2p4 - 
     $      (32*mt**2*p1p2*p1p3*p2p4)/s - 
     $      (32*p1p2*p1p3**2*p2p4)/s + 
     $      (16*mt**2*p1p3*p1p4*p2p4)/s - 
     $      (16*p1p2*p1p3*p1p4*p2p4)/s - 
     $      (128*mt**2*p1p3**2*p1p4*p2p4)/s**2 + 
     $      (128*p1p2*p1p3**2*p1p4*p2p4)/s**2 - 
     $      8*mt**2*p2p3*p2p4 + 
     $      (16*mt**2*p1p3*p2p3*p2p4)/s + 
     $      (16*mt**2*p1p4*p2p3*p2p4)/s + 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s + 
     $      (16*mt**2*p1p3*p2p4**2)/s - (16*p1p3**2*p2p4**2)/s + 
     $      (32*p1p3**3*p2p4**2)/s**2 - 
     $      (32*p1p3**2*p1p4*p2p4**2)/s**2) + 
     $      cu3*(2*mt**4*s + 2*mt**2*s*p1p2 + 16*mt**2*p1p2**2 - 
     $      4*s*p1p2**2 - 16*p1p2**3 + 8*p1p2**2*p1p3 - 
     $      8*mt**2*p1p2*p1p4 - 8*mt**2*p1p3*p1p4 - 
     $      8*mt**2*p1p2*p2p3 - 4*mt**2*p1p4*p2p3 + 
     $      8*p1p2*p1p4*p2p3 - 
     $      (64*mt**2*p1p2*p1p4*p2p3)/s + 
     $      (64*p1p2**2*p1p4*p2p3)/s + 
     $      (16*mt**2*p1p4**2*p2p3)/s + 
     $      (16*mt**2*p1p4*p2p3**2)/s + 
     $      (64*mt**2*p1p4**2*p2p3**2)/s**2 - 
     $      (64*p1p2*p1p4**2*p2p3**2)/s**2 - 
     $      (32*p1p3*p1p4**2*p2p3**2)/s**2 + 
     $      8*p1p2**2*p2p4 - 4*mt**2*p1p3*p2p4 + 
     $      16*p1p2*p1p3*p2p4 - 
     $      (64*mt**2*p1p2*p1p3*p2p4)/s + 
     $      (64*p1p2**2*p1p3*p2p4)/s - 
     $      (32*p1p2*p1p3**2*p2p4)/s + 
     $      (16*mt**2*p1p3*p1p4*p2p4)/s - 
     $      8*mt**2*p2p3*p2p4 + 
     $      (16*mt**2*p1p3*p2p3*p2p4)/s + 
     $      (128*mt**2*p1p3*p1p4*p2p3*p2p4)/s**2 + 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (128*p1p2*p1p3*p1p4*p2p3*p2p4)/s**2 - 
     $      (32*p1p4**2*p2p3**2*p2p4)/s**2 - 
     $      (32*p1p2*p1p3*p2p4**2)/s + 
     $      (64*mt**2*p1p3**2*p2p4**2)/s**2 - 
     $      (16*p1p3**2*p2p4**2)/s - 
     $      (64*p1p2*p1p3**2*p2p4**2)/s**2 + 
     $      (32*p1p3**3*p2p4**2)/s**2 + (32*p1p3**2*p2p4**3)/s**2) + 
     $      cu2*(2*mt**4*s - 16*mt**4*p1p2 + 2*mt**2*s*p1p2 + 
     $      16*mt**2*p1p2**2 - 4*s*p1p2**2 - 8*mt**2*p1p2*p1p3 - 
     $      8*mt**2*p1p2*p1p4 - 8*mt**2*p1p3*p1p4 - 
     $      8*mt**2*p1p2*p2p3 + 16*p1p2**2*p2p3 - 
     $      4*mt**2*p1p4*p2p3 + (32*mt**4*p1p4*p2p3)/s + 
     $      8*p1p2*p1p4*p2p3 - 
     $      (32*mt**2*p1p2*p1p4*p2p3)/s + 
     $      (16*mt**2*p1p3*p1p4*p2p3)/s + 
     $      (16*mt**2*p1p4**2*p2p3)/s + 
     $      (16*mt**2*p1p4*p2p3**2)/s - 
     $      (48*p1p2*p1p4*p2p3**2)/s + 
     $      (32*p1p4**2*p2p3**3)/s**2 + 8*p1p2**2*p2p4 - 
     $      4*mt**2*p1p3*p2p4 + (32*mt**4*p1p3*p2p4)/s + 
     $      16*p1p2*p1p3*p2p4 - 
     $      (32*mt**2*p1p2*p1p3*p2p4)/s + 
     $      (16*mt**2*p1p3**2*p2p4)/s + 
     $      (16*mt**2*p1p3*p1p4*p2p4)/s - 8*mt**2*p2p3*p2p4 + 
     $      (64*mt**2*p1p2*p2p3*p2p4)/s - 
     $      (64*p1p2**2*p2p3*p2p4)/s + 
     $      (16*mt**2*p1p3*p2p3*p2p4)/s - 
     $      (16*p1p2*p1p3*p2p3*p2p4)/s + 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (128*mt**2*p1p4*p2p3**2*p2p4)/s**2 + 
     $      (128*p1p2*p1p4*p2p3**2*p2p4)/s**2 - 
     $      (32*p1p4**2*p2p3**2*p2p4)/s**2 - 
     $      (32*p1p2*p1p3*p2p4**2)/s - 
     $      (16*p1p3**2*p2p4**2)/s - 
     $      (128*mt**2*p1p3*p2p3*p2p4**2)/s**2 + 
     $      (128*p1p2*p1p3*p2p3*p2p4**2)/s**2 - 
     $      (32*p1p3**2*p2p3*p2p4**2)/s**2 + 
     $      (32*p1p3**2*p2p4**3)/s**2)

       do i=1,3
          mat(i,66)=0d0
       end do

       mat(1,67)=cs1*(-8*mt**2*p1p2*p1p3 + 8*p1p2**2*p1p3 + 
     $      8*mt**2*p1p2*p1p4 - 8*p1p2**2*p1p4 + 
     $      8*p1p2*p1p4*p2p3 + 
     $      (16*mt**2*p1p3*p1p4*p2p3)/s - 
     $      (16*p1p2*p1p3*p1p4*p2p3)/s - 
     $      (16*mt**2*p1p4**2*p2p3)/s + 
     $      (16*p1p2*p1p4**2*p2p3)/s - 
     $      (16*p1p4**2*p2p3**2)/s - 8*p1p2*p1p3*p2p4 + 
     $      (16*mt**2*p1p3**2*p2p4)/s - 
     $      (16*p1p2*p1p3**2*p2p4)/s - 
     $      (16*mt**2*p1p3*p1p4*p2p4)/s + 
     $      (16*p1p2*p1p3*p1p4*p2p4)/s + 
     $      (16*p1p3**2*p2p4**2)/s) + 
     $      cs2*(8*mt**2*p1p2*p2p3 - 8*p1p2**2*p2p3 + 
     $      8*p1p2*p1p4*p2p3 - (16*mt**2*p1p4*p2p3**2)/s + 
     $      (16*p1p2*p1p4*p2p3**2)/s - 
     $      (16*p1p4**2*p2p3**2)/s - 8*mt**2*p1p2*p2p4 + 
     $      8*p1p2**2*p2p4 - 8*p1p2*p1p3*p2p4 - 
     $      (16*mt**2*p1p3*p2p3*p2p4)/s + 
     $      (16*p1p2*p1p3*p2p3*p2p4)/s + 
     $      (16*mt**2*p1p4*p2p3*p2p4)/s - 
     $      (16*p1p2*p1p4*p2p3*p2p4)/s + 
     $      (16*mt**2*p1p3*p2p4**2)/s - 
     $      (16*p1p2*p1p3*p2p4**2)/s + (16*p1p3**2*p2p4**2)/s)
       mat(2,67)=ct1*(2*mt**4*s - 16*mt**4*p1p2 + 2*mt**2*s*p1p2 + 
     $      16*mt**2*p1p2**2 - 4*s*p1p2**2 - 8*mt**2*p1p2*p1p3 + 
     $      16*p1p2**2*p1p3 + 8*p1p2**2*p1p4 - 
     $      8*mt**2*p1p3*p1p4 + 
     $      (64*mt**2*p1p2*p1p3*p1p4)/s - 
     $      (64*p1p2**2*p1p3*p1p4)/s - 8*mt**2*p1p2*p2p3 - 
     $      4*mt**2*p1p4*p2p3 + (32*mt**4*p1p4*p2p3)/s + 
     $      16*p1p2*p1p4*p2p3 - 
     $      (32*mt**2*p1p2*p1p4*p2p3)/s + 
     $      (16*mt**2*p1p3*p1p4*p2p3)/s - 
     $      (16*p1p2*p1p3*p1p4*p2p3)/s - 
     $      (32*p1p2*p1p4**2*p2p3)/s - 
     $      (128*mt**2*p1p3*p1p4**2*p2p3)/s**2 + 
     $      (128*p1p2*p1p3*p1p4**2*p2p3)/s**2 + 
     $      (16*mt**2*p1p4*p2p3**2)/s - (16*p1p4**2*p2p3**2)/s - 
     $      (32*p1p3*p1p4**2*p2p3**2)/s**2 + 
     $      (32*p1p4**3*p2p3**2)/s**2 - 8*mt**2*p1p2*p2p4 - 
     $      4*mt**2*p1p3*p2p4 + (32*mt**4*p1p3*p2p4)/s + 
     $      8*p1p2*p1p3*p2p4 - 
     $      (32*mt**2*p1p2*p1p3*p2p4)/s + 
     $      (16*mt**2*p1p3**2*p2p4)/s - 
     $      (48*p1p2*p1p3**2*p2p4)/s - 
     $      (128*mt**2*p1p3**2*p1p4*p2p4)/s**2 + 
     $      (128*p1p2*p1p3**2*p1p4*p2p4)/s**2 - 
     $      8*mt**2*p2p3*p2p4 + 
     $      (16*mt**2*p1p3*p2p3*p2p4)/s + 
     $      (16*mt**2*p1p4*p2p3*p2p4)/s + 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s + 
     $      (16*mt**2*p1p3*p2p4**2)/s + 
     $      (32*p1p3**3*p2p4**2)/s**2 - 
     $      (32*p1p3**2*p1p4*p2p4**2)/s**2) + 
     $      ct3*(2*mt**4*s + 2*mt**2*s*p1p2 + 16*mt**2*p1p2**2 - 
     $      4*s*p1p2**2 - 16*p1p2**3 - 8*mt**2*p1p2*p1p3 + 
     $      8*p1p2**2*p1p4 - 8*mt**2*p1p3*p1p4 + 
     $      8*p1p2**2*p2p3 - 4*mt**2*p1p4*p2p3 + 
     $      16*p1p2*p1p4*p2p3 - 
     $      (64*mt**2*p1p2*p1p4*p2p3)/s + 
     $      (64*p1p2**2*p1p4*p2p3)/s + 
     $      (16*mt**2*p1p3*p1p4*p2p3)/s - 
     $      (32*p1p2*p1p4**2*p2p3)/s - 
     $      (32*p1p2*p1p4*p2p3**2)/s + 
     $      (64*mt**2*p1p4**2*p2p3**2)/s**2 - 
     $      (16*p1p4**2*p2p3**2)/s - 
     $      (64*p1p2*p1p4**2*p2p3**2)/s**2 + 
     $      (32*p1p4**3*p2p3**2)/s**2 + 
     $      (32*p1p4**2*p2p3**3)/s**2 - 8*mt**2*p1p2*p2p4 - 
     $      4*mt**2*p1p3*p2p4 + 8*p1p2*p1p3*p2p4 - 
     $      (64*mt**2*p1p2*p1p3*p2p4)/s + 
     $      (64*p1p2**2*p1p3*p2p4)/s + 
     $      (16*mt**2*p1p3**2*p2p4)/s - 8*mt**2*p2p3*p2p4 + 
     $      (16*mt**2*p1p4*p2p3*p2p4)/s + 
     $      (128*mt**2*p1p3*p1p4*p2p3*p2p4)/s**2 + 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (128*p1p2*p1p3*p1p4*p2p3*p2p4)/s**2 + 
     $      (16*mt**2*p1p3*p2p4**2)/s + 
     $      (64*mt**2*p1p3**2*p2p4**2)/s**2 - 
     $      (64*p1p2*p1p3**2*p2p4**2)/s**2 - 
     $      (32*p1p3**2*p1p4*p2p4**2)/s**2 - 
     $      (32*p1p3**2*p2p3*p2p4**2)/s**2) + 
     $      ct2*(2*mt**4*s - 16*mt**4*p1p2 + 2*mt**2*s*p1p2 + 
     $      16*mt**2*p1p2**2 - 4*s*p1p2**2 - 8*mt**2*p1p2*p1p3 - 
     $      8*mt**2*p1p2*p1p4 - 8*mt**2*p1p3*p1p4 + 
     $      8*p1p2**2*p2p3 - 4*mt**2*p1p4*p2p3 + 
     $      (32*mt**4*p1p4*p2p3)/s + 16*p1p2*p1p4*p2p3 - 
     $      (32*mt**2*p1p2*p1p4*p2p3)/s + 
     $      (16*mt**2*p1p3*p1p4*p2p3)/s + 
     $      (16*mt**2*p1p4**2*p2p3)/s - 
     $      (32*p1p2*p1p4*p2p3**2)/s - 
     $      (16*p1p4**2*p2p3**2)/s + (32*p1p4**2*p2p3**3)/s**2 - 
     $      8*mt**2*p1p2*p2p4 + 16*p1p2**2*p2p4 - 
     $      4*mt**2*p1p3*p2p4 + (32*mt**4*p1p3*p2p4)/s + 
     $      8*p1p2*p1p3*p2p4 - 
     $      (32*mt**2*p1p2*p1p3*p2p4)/s + 
     $      (16*mt**2*p1p3**2*p2p4)/s + 
     $      (16*mt**2*p1p3*p1p4*p2p4)/s - 8*mt**2*p2p3*p2p4 + 
     $      (64*mt**2*p1p2*p2p3*p2p4)/s - 
     $      (64*p1p2**2*p2p3*p2p4)/s + 
     $      (16*mt**2*p1p4*p2p3*p2p4)/s - 
     $      (16*p1p2*p1p4*p2p3*p2p4)/s + 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (128*mt**2*p1p4*p2p3**2*p2p4)/s**2 + 
     $      (128*p1p2*p1p4*p2p3**2*p2p4)/s**2 - 
     $      (32*p1p4**2*p2p3**2*p2p4)/s**2 + 
     $      (16*mt**2*p1p3*p2p4**2)/s - 
     $      (48*p1p2*p1p3*p2p4**2)/s - 
     $      (128*mt**2*p1p3*p2p3*p2p4**2)/s**2 + 
     $      (128*p1p2*p1p3*p2p3*p2p4**2)/s**2 - 
     $      (32*p1p3**2*p2p3*p2p4**2)/s**2 + 
     $      (32*p1p3**2*p2p4**3)/s**2)
       mat(3,67)=cu1*(-2*mt**4*s - 16*mt**4*p1p2 + 2*mt**2*s*p1p2 + 
     $      16*mt**2*p1p2**2 + 8*mt**4*p1p3 + 8*mt**4*p1p4 - 
     $      8*mt**2*p1p2*p1p4 + 8*p1p2**2*p1p4 + 
     $      8*mt**2*p1p3*p1p4 + 
     $      (64*mt**2*p1p2*p1p3*p1p4)/s - 
     $      (64*p1p2**2*p1p3*p1p4)/s - 
     $      (32*mt**2*p1p3**2*p1p4)/s - 
     $      (32*mt**2*p1p3*p1p4**2)/s - 8*mt**2*p1p2*p2p3 - 
     $      4*mt**2*p1p4*p2p3 + (32*mt**4*p1p4*p2p3)/s - 
     $      8*p1p2*p1p4*p2p3 - 
     $      (32*mt**2*p1p2*p1p4*p2p3)/s + 
     $      (32*p1p2*p1p3*p1p4*p2p3)/s + 
     $      (16*mt**2*p1p4**2*p2p3)/s - 
     $      (16*p1p2*p1p4**2*p2p3)/s - 
     $      (128*mt**2*p1p3*p1p4**2*p2p3)/s**2 + 
     $      (128*p1p2*p1p3*p1p4**2*p2p3)/s**2 + 
     $      (16*mt**2*p1p4*p2p3**2)/s + (16*p1p4**2*p2p3**2)/s - 
     $      (64*p1p3*p1p4**2*p2p3**2)/s**2 - 
     $      8*mt**2*p1p2*p2p4 - 4*mt**2*p1p3*p2p4 + 
     $      (32*mt**4*p1p3*p2p4)/s - 
     $      (32*mt**2*p1p2*p1p3*p2p4)/s + 
     $      (16*mt**2*p1p3*p1p4*p2p4)/s + 
     $      (16*p1p2*p1p3*p1p4*p2p4)/s - 
     $      (128*mt**2*p1p3**2*p1p4*p2p4)/s**2 + 
     $      (128*p1p2*p1p3**2*p1p4*p2p4)/s**2 + 
     $      8*mt**2*p2p3*p2p4 - 
     $      (16*mt**2*p1p3*p2p3*p2p4)/s - 
     $      (16*mt**2*p1p4*p2p3*p2p4)/s - 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s + 
     $      (64*p1p3**2*p1p4*p2p3*p2p4)/s**2 + 
     $      (64*p1p3*p1p4**2*p2p3*p2p4)/s**2 + 
     $      (16*mt**2*p1p3*p2p4**2)/s - 
     $      (64*p1p3**2*p1p4*p2p4**2)/s**2) + 
     $      cu3*(16*mt**6 - 2*mt**4*s - 16*mt**4*p1p2 +2*mt**2*s*p1p2+ 
     $      8*mt**4*p1p3 - 8*mt**2*p1p2*p1p4 + 
     $      8*mt**2*p1p3*p1p4 - (64*mt**4*p1p3*p1p4)/s + 
     $      (64*mt**2*p1p2*p1p3*p1p4)/s - 
     $      (32*mt**2*p1p3**2*p1p4)/s - 8*mt**2*p1p2*p2p3 - 
     $      4*mt**2*p1p4*p2p3 - 8*p1p2*p1p4*p2p3 + 
     $      (32*p1p2*p1p3*p1p4*p2p3)/s + 
     $      (16*mt**2*p1p4**2*p2p3)/s + 
     $      (16*mt**2*p1p4*p2p3**2)/s + (16*p1p4**2*p2p3**2)/s - 
     $      (64*p1p3*p1p4**2*p2p3**2)/s**2 + 8*mt**4*p2p4 - 
     $      4*mt**2*p1p3*p2p4 - 
     $      (16*mt**2*p1p3*p1p4*p2p4)/s + 
     $      8*mt**2*p2p3*p2p4 - (64*mt**4*p2p3*p2p4)/s + 
     $      (64*mt**2*p1p2*p2p3*p2p4)/s - 
     $      (16*mt**2*p1p3*p2p3*p2p4)/s + 
     $      (32*p1p2*p1p4*p2p3*p2p4)/s + 
     $      (256*mt**2*p1p3*p1p4*p2p3*p2p4)/s**2 - 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (256*p1p2*p1p3*p1p4*p2p3*p2p4)/s**2 + 
     $      (64*p1p3**2*p1p4*p2p3*p2p4)/s**2 - 
     $      (64*p1p4**2*p2p3**2*p2p4)/s**2 - 
     $      (32*mt**2*p2p3*p2p4**2)/s + 
     $      (64*p1p3*p1p4*p2p3*p2p4**2)/s**2) + 
     $      cu2*(-2*mt**4*s - 16*mt**4*p1p2 + 2*mt**2*s*p1p2 + 
     $      16*mt**2*p1p2**2 - 8*mt**2*p1p2*p1p3 - 
     $      8*mt**2*p1p2*p1p4 + 8*mt**2*p1p3*p1p4 + 
     $      8*mt**4*p2p3 - 8*mt**2*p1p2*p2p3 + 
     $      8*p1p2**2*p2p3 - 4*mt**2*p1p4*p2p3 + 
     $      (32*mt**4*p1p4*p2p3)/s - 8*p1p2*p1p4*p2p3 - 
     $      (32*mt**2*p1p2*p1p4*p2p3)/s - 
     $      (16*mt**2*p1p3*p1p4*p2p3)/s + 
     $      (16*mt**2*p1p4**2*p2p3)/s + 
     $      (16*mt**2*p1p4*p2p3**2)/s - 
     $      (16*p1p2*p1p4*p2p3**2)/s + 
     $      (16*p1p4**2*p2p3**2)/s + 8*mt**4*p2p4 - 
     $      4*mt**2*p1p3*p2p4 + (32*mt**4*p1p3*p2p4)/s - 
     $      (32*mt**2*p1p2*p1p3*p2p4)/s + 
     $      (16*mt**2*p1p3**2*p2p4)/s - 
     $      (16*mt**2*p1p3*p1p4*p2p4)/s + 8*mt**2*p2p3*p2p4 + 
     $      (64*mt**2*p1p2*p2p3*p2p4)/s - 
     $      (64*p1p2**2*p2p3*p2p4)/s + 
     $      (16*mt**2*p1p3*p2p3*p2p4)/s + 
     $      (16*p1p2*p1p3*p2p3*p2p4)/s + 
     $      (32*p1p2*p1p4*p2p3*p2p4)/s - 
     $      (16*p1p3*p1p4*p2p3*p2p4)/s - 
     $      (32*mt**2*p2p3**2*p2p4)/s - 
     $      (128*mt**2*p1p4*p2p3**2*p2p4)/s**2 + 
     $      (128*p1p2*p1p4*p2p3**2*p2p4)/s**2 + 
     $      (64*p1p3*p1p4*p2p3**2*p2p4)/s**2 - 
     $      (64*p1p4**2*p2p3**2*p2p4)/s**2 - 
     $      (32*mt**2*p2p3*p2p4**2)/s - 
     $      (128*mt**2*p1p3*p2p3*p2p4**2)/s**2 + 
     $      (128*p1p2*p1p3*p2p3*p2p4**2)/s**2 - 
     $      (64*p1p3**2*p2p3*p2p4**2)/s**2 + 
     $      (64*p1p3*p1p4*p2p3*p2p4**2)/s**2)

       mat(1,68)=cs1*(-8*mt**4*p1p3 + 8*mt**2*p1p2*p1p3 + 
     $      8*mt**4*p1p4 - 8*mt**2*p1p2*p1p4 + 
     $      (32*mt**2*p1p3**2*p1p4)/s - 
     $      (32*p1p2*p1p3**2*p1p4)/s - 
     $      (32*mt**2*p1p3*p1p4**2)/s + 
     $      (32*p1p2*p1p3*p1p4**2)/s + 8*mt**2*p1p4*p2p3 - 
     $      (32*p1p3*p1p4**2*p2p3)/s - 8*mt**2*p1p3*p2p4 + 
     $      (32*p1p3**2*p1p4*p2p4)/s) + 
     $      cs2*(8*mt**4*p2p3 - 8*mt**2*p1p2*p2p3 + 
     $      8*mt**2*p1p4*p2p3 - (32*mt**2*p1p3*p1p4*p2p3)/s + 
     $      (32*p1p2*p1p3*p1p4*p2p3)/s - 
     $      (32*p1p3*p1p4**2*p2p3)/s - 8*mt**4*p2p4 + 
     $      8*mt**2*p1p2*p2p4 - 8*mt**2*p1p3*p2p4 + 
     $      (32*mt**2*p1p3*p1p4*p2p4)/s - 
     $      (32*p1p2*p1p3*p1p4*p2p4)/s + 
     $      (32*p1p3**2*p1p4*p2p4)/s)
       mat(2,68)=ct1*(-16*mt**6 + 2*mt**4*s + 16*mt**4*p1p2 - 
     $      2*mt**2*s*p1p2 - 8*mt**4*p1p3 + 16*mt**2*p1p2*p1p3 + 
     $      8*mt**2*p1p2*p1p4 - 8*mt**2*p1p3*p1p4 + 
     $      (128*mt**4*p1p3*p1p4)/s + 8*p1p2*p1p3*p1p4 - 
     $      (128*mt**2*p1p2*p1p3*p1p4)/s + 
     $      (32*mt**2*p1p3**2*p1p4)/s - 
     $      (64*p1p2*p1p3**2*p1p4)/s - 
     $      (32*p1p2*p1p3*p1p4**2)/s - 
     $      (256*mt**2*p1p3**2*p1p4**2)/s**2 + 
     $      (256*p1p2*p1p3**2*p1p4**2)/s**2 - 8*mt**4*p2p3 + 
     $      4*mt**2*p1p4*p2p3 + 
     $      (48*mt**2*p1p3*p1p4*p2p3)/s - 
     $      (16*mt**2*p1p4**2*p2p3)/s - 
     $      (16*p1p3*p1p4**2*p2p3)/s - 
     $      (64*p1p3**2*p1p4**2*p2p3)/s**2 + 
     $      (64*p1p3*p1p4**3*p2p3)/s**2 - 8*mt**4*p2p4 - 
     $      4*mt**2*p1p3*p2p4 - (16*mt**2*p1p3**2*p2p4)/s + 
     $      (48*mt**2*p1p3*p1p4*p2p4)/s + 
     $      (16*p1p3**2*p1p4*p2p4)/s + 
     $      (64*p1p3**3*p1p4*p2p4)/s**2 - 
     $      (64*p1p3**2*p1p4**2*p2p4)/s**2) + 
     $      ct3*(2*mt**4*s + 16*mt**4*p1p2 - 2*mt**2*s*p1p2 - 
     $      16*mt**2*p1p2**2 - 8*p1p2**2*p1p3 + 
     $      8*mt**2*p1p2*p1p4 - 8*mt**2*p1p3*p1p4 + 
     $      8*p1p2*p1p3*p1p4 - 
     $      (64*mt**2*p1p2*p1p3*p1p4)/s + 
     $      (64*p1p2**2*p1p3*p1p4)/s - 
     $      (32*p1p2*p1p3*p1p4**2)/s + 8*mt**2*p1p2*p2p3 + 
     $      4*mt**2*p1p4*p2p3 - (32*mt**4*p1p4*p2p3)/s + 
     $      (32*mt**2*p1p2*p1p4*p2p3)/s - 
     $      (16*mt**2*p1p4**2*p2p3)/s + 
     $      (128*mt**2*p1p3*p1p4**2*p2p3)/s**2 - 
     $      (16*p1p3*p1p4**2*p2p3)/s - 
     $      (128*p1p2*p1p3*p1p4**2*p2p3)/s**2 + 
     $      (64*p1p3*p1p4**3*p2p3)/s**2 - 
     $      (16*mt**2*p1p4*p2p3**2)/s + 
     $      (32*p1p3*p1p4**2*p2p3**2)/s**2 - 8*mt**4*p2p4 - 
     $      4*mt**2*p1p3*p2p4 - (32*mt**4*p1p3*p2p4)/s + 
     $      (32*mt**2*p1p2*p1p3*p2p4)/s + 
     $      (32*p1p2*p1p3**2*p2p4)/s + 
     $      (48*mt**2*p1p3*p1p4*p2p4)/s + 
     $      (128*mt**2*p1p3**2*p1p4*p2p4)/s**2 + 
     $      (16*p1p3**2*p1p4*p2p4)/s - 
     $      (128*p1p2*p1p3**2*p1p4*p2p4)/s**2 - 
     $      (64*p1p3**2*p1p4**2*p2p4)/s**2 - 
     $      (16*mt**2*p1p3*p2p3*p2p4)/s - 
     $      (32*p1p3**3*p2p4**2)/s**2) + 
     $      ct2*(2*mt**4*s - 2*mt**2*s*p1p2 - 16*mt**2*p1p2**2 + 
     $      16*p1p2**3 - 8*p1p2**2*p1p3 - 8*p1p2**2*p1p4 - 
     $      8*mt**2*p1p3*p1p4 + 8*p1p2*p1p3*p1p4 + 
     $      8*mt**2*p1p2*p2p3 + 4*mt**2*p1p4*p2p3 + 
     $      (64*mt**2*p1p2*p1p4*p2p3)/s - 
     $      (64*p1p2**2*p1p4*p2p3)/s + 
     $      (32*p1p2*p1p4**2*p2p3)/s - 
     $      (16*p1p3*p1p4**2*p2p3)/s - 
     $      (16*mt**2*p1p4*p2p3**2)/s - 
     $      (64*mt**2*p1p4**2*p2p3**2)/s**2 + 
     $      (64*p1p2*p1p4**2*p2p3**2)/s**2 + 
     $      (32*p1p3*p1p4**2*p2p3**2)/s**2 - 
     $      (32*p1p4**3*p2p3**2)/s**2 - 8*mt**4*p2p4 + 
     $      16*mt**2*p1p2*p2p4 - 4*mt**2*p1p3*p2p4 + 
     $      (64*mt**2*p1p2*p1p3*p2p4)/s - 
     $      (64*p1p2**2*p1p3*p2p4)/s + 
     $      (32*p1p2*p1p3**2*p2p4)/s + 
     $      (32*mt**2*p1p3*p1p4*p2p4)/s - 
     $      (32*p1p2*p1p3*p1p4*p2p4)/s + 
     $      (16*p1p3**2*p1p4*p2p4)/s - 
     $      (16*mt**2*p1p3*p2p3*p2p4)/s - 
     $      (16*mt**2*p1p4*p2p3*p2p4)/s - 
     $      (128*mt**2*p1p3*p1p4*p2p3*p2p4)/s**2 + 
     $      (128*p1p2*p1p3*p1p4*p2p3*p2p4)/s**2 - 
     $      (16*mt**2*p1p3*p2p4**2)/s - 
     $      (64*mt**2*p1p3**2*p2p4**2)/s**2 + 
     $      (64*p1p2*p1p3**2*p2p4**2)/s**2 - 
     $      (32*p1p3**3*p2p4**2)/s**2 + 
     $      (32*p1p3**2*p1p4*p2p4**2)/s**2)
       mat(3,68)=cu1*(-16*mt**6 + 2*mt**4*s + 16*mt**4*p1p2 - 
     $      2*mt**2*s*p1p2 + 8*mt**2*p1p2*p1p3 - 8*mt**4*p1p4 + 
     $      16*mt**2*p1p2*p1p4 - 8*mt**2*p1p3*p1p4 + 
     $      (128*mt**4*p1p3*p1p4)/s + 8*p1p2*p1p3*p1p4 - 
     $      (128*mt**2*p1p2*p1p3*p1p4)/s - 
     $      (32*p1p2*p1p3**2*p1p4)/s + 
     $      (32*mt**2*p1p3*p1p4**2)/s - 
     $      (64*p1p2*p1p3*p1p4**2)/s - 
     $      (256*mt**2*p1p3**2*p1p4**2)/s**2 + 
     $      (256*p1p2*p1p3**2*p1p4**2)/s**2 - 8*mt**4*p2p3 - 
     $      4*mt**2*p1p4*p2p3 + 
     $      (48*mt**2*p1p3*p1p4*p2p3)/s - 
     $      (16*mt**2*p1p4**2*p2p3)/s + 
     $      (16*p1p3*p1p4**2*p2p3)/s - 
     $      (64*p1p3**2*p1p4**2*p2p3)/s**2 + 
     $      (64*p1p3*p1p4**3*p2p3)/s**2 - 8*mt**4*p2p4 + 
     $      4*mt**2*p1p3*p2p4 - (16*mt**2*p1p3**2*p2p4)/s + 
     $      (48*mt**2*p1p3*p1p4*p2p4)/s - 
     $      (16*p1p3**2*p1p4*p2p4)/s + 
     $      (64*p1p3**3*p1p4*p2p4)/s**2 - 
     $      (64*p1p3**2*p1p4**2*p2p4)/s**2) + 
     $      cu3*(2*mt**4*s + 16*mt**4*p1p2 - 2*mt**2*s*p1p2 - 
     $      16*mt**2*p1p2**2 + 8*mt**2*p1p2*p1p3 - 
     $      8*p1p2**2*p1p4 - 8*mt**2*p1p3*p1p4 + 
     $      8*p1p2*p1p3*p1p4 - 
     $      (64*mt**2*p1p2*p1p3*p1p4)/s + 
     $      (64*p1p2**2*p1p3*p1p4)/s - 
     $      (32*p1p2*p1p3**2*p1p4)/s - 8*mt**4*p2p3 - 
     $      4*mt**2*p1p4*p2p3 - (32*mt**4*p1p4*p2p3)/s + 
     $      (32*mt**2*p1p2*p1p4*p2p3)/s + 
     $      (48*mt**2*p1p3*p1p4*p2p3)/s + 
     $      (32*p1p2*p1p4**2*p2p3)/s + 
     $      (128*mt**2*p1p3*p1p4**2*p2p3)/s**2 + 
     $      (16*p1p3*p1p4**2*p2p3)/s - 
     $      (128*p1p2*p1p3*p1p4**2*p2p3)/s**2 - 
     $      (64*p1p3**2*p1p4**2*p2p3)/s**2 - 
     $      (32*p1p4**3*p2p3**2)/s**2 + 8*mt**2*p1p2*p2p4 + 
     $      4*mt**2*p1p3*p2p4 - (32*mt**4*p1p3*p2p4)/s + 
     $      (32*mt**2*p1p2*p1p3*p2p4)/s - 
     $      (16*mt**2*p1p3**2*p2p4)/s + 
     $      (128*mt**2*p1p3**2*p1p4*p2p4)/s**2 - 
     $      (16*p1p3**2*p1p4*p2p4)/s - 
     $      (128*p1p2*p1p3**2*p1p4*p2p4)/s**2 + 
     $      (64*p1p3**3*p1p4*p2p4)/s**2 - 
     $      (16*mt**2*p1p4*p2p3*p2p4)/s - 
     $      (16*mt**2*p1p3*p2p4**2)/s + 
     $      (32*p1p3**2*p1p4*p2p4**2)/s**2) + 
     $      cu2*(2*mt**4*s - 2*mt**2*s*p1p2 - 16*mt**2*p1p2**2 + 
     $      16*p1p2**3 - 8*p1p2**2*p1p3 - 8*p1p2**2*p1p4 - 
     $      8*mt**2*p1p3*p1p4 + 8*p1p2*p1p3*p1p4 - 
     $      8*mt**4*p2p3 + 16*mt**2*p1p2*p2p3 - 
     $      4*mt**2*p1p4*p2p3 + (64*mt**2*p1p2*p1p4*p2p3)/s - 
     $      (64*p1p2**2*p1p4*p2p3)/s + 
     $      (32*mt**2*p1p3*p1p4*p2p3)/s - 
     $      (32*p1p2*p1p3*p1p4*p2p3)/s + 
     $      (32*p1p2*p1p4**2*p2p3)/s + 
     $      (16*p1p3*p1p4**2*p2p3)/s - 
     $      (16*mt**2*p1p4*p2p3**2)/s - 
     $      (64*mt**2*p1p4**2*p2p3**2)/s**2 + 
     $      (64*p1p2*p1p4**2*p2p3**2)/s**2 + 
     $      (32*p1p3*p1p4**2*p2p3**2)/s**2 - 
     $      (32*p1p4**3*p2p3**2)/s**2 + 8*mt**2*p1p2*p2p4 + 
     $      4*mt**2*p1p3*p2p4 + (64*mt**2*p1p2*p1p3*p2p4)/s - 
     $      (64*p1p2**2*p1p3*p2p4)/s + 
     $      (32*p1p2*p1p3**2*p2p4)/s - 
     $      (16*p1p3**2*p1p4*p2p4)/s - 
     $      (16*mt**2*p1p3*p2p3*p2p4)/s - 
     $      (16*mt**2*p1p4*p2p3*p2p4)/s - 
     $      (128*mt**2*p1p3*p1p4*p2p3*p2p4)/s**2 + 
     $      (128*p1p2*p1p3*p1p4*p2p3*p2p4)/s**2 - 
     $      (16*mt**2*p1p3*p2p4**2)/s - 
     $      (64*mt**2*p1p3**2*p2p4**2)/s**2 + 
     $      (64*p1p2*p1p3**2*p2p4**2)/s**2 - 
     $      (32*p1p3**3*p2p4**2)/s**2 + 
     $      (32*p1p3**2*p1p4*p2p4**2)/s**2)
*
      do i=1,3
         do j=1,68
            mat(i,j)=mat(i,j)/4d0
         end do
      end do
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine matriu(s,t1,t2,u1,u2,mt,mh,matu)
c----------------------------------------------------------------------
c mat(i,j) = (summe ueber polarisat. der gluonen)*trace(i-channel-
c            born-matrixelem.*(p1s-mt)*j-entwicklungskoeff.*(p2+mt))/4
c  i: s,t,u-channel, j: 1-25 u-channel SMEs 
c----------------------------------------------------------------------
      implicit none
      integer i,j
      real*8 s,s1,s2,s3,t1,t2,u1,u2,mt,mh
      real*8 p1p2,p1p3,p1p4,p2p3,p2p4
      real*8 cs1,cs2,ct1,ct2,ct3,cu1,cu2,cu3
      real*8 matu(3,68)
c initialization
      do i=1,3
         do j=1,68
            matu(i,j)=0d0
         end do
      end do
      s3=4d0*mt**2+mh**2-s-t1-t2-u1-u2
      p1p2=(s3-2d0*mt**2)/2d0
      p1p3=(mt**2-t2)/2d0
      p1p4=(mt**2-u2)/2d0
      p2p3=(mt**2-u1)/2d0
      p2p4=(mt**2-t1)/2d0
      s1=s+t2+u2-2d0*mt**2-mh**2
      s2=s+t1+u1-2d0*mt**2-mh**2
      cs1=1d0/(s1-mt**2)/s
      cs2=1d0/(s2-mt**2)/s
      ct1=1d0/(s1-mt**2)/(t2-mt**2)
      ct2=1d0/(s2-mt**2)/(t1-mt**2)
      ct3=1d0/(t1-mt**2)/(t2-mt**2)
      cu1=1d0/(s1-mt**2)/(u2-mt**2)
      cu2=1d0/(s2-mt**2)/(u1-mt**2)
      cu3=1d0/(u1-mt**2)/(u2-mt**2)

      do i=1,3
         do j=1,68
            matu(i,j)=matu(i,j)/4d0
         end do
      end do
      return
      end
