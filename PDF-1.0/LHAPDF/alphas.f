      real*8 function alphasPDF(Q)
      implicit none
      real*8 Q,a
*
      call evolveAs(Q,a)
      alphasPDF=a
      return
*
      end
*
      subroutine evolveAs(Q,alphas)
      implicit none
      character*16 s1,s2,s3
      real*8 Q,Q0,alphas,alfas0,alfasQ,scale0,Mc,Mb,Mt
      real*8 b0,b1,b2,L,As,alphasPDF,mass
      parameter (b0=1.2202,b1=0.4897,b2=0.1913)
      integer order,EvlOrd,parm,Etype,Method,nf,n
      save EvlOrd,parm,Etype,alfasQ,Q0,Method,Mc,Mb,Mt
*
      if (method.eq.0) then
         if ((Etype.eq.1).or.(Etype.eq.2)) then
            L=log(Q/Q0)
            As=alfasQ
            if (Etype.eq.2) call GetParmPDF(parm,As)
            if (EvlOrd.eq.0) L=b0*L
            if (EvlOrd.eq.1) L=(b0+As*b1)*L
            if (EvlOrd.eq.2) L=(b0+As*b1+As**2*b2)*L
     .                        -0.5*As**2*b0*b1/2d0*L**2
            alphas=As/(1.0+As*L)
         endif
      endif
      if (method.eq.1) then
         call alfasevolve(alphas,Q)
      endif
      return
*
      entry GetQmass(nf,mass)
      n=abs(nf)
      mass=0d0
      if (n.eq.4) mass=Mc
      if (n.eq.5) mass=Mb
      if (n.eq.6) mass=Mt
      return

      entry GetAlfas(alfas0,scale0)
      scale0=Q0
      alfas0=alfasQ
      if (Etype.eq.2) call GetParmPDF(parm,alfas0)
      return
*
      entry GetOrderAs(order)
      order=EvlOrd
      return
*
      entry InitAlphasPDF
      Etype=-1
      EvlOrd=-1
      read(1,*) s1,s2,s3
      if (index(s2,'lo').eq.1) EvlOrd=0
      if (index(s2,'nlo').eq.1) EvlOrd=1
      if (index(s2,'nnlo').eq.1) EvlOrd=2
      if (EvlOrd.lt.0) then
         write(*,*) 'File description error:'
         write(*,*) 'Unknown alpha_s evolution order ',s2
         stop
      endif
      if (index(s1,'Fixed').eq.1) then
         Etype=1
         parm=-1
         read(1,*) alfasQ,Q0,Mc,Mb,Mt
      endif
      if (index(s1,'Variable').eq.1) then
         Etype=2
         alfasQ=0d0
         read(1,*) parm,Q0,Mc,Mb,Mt
      endif
      if (Etype.lt.0) then
         write(*,*) 'File description error:'
         write(*,*) 'Unknown alpha_S evolution method ',s1
         stop
      endif
      Method=-1
      if (index(s3,'Internal').eq.1) Method=0
      if (index(s3,'EvolCode').eq.1) Method=1
      if (Method.lt.0) then
         write(*,*) 'File description error:'
         write(*,*) 'Unknown alpha_S method ',s3
         stop
      endif
      return
*
      end

