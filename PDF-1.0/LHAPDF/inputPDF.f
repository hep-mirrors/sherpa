      subroutine parmPDF(x,pdf)
      implicit none
      include 'parmsetup.f'
      character*16 s1,s2
      integer i,j,iset,nop,parmN,Nfunc,Fw,nfmax
      real*8 x,N,b0,Poly,pdf(-6:6),Fparm(nopmax),F(nofmax)
      real*8 Ccoef(-6:6,nofmax),Fpow(nofmax)
      real*8 Q,Treshold(-6:6)/13*0d0/
      integer Fmap(nofmax,npfmax)
      integer Ftype(nofmax),Fn(nofmax),Ctype(-6:6)
      save Nfunc,Fn,Fw,Fpow,Fmap,Ccoef,Fparm,Ftype,Ctype
*
      do i=1,Nfunc
         if (Ftype(i).eq.1) then
            Poly=1.0
            do j=4,Fn(i)
               Poly=Poly+Fparm(Fmap(i,j))*x**(float(j-3)/Fpow(i))
            enddo
            Poly=Fparm(Fmap(i,1))*Poly
            F(i)=x**Fparm(Fmap(i,2))*(1.0-x)**Fparm(Fmap(i,3))*Poly
         endif
         if (Ftype(i).eq.2) then
            if (x.lt.0.9999999) then
               Poly=Fparm(Fmap(i,2))*log(x)+Fparm(Fmap(i,3))*log(1.0-x)
     .             +Fparm(Fmap(i,4))*x
     .             +Fparm(Fmap(i,6))*log(1.0+x*exp(Fparm(Fmap(i,5))))
               F(i)=Fparm(Fmap(i,1))*exp(Poly)
            else
               F(i)=0d0
            endif
         endif
         if (Ftype(i).eq.101) then
            Poly=exp(Fparm(Fmap(i,1)))*x**(Fparm(Fmap(i,2))-1)
     .          *(1d0-x)**Fparm(Fmap(i,3))
            Poly=Poly+(1d0+Fparm(Fmap(i,4))*x)*(1d0-x)**Fparm(Fmap(i,5))
            b0=10d0
            if (Poly.gt.b0) then
               F(i)=Poly
            elseif (Poly.lt.-b0) then
               F(i)=0d0
            else
               F(i)=Poly+log(1d0+exp(-b0*Poly)-exp(-b0))/b0
            endif
         endif
      enddo
      do i=-6,6
         pdf(i)=0.0
         if (Ctype(i).gt.0) then
            if (Ctype(i).eq.1) then
               do j=1,Nfunc
                  pdf(i)=pdf(i)+Ccoef(i,j)*F(j)
               enddo
            endif
            if (Ctype(i).eq.101) then
               if (i.eq.-2) then
                  pdf(i)=F(Ccoef(i,1))/(F(Ccoef(i,2))+1d0)
               endif
               if (i.eq.-1) then
                  pdf(i)=F(Ccoef(i,1))/(1d0/F(Ccoef(i,2))+1d0)
               endif
               if (i.eq.1) then
                  pdf(i)=F(1)+pdf(-1)
               endif
               if (i.eq.2) then
                  pdf(i)=F(2)+pdf(-2)
               endif
            endif
         endif
      enddo
      return
*
      entry weightPDF(x)
      if (Fw.ge.0) then
         x=Fparm(Fw)
      else
         call numberPDF(nop)
         x=1.0/float(nop)
      endif
      return
*
      entry GetParmPDF(iset,x)
      x=Fparm(iset)
      return
*
      entry GetNf(nfmax)
      nfmax=0
      do i=1,6
         if (Treshold(-i).ge.0d0) nfmax=nfmax+1
         if (Treshold(i).ge.0d0) nfmax=nfmax+1
      enddo
      nfmax=nfmax/2
      return
*
      entry GetThreshold(iset,Q)
      Q=Treshold(iset)
      return
*
      entry InitPDF(iset)
      call listPDF(iset,Fparm)
      call InitEvolvePDF
      return
*
      entry initInputPDF
      read(1,*) s1,Fw,Nfunc
      write(*,*) 'Parametrization: ',s1
      write(*,*)
      do i=1,Nfunc
         Ftype(i)=-1
         read(1,*) s1,s2
         if (index(s2,'x-taylor').eq.1) then
            Ftype(i)=1
            read(1,*) FPow(i),Fn(i)
            read(1,*) (Fmap(i,j),j=1,Fn(i))
         endif
         if (index(s2,'log-pade').eq.1) then
            Ftype(i)=2
            FPow(i)=0d0
            read(1,*) Fn(i)
            read(1,*) (Fmap(i,j),j=1,Fn(i))
         endif
         if (index(s2,'cteq6-ratio').eq.1) then
            Ftype(i)=101
            Fpow(i)=0d0
            Fn(i)=5
            read(1,*) (Fmap(i,j),j=1,Fn(i))
         endif
         if (Ftype(i).lt.0) then
            write(*,*) 'File description error:'
            write(*,*) 'Unknown functional ',s2
            stop
         endif
      enddo
      read(1,*) s1
      do i=-6,6
         Ctype(i)=-1
         read(1,*) s1,s2
         if (index(s2,'none').eq.1) then
            Ctype(i)=0
            Treshold(i)=-1d0
         endif
         if (index(s2,'treshold').eq.1) then
            Ctype(i)=0
            read(1,*) Treshold(i)
         endif
         if (index(s2,'composite').eq.1) then
            Ctype(i)=1
            Treshold(i)=0d0
            read(1,*) (Ccoef(i,j),j=1,Nfunc)
         endif
         if (index(s2,'cteq6-ratio').eq.1) then
            Ctype(i)=101
            Treshold(i)=0d0
            read(1,*) (Ccoef(i,j),j=1,3)
         endif
         if (Ctype(i).lt.0) then
            write(*,*) 'File description error:'
            write(*,*) 'Unknown composit type ',s2
            stop
         endif
      enddo
      if (Fw.ge.0) then
         write(*,*) '***********************************************'
         write(*,*) '* Note that this is a weigthed PDF set.       *'
         write(*,*) '* See manual for proper use.                  *'
         write(*,*) '***********************************************'
      endif
      return
*
      end
