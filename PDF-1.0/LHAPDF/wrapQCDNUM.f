      subroutine evolvePDF(x,Q,pdf)
      implicit none
      character*64 gridname
      character*16 s1,dummy
      integer Eorder,index
      real*8 x,Q,Qlam,Q2,Q2fit,alfas0,scale0,alfas
      real*8 mc,mc2,mb,mb2,mt,mt2,tc,tc2,tb,tb2
      real*8 singlet,dm,um,dp,up,sp,ub,db,sb,cb,bb,bp,cp,gl,xx
      real*8 QALFAS,QPDFXQ,XFROMIX
      real*8 f(-6:6),pdf(-6:6)
      integer iq0,iqc,iqb,nf,nf0,qnerr,NFLGET,iflag,ix,i,IQFROMQ
      integer nx,nq
      real*8 xmin,xmax,qmin,qmax,S
      save iq0,iqc,iqb,nf0,mc2,mb2,tc2,tb2,mt2
      save nx,xmin,xmax,nq,qmin,qmax,gridname
*
      Q2=Q*Q
      nf=5
      if (Q2.lt.tb2) nf=4
      if (Q2.lt.tc2) nf=3
      pdf(0)=QPDFXQ('GLUON'  ,x,Q2,iflag)
      singlet=QPDFXQ('SINGLET',x,Q2,iflag)
      dm= QPDFXQ('DM',x,Q2,IFLAG)
      um= QPDFXQ('UM',x,Q2,IFLAG)
      dp= QPDFXQ('DP',x,Q2,IFLAG)
      up= QPDFXQ('UP',x,Q2,IFLAG)
      sp= QPDFXQ('SP',x,Q2,IFLAG)
      ub=0.5d0*(up-um+singlet/dble(nf))
      db=0.5d0*(dp-dm+singlet/dble(nf))
      sb=0.5d0*(sp+singlet/dble(nf))
      cb=0.d0
      if (nf.ge.4) then
         cp= QPDFXQ('CP',X,Q2,IFLAG)
         cb=0.5d0*(cp+singlet/dble(nf))
      end if
      bb=0.d0
      if (nf.ge.5) then
         bp= QPDFXQ('BP',X,Q2,IFLAG)            
         bb=0.5d0*(bp+singlet/dble(nf))
      end if

      pdf(1)=dm+db
      pdf(2)=um+ub
      pdf(3)=sb
      pdf(4)=cb
      pdf(5)=bb
      pdf(6)=0d0
      pdf(-1)=db
      pdf(-2)=ub
      pdf(-3)=sb
      pdf(-4)=cb
      pdf(-5)=bb
      pdf(-6)=0d0
      return
*
      entry alfasevolve(alfas,Q)
      Q2=Q*Q
      nf=6
      if (Q2.lt.mt2) nf=5
      if (Q2.lt.mb2) nf=4
      if (Q2.lt.mc2) nf=3
      alfas=QALFAS(Q2,Qlam,nf,iflag)
      return
*
      entry readevolve
      read(1,*) s1
      if (index(s1,'QCDNUM').ne.1) then
         write(*,*) s1,' evolution code requested.'
         write(*,*) 'LHAPDF is compiled with QCDNUM evolution code.'
         write(*,*) 'Using QCDNUM evolution with default large grid.'
         gridname='large.grid'
         nx=400
         xmin=1d-6
         xmax=1d0
         nq=112
         qmin=1d0
         qmax=1d10
         read(1,*) dummy
      else
         read(1,*) gridname,nx,xmin,xmax,nq,qmin,qmax
      endif
      return
*
      entry initevolution(Eorder,Q2fit)
      call qninit
      call QNISET('ORDER',Eorder+1)
      call grxdef(nx,xmin)
      call grqdef(nq,qmin,qmax)
      call grqinp(Q2fit,1)
      call getQMass(4,mc)
      mc2=mc*mc
      CALL QNRSET('CMASS',mc)
      CALL QNRSET('MCALF',mc)
      call getThreshold(4,tc)
      tc2=tc*tc
      call GRQINP(tc2,1)
      CALL GRQINP(tc2-1.0d-4,1)
      call getQMass(5,mb)
      mb2=mb*mb
      CALL QNRSET('BMASS',mb)
      CALL QNRSET('MBALF',mb)
      call getThreshold(5,tb)
      tb2=tb*tb
      call GRQINP(tb2,1)
      CALL GRQINP(tb2-1.0d-4,1)
      call getQMass(6,mt)
      mt2=mt*mt
      CALL QNRSET('TMASS',mt)
      CALL QNRSET('MTALF',mt)
      call qthres(tc2,tb2)
      iq0=IQFROMQ(Q2fit)
      iqc=IQFROMQ(tc2)
      iqb=IQFROMQ(tb2)
      nf0= NFLGET(iq0)
      CALL QNBOOK(2,'dm')
      CALL QNBOOK(3,'um')
      CALL QNBOOK(4,'dp')
      CALL QNBOOK(5,'up')
      CALL QNBOOK(6,'sp')
      CALL QNBOOK(7,'cp')
      CALL QNBOOK(8,'bp')
      CALL QNLSET('W1ANA',.TRUE.)
      CALL QNLSET('W2NUM',.TRUE.)
      CALL QNLSET('W2STF',.FALSE.)
      if (index(gridname,'none').eq.1) then
         call qnfilw(0,0)
      else
         qnerr=-1
         open(unit=2,status='old',file=gridname,
     .        form='unformatted',err=1)
         call QNREAD(2,1,qnerr)
 1       close(2)
         if (qnerr.ne.0) then
            write(*,*) 'Grid file problem: ',gridname
            if (qnerr.lt.0) then 
               write(*,*) 'Grid file does not exist'
               write(*,*) 'Calculating and creating grid file'
               call qnfilw(0,0)
               open(unit=2,status='unknown',file=gridname,
     .              form='unformatted')
               call QNDUMP(2)
               close(2)
            else
               write(*,*) 'Existing grid file is inconsistent'
               if (qnerr.eq.1)
     .              write(*,*) 'Defined grid different'
               if (qnerr.eq.2)
     .              write(*,*) 'Heavy quark weight table different'
               if (qnerr.eq.3)
     .              write(*,*) 'Charm mass different'
               if (qnerr.eq.4)
     .              write(*,*) 'Bottom mass different'
               stop
            endif
         endif
      endif
      return
*
      entry initevolvePDF
      call GetAlfas(alfas0,scale0)
      Q2=scale0*scale0
      CALL QNRSET('ALFAS',alfas0)      
      CALL QNRSET('ALFQ0',Q2)
      DO ix = 1,nx
         xx = XFROMIX(ix)
         call parmPDF(xx,f)
         singlet=0.d0
         do i=1,nf0
            singlet=singlet+f(i)+f(-i)
         end do
         gl=f(0)
         dm=f(1)-f(-1)
         um=f(2)-f(-2)
         dp=f(1)+f(-1)-singlet/dble(nf0)
         up=f(2)+f(-2)-singlet/dble(nf0)
         sp=f(3)+f(-3)-singlet/dble(nf0)
         CALL QNPSET('SINGLET',ix,iq0,singlet)
         CALL QNPSET('GLUON',ix,iq0,gl)
         CALL QNPSET('DM',ix,iq0,DM)
         CALL QNPSET('UM',ix,iq0,UM)
         CALL QNPSET('DP',ix,iq0,DP)
         CALL QNPSET('UP',ix,iq0,UP)
         CALL QNPSET('SP',ix,iq0,SP)
      ENDDO
      CALL EVOLSG(iq0,1,nq)
      CALL EVOLNM('DM',iq0,1,nq)
      CALL EVOLNM('UM',iq0,1,nq)
      CALL EVLSEA('dp',iq0,iqc,iqb,nq)
      CALL EVLSEA('up',iq0,iqc,iqb,nq)
      CALL EVLSEA('sp',iq0,iqc,iqb,nq)

C--   Heavy quark evolution

      CALL EVOLCP('cp',iqc,iqb,nq)
      CALL EVOLBP('bp',iqb,nq)

      return
*
      end
      subroutine EVLSEA(name,IQ0,IQC,IQB,NQGRI)
      implicit none
      
      CHARACTER*(*) name
      integer iq0,iqc,iqb,nqgri
      real*8 f34,f45,f43,f54,factor
      parameter(f34=1.d0/12.d0,f45=1.d0/20.d0,f43=-1.d0/12.d0,
     >          f54=-1.d0/20.d0)


      If(IQ0.le.IQC)then
         CAll EVPLUS(name,IQ0,1,IQC)
         CALL QADDSI(name,IQC,f34)
         CAll EVPLUS(name,IQC,IQC,IQB)
         CALL QADDSI(name,IQB,f45)
         CAll EVPLUS(name,IQB,IQB,NQGRI)
      else if(IQ0.gt.IQC.and.IQ0.le.IQB)then
         CAll EVPLUS(name,IQ0,IQC,IQB)
         CALL QADDSI(name,IQC,f43)
         CAll EVPLUS(name,IQC,1,IQC)
         CALL QADDSI(name,IQC,f34)
         CALL QADDSI(name,IQB,f45)
         CAll EVPLUS(name,IQB,IQB,NQGRI)
      else if(IQ0.gt.IQB)then
         CAll EVPLUS(name,IQ0,IQB,NQGRI)
         CALL QADDSI(name,IQB,f54)
         CAll EVPLUS(name,IQB,IQC,IQB)
         CALL QADDSI(name,IQB,f45)
         CALL QADDSI(name,IQC,f43)
         CAll EVPLUS(name,IQC,1,IQC)
         CALL QADDSI(name,IQC,f34)
      end if         

      end


      subroutine EVOLCP(name,IQC,IQB,NQGRI)
      implicit none
      
      CHARACTER*(*) name
      integer iq0,iqc,iqb,nqgri
      real*8 f4,f45
      parameter(f4=-1.d0/4.d0,f45=1.d0/20.d0)

c
c     First set to zero to avoid adding -1/4Singl at each iteration
c

      CAll QNPNUL(name)
      CALL QADDSI(name,IQC,f4)      
      CAll EVPLUS(name,IQC,IQC,IQB)
      CALL QADDSI(name,IQB,f45)
      CAll EVPLUS(name,IQB,IQB,NQGRI)
      
      end

      subroutine EVOLBP(name,IQB,NQGRI)
      implicit none
      
      CHARACTER*(*) name
      integer iqb,nqgri
      real*8 f5
      parameter(f5=-1.d0/5.d0)

c
c     First set to zero to avoid adding -1/5Singl at each iteration
c

      CAll QNPNUL(name)
      CALL QADDSI(name,IQB,f5)      
      CAll EVPLUS(name,IQB,IQB,NQGRI)
      
      end

