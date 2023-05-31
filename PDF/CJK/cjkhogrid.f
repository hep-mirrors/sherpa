*********************************************************************
****   NLO parametrization of the parton densities in the real   ****
****            photon for the CJK NLO fit based on              ****
****    "A New 5 Flavour NLO Analysis and Parametrizations       ****
****       of Parton Distributions of the Real Photon"           ****
****           by F.Cornet, P.Jankowski & M.Krawczyk             ****
****                     hep-ph/0404063                          ****
*********************************************************************
****                                                             ****
****   valid for 10^(-5) < x < 1 and 1 < Q^2 < 2*10^5 GeV^2      ****
****                                                             ****
****       x      - Bjorken x variable                           ****
****       Q2     - square of momentum scale (in GeV**2)         ****
****   XPDF(-5:5) - matrix containing x*f(x,Q2)/alfa             ****
****                                                             ****
****   PDF =   -5,  -4,    -3 ,  -2 ,  -1 ,0 ,1,2,3,4,5          ****
****         b_bar,c_bar,s_bar,u_bar,d_bar,gl,d,u,s,c,b          ****
****                                                             ****
****  All antiquark and corresponding quark densities are equal. ****
****                                                             ****
****   heavy-quark masses:                                       ****
****             M_charm=1.3, M_beauty=4.3 GeV                   ****
****                                                             ****
****   Lambda_QCD values for active N_f falvours:                ****
****   N_f =   3      4      5                                   ****
****         0.323  0.280  0.200 GeV                             ****
****                                                             ****
****  Grid parametrization utilizing the bicubic interpolation   ****
****  in the Hermite polynomials basis.                          ****
****                                                             ****
****  To use it one must add in ones main program:               ****
****        INTEGER IREAD                                        ****
****        common /IREADVH2/ IREAD                              ****
****        IREAD = 0                                            ****
****  This allows for fast multiple use of the parametrization.  ****
****                                                             ****
****                                                             ****
****   IOPT = 1 : light partons (gl,up,dn,str)                   ****
****        = 2 : all partons   (gl,up,dn,str,chm,bot)           ****
****                                                             ****
*********************************************************************
****  Evolution, parametrization, checks performed and programs  ****
****  written by                                                 ****
****     Pawel Jankowski, (pjank@fuw.edu.pl)                     ****
****                      Institute of Theoretical Physics,      ****
****                      Warsaw University                      ****
****                      ul. Hoza 69                            ****
****                      00-681 Warsaw                          ****
****                      Poland                                 ****
****                                                             ****
****  Last changes - 07 April 2004                               ****
*********************************************************************

      SUBROUTINE CJKHOGRID(IOPT,X,Q2,XPDF)
      IMPLICIT DOUBLE PRECISION (a-z)
      INTEGER IOPT
      parameter (Pi=3.1415926535897932d0,alfa=7.29735308D-3,
     %           e1=1.d0/3.d0,e2=2.d0/3.d0)
      logical t
      dimension XPDF(-5:5)

      XX = X
      QQ2 = Q2

      if ((XX.LE.1.d-5).OR.(XX.GE.1.d0)) then
         print *,'X out of range: ',XX
         stop
      endif
      if ((QQ2.LE.5.d-1).OR.(QQ2.GE.5.d5)) then
         print *,'Q2 out of range: ',QQ2
         stop
      endif

      if (IOPT.EQ.1) then

         call GRIDVH2(1,XX,QQ2,XGLU,XDN,XUP,XSTR,XCHM,XBOT)

         XPDF(0) = xglu/alfa
         XPDF(1) = xdn/alfa
         XPDF(2) = xup/alfa
         XPDF(3) = xstr/alfa
         XPDF(4) = 0.d0
         XPDF(5) = 0.d0

         goto 1000

      endif

      if (IOPT.EQ.2) then

         call GRIDVH2(2,XX,QQ2,XGLU,XDN,XUP,XSTR,XCHM,XBOT)

         XPDF(0) = xglu/alfa
         XPDF(1) = xdn/alfa
         XPDF(2) = xup/alfa
         XPDF(3) = xstr/alfa
         XPDF(4) = xchm/alfa
         XPDF(5) = xbot/alfa

      endif

 1000 continue

      do 10 I=1,5
 10      XPDF(-I) = XPDF(I)

      RETURN
      END

*********************************************************************
*********************************************************************

      SUBROUTINE GRIDVH2(IOPT,XIN,Q2IN,XGLU,XDN,XUP,XSTR,XCHM,XBOT)
      IMPLICIT DOUBLE PRECISION (a-z)
      INTEGER I,J,JH,K,L,M,N,bord,bordc,bordb,Imaxc,Imaxb,
     %        NX,NQ2,NQ2H,NQQ2,NQQ2H,IOPT,IREAD
      parameter (NX=52,NQ2=32,NQ2H=48,alfa=7.29735308D-3)
      dimension gl(0:NQ2+1,0:NX+1),dn(0:NQ2+1,0:NX+1),
     %          up(0:NQ2+1,0:NX+1),st(0:NQ2+1,0:NX+1),
     %          ch(0:NQ2H+1,0:NX+1),bt(0:NQ2H+1,0:NX+1),
     %          chx(0:NQ2H+1,24),btx(0:NQ2H+1,24),
     %          x(4),xc(4),xb(4),q2(4),q2h(4),
     %          glh(4,4),dnh(4,4),uph(4,4),sth(4,4),
     %          chh(4,4),bth(4,4),Imaxc(NQ2H),Imaxb(NQ2H),
     %          xdata(0:NX+1),q2data(0:NQ2+1),q2hdata(0:NQ2H+1)
      data xdata/0d0,1d-5,2d-5,4d-5,6d-5,8d-5,
     %	         1d-4,2d-4,4d-4,6d-4,8d-4,
     %	         1d-3,2d-3,4d-3,6d-3,8d-3,
     %	         1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     %	        .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     % 	        .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     %	        .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     %	        .8d0,.85d0,.9d0,.95,.98,1d0,0d0/
      data q2data/0d0,0.5d0,0.75d0,1.d0,1.25d0,1.5d0,2d0,2.5d0,3.2d0,
     %            4d0,5d0,6.4d0,8d0,1d1,1.2d1,1.8d1,2.6d1,4d1,6.4d1,
     %            1d2,1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     %            1.8d4,3.2d4,5.6d4,1d5,0d0/
      data q2hdata/0d0,0.5d0,0.75d0,1.d0,1.25d0,1.5d0,1.75d0,2d0,2.25,
     %             2.5d0,2.75d0,3.2d0,3.6d0,4d0,4.5d0,5d0,5.4d0,6d0,
     %             6.4d0,7.2d0,8d0,1d1,1.2d1,1.5d1,1.8d1,2.2d1,2.6d1,
     %             3.3d1,4d1,5.2d1,6.4d1,7.2d1,8.5d1,1d2,1.3d2,1.6d2,
     %             2d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,1.8d4,
     %             3.2d4,5.6d4,1d5,0d0/
      data Imaxc/21,22,23,24,25,26,27,28,29,29,30,31,32,33,34,35,36,37,
     %           38,39,41,42,43,44,45,45,46,46,47,47,48,48,48,48,48,48,
     %           49,49,49,49,49,49,49,49,49,49,49,49/
      data Imaxb/14,15,16,17,17,18,18,18,19,19,19,20,20,20,20,21,21,21,
     %           22,22,23,24,25,26,27,28,30,32,34,36,37,38,40,42,43,44,
     %           45,46,47,48,48,49,49,49,49,49,49,49/
      common /IREADVH2/ IREAD
      common /PARTVH2/ gl,dn,up,st,ch,bt,chx,btx

      XX = XIN
      QQ2 = Q2IN

      mc = 1.3d0
      mc2 = mc*mc
      mb = 4.3d0
      mb2 = mb*mb

      xmaxc = 1.d0/(1.d0+4.d0*mc2/QQ2)
      xmaxb = 1.d0/(1.d0+4.d0*mb2/QQ2)

      xchm = 0.d0
      xbot = 0.d0

      bord = 1
      bordc = 1
      bordb = 1

*****************************************************************
****                  reading grid data                      ****
*****************************************************************

      if (IREAD.EQ.0) then
         call readtabvh2(IOPT)
      endif
      IREAD = 100

*****************************************************************
****   searching for the J such that: Q2(J) < QQ2 < Q2(J+1)  ****
****    searching for the I such that: x(I) < XX < x(I+1)    ****
****            for the light quarks and gluon               ****
*****************************************************************

      NQQ2 = NQ2
      call findq2vh2(NQQ2,q2data,QQ2,J)
      call findxvh2(xdata,XX,I)
      if (I.EQ.1.OR.I.EQ.NX-1) bord = 0
      if (J.EQ.1.OR.J.EQ.NQ2-1) bord = 0

*****************************************************************
****     *x1(I-1)   $x2(I)   $$xx   $x3(I+1)   *x4(I+2)      ****
*****************************************************************
****              Only 3 points at borders!!!                ****
*****************************************************************

      do 10 K=1,4
         x(K)  = xdata(I+K-2)
 10      q2(K) = q2data(J+K-2)

*****************************************************************
****                  tbh(1,1) = tb(J-1,I-1)                 ****
****                  tbh(1,2) = tb(J-1,I)                   ****
****                  tbh(1,3) = tb(J-1,I+1)                 ****
****                  tbh(1,4) = tb(J-1,I+2) ...             ****
*****************************************************************

      do 21 K=J-1,J+2
         do 21 L=I-1,I+2
            M = K+2-J
            N = L+2-I
            glh(M,N) = gl(K,L)
            dnh(M,N) = dn(K,L)
            uph(M,N) = up(K,L)
 21         sth(M,N) = st(K,L)
      call fitvh2(bord,QQ2,XX,q2,x,glh,xglu)
      call fitvh2(bord,QQ2,XX,q2,x,dnh,xdn)
      call fitvh2(bord,QQ2,XX,q2,x,uph,xup)
      call fitvh2(bord,QQ2,XX,q2,x,sth,xstr)
      xglu = alfa*xglu
      xdn  = alfa*xdn
      xup  = alfa*xup
      xstr = alfa*xstr

      if (IOPT.EQ.1) goto 1000

*****************************************************************
****   searching for the J such that: Q2(J) < QQ2 < Q2(J+1)  ****
****    searching for the I such that: x(I) < XX < x(I+1)    ****
****                  for the heavy quarks                   ****
*****************************************************************

      NQQ2H = NQ2H
      call findq2vh2(NQQ2H,q2hdata,QQ2,JH)

      if (JH.EQ.1.OR.JH.EQ.NQ2H-1) then
         bordc = 0
         bordb = 0
      endif

      do 30 K=1,4
 30      q2h(K) = q2hdata(JH+K-2)

      if (XX.LT.xmaxc) then
         xmxc = 1.d0/(1.d0+4.d0*mc2/q2h(2))
         call findxhvh2(xmxc,xdata,ch,chx,XX,Imaxc(JH),JH,bordc,xc,chh)
         call fitvh2(bordc,QQ2,XX,q2h,xc,chh,xchm)
         xchm = alfa*xchm
      endif

      if (XX.LT.xmaxb) then
         xmxb = 1.d0/(1.d0+4.d0*mb2/q2h(2))
         call findxhvh2(xmxb,xdata,bt,btx,XX,Imaxb(JH),JH,bordb,xb,bth)
         call fitvh2(bordb,QQ2,XX,q2h,xb,bth,xbot)
         xbot = alfa*xbot
      endif

 1000 continue

      RETURN
      END

*****************************************************************

      SUBROUTINE READTABVH2(IOPT)
      DOUBLE PRECISION gl,dn,up,st,ch,bt,chx,btx
      parameter (NX=52,NQ2=32,NQ2H=48)
      dimension gl(0:NQ2+1,0:NX+1),dn(0:NQ2+1,0:NX+1),
     %          up(0:NQ2+1,0:NX+1),st(0:NQ2+1,0:NX+1),
     %          ch(0:NQ2H+1,0:NX+1),bt(0:NQ2H+1,0:NX+1),
     %          chx(0:NQ2H+1,24),btx(0:NQ2H+1,24)
      character*15 name
      common /PARTVH2/ gl,dn,up,st,ch,bt,chx,btx

      open(10,file='cjkhobest.dat',status='old')

      do 1 I=0,NQ2+1
         do 1 J=0,NX+1
            gl(I,J) = 0.d0
            dn(I,J) = 0.d0
            up(I,J) = 0.d0
 1          st(I,J) = 0.d0

         do 10 I=1,NQ2
            do 10 J=1,NX
 10            read(10,100) gl(I,J),dn(I,J),
     %                      up(I,J),st(I,J)
 100          format (F11.7,3(2X,F10.7))

         if (IOPT.EQ.1) goto 1000

         do 2 I=0,NQ2H+1
            do 2 J=0,NX+1
               ch(I,J) = 0.d0
 2             bt(I,J) = 0.d0

         do 20 I=1,NQ2H
            do 20 J=1,NX
 20            read(10,200) ch(I,J),bt(I,J)
 200           format (F10.7,2X,F10.7)

         do 30 I=1,NQ2H
            read(10,300) (chx(I,J),J=1,24)
 30      continue
         do 40 I=1,NQ2H
            read(10,300) (btx(I,J),J=1,24)
 40      continue
 300     format (23(F10.7,2X),F10.7)
         do 31 I=0,NQ2H+1,NQ2H+1
            do 31 J=1,24
               chx(I,J) = 0.d0
 31            btx(I,J) = 0.d0

 1000    continue

         close(10)

 2000 continue

      RETURN
      END

*********************************************************************
****            Here I use the bisection method                  ****
*********************************************************************

      SUBROUTINE FINDQ2VH2(NQ2,q2data,QQ2,I)
      DOUBLE PRECISION QQ2,Q2,q2data
      INTEGER I,iu,ul,NQ2
      dimension q2data(0:NQ2+1)

      il = 1
      iu = NQ2

 10   if (iu-il.GT.1) then
         I = (iu+il)/2
         Q2 = q2data(I)
         if (QQ2.GE.Q2) then
            il = I
         else
            iu = I
         endif
      goto 10
      endif
      I = il

 100  continue

      RETURN
      END

*********************************************************************
****            Here I use the bisection method                  ****
*********************************************************************

      SUBROUTINE FINDXVH2(xdata,XX,I)
      DOUBLE PRECISION XX,x,xdata
      INTEGER I,il,iu
      parameter (NX=52)
      dimension xdata(0:NX+1)

      il = 1
      iu = NX

 10   if (iu-il.GT.1) then
         I = (iu+il)/2
         x = xdata(I)
         if (XX.GE.x) then
            il = I
         else
            iu = I
         endif
      goto 10
      endif
      I = il

 100  continue

      RETURN
      END

*********************************************************************
****            Here I use the bisection method                  ****
*********************************************************************

      SUBROUTINE FINDXHVH2(xmax,xdata,hq,hqx,XX,Imax,J,bord,x,hqh)
      IMPLICIT DOUBLE PRECISION (a-z)
      INTEGER I,J,K,L,M,N,il,iu,Imax,bord,NX,NQ2
      parameter (NX=52,NQ2=48)
      dimension xdata(0:NX+1),hq(0:NQ2+1,0:NX+1),
     % hqx(0:NQ2+1,24),hqh(4,4),x(4),per(7),xdatah(0:Imax+7)
      data per/0.93d0,0.95d0,0.96d0,0.97d0,0.99d0,0.999d0,0d0/

      do 1 I=0,Imax
 1       xdatah(I) = xdata(I)
      do 2 I=1,7
 2       xdatah(Imax+I) = per(I)*xmax

      il = 1
      iu = Imax+6

 10   if (iu-il.GT.1) then
         I = (iu+il)/2
         xb = xdatah(I)
         if (XX.GE.xb) then
            il = I
         else
            iu = I
         endif
         goto 10
       endif
       I = il

       K = I-Imax+1
       if (K.LT.0) K = 0

*****************************************************************
****     *x1(I-1)   $x2(I)   $$xx   $x3(I+1)   *x4(I+2)      ****
*****************************************************************
****              Only 3 points at borders!!!                ****
*****************************************************************

      do 20 L=1,4
 20      x(L) = xdatah(I+L-2)


      if (K.EQ.0) then

         if (I+2.LE.Imax) then
            do 30 M=J-1,J+2
               do 30 L=I-1,I+2
 30               hqh(M+2-J,L+2-I) = hq(M,L)
         else
            do 40 M=J-1,J+2
               N = M+2-J
               hqh(N,4) = hqx(M,6*N-5)
               do 40 L=I-1,I+1
 40               hqh(N,L+2-I) = hq(M,L)
         endif

      elseif (K.EQ.1) then

         do 50 M=J-1,J+2
            N = M+2-J
            hqh(N,1) = hq(M,Imax-1)
            hqh(N,2) = hq(M,Imax)
            hqh(N,3) = hqx(M,6*N-5)
 50         hqh(N,4) = hqx(M,6*N-4)

      elseif (K.EQ.2) then

         do 60 M=J-1,J+2
            N = M+2-J
            hqh(N,1) = hq(M,Imax)
            do 60 I=2,4
 60            hqh(N,I) = hqx(M,6*N+I-7)

         if (hqh(1,3).LT.1.d-10) bord = 0

      elseif (K.EQ.3) then

         do 70 M=J-1,J+2
            N = M+2-J
            do 70 I=1,4
 70            hqh(N,I) = hqx(M,6*N+I-6)

         if (hqh(1,3).LT.1.d-10.OR.hqh(2,4).LT.1.d-10) bord = 0

      elseif (K.EQ.4) then

         do 80 M=J-1,J+2
            N = M+2-J
            do 80 I=1,4
 80            hqh(N,I) = hqx(M,6*N+I-5)

         if (hqh(1,3).LT.1.d-10.OR.hqh(2,4).LT.1.d-10) bord = 0

      elseif (K.EQ.5) then

         do 90 M=J-1,J+2
            N = M+2-J
            do 90 I=1,4
 90            hqh(N,I) = hqx(M,6*N+I-4)

         if (hqh(1,3).LT.1.d-10.OR.hqh(2,4).LT.1.d-10) bord = 0

      elseif (K.EQ.6) then

         bord = 0
         do 100 M=J-1,J+2
            N = M+2-J
            hqh(N,4) = 0.d0
            do 100 I=1,3
 100           hqh(N,I) = hqx(M,6*N+I-3)

      endif

      RETURN
      END

*********************************************************************
****         Here I use the bicubic interpolation in             ****
****              the Hermit polynomials basis                   ****
*********************************************************************

      SUBROUTINE FITVH2(bord,xx1,xx2,x1,x2,yg,result)
      IMPLICIT DOUBLE PRECISION (a-z)
      dimension x1(4),x2(4),yg(4,4),y(4),y1(4),y2(4),y12(4)
      integer bord
      external d1fvh2,d2fvh2

*****************************************************************
****              4 *(x1l,x2u)    3 *(x1u,x2u)               ****
****              1 *(x1l,x2l)    2 *(x1u,x2l)               ****
*****************************************************************

      x1l = x1(2)
      x1u = x1(3)
      x2l = x2(2)
      x2u = x2(3)

*****************************************************************
****      Function values, first and cross-derivatives       ****
****              at 4 corners of a grid cell                ****
*****************************************************************

      y(1) = yg(2,2)
      y(2) = yg(3,2)
      y(3) = yg(3,3)
      y(4) = yg(2,3)

      if (bord.EQ.1) then

         y1(1) = d1fvh2(x1(1),x1(3),yg(1,2),yg(3,2))
         y1(2) = d1fvh2(x1(2),x1(4),yg(2,2),yg(4,2))
         y1(3) = d1fvh2(x1(2),x1(4),yg(2,3),yg(4,3))
         y1(4) = d1fvh2(x1(1),x1(3),yg(1,3),yg(3,3))

         y2(1) = d1fvh2(x2(1),x2(3),yg(2,1),yg(2,3))
         y2(2) = d1fvh2(x2(1),x2(3),yg(3,1),yg(3,3))
         y2(3) = d1fvh2(x2(2),x2(4),yg(3,2),yg(3,4))
         y2(4) = d1fvh2(x2(2),x2(4),yg(2,2),yg(2,4))

         y12(1) = d2fvh2(x1(1),x1(3),x2(1),x2(3),
     %                yg(1,1),yg(1,3),yg(3,1),yg(3,3))
         y12(2) = d2fvh2(x1(2),x1(4),x2(1),x2(3),
     %                yg(2,1),yg(2,3),yg(4,1),yg(4,3))
         y12(3) = d2fvh2(x1(2),x1(4),x2(2),x2(4),
     %                yg(2,2),yg(2,4),yg(4,2),yg(4,4))
         y12(4) = d2fvh2(x1(1),x1(3),x2(2),x2(4),
     %                yg(1,2),yg(1,4),yg(3,2),yg(3,4))

      else

         y1(1) = d1fvh2(x1(2),x1(3),yg(2,2),yg(3,2))
         y1(2) = y1(1)
         y1(3) = d1fvh2(x1(2),x1(3),yg(2,3),yg(3,3))
         y1(4) = y1(3)

         y2(1) = d1fvh2(x2(2),x2(3),yg(2,2),yg(2,3))
         y2(2) = d1fvh2(x2(2),x2(3),yg(3,2),yg(3,3))
         y2(3) = y2(2)
         y2(4) = y2(1)

         y12(1) = d2fvh2(x1(2),x1(3),x2(2),x2(3),
     %                yg(2,2),yg(2,3),yg(3,2),yg(3,3))
         y12(2) = y12(1)
         y12(3) = y12(1)
         y12(4) = y12(1)

      endif

      call itervh2(y,y1,y2,y12,x1l,x1u,x2l,x2u,xx1,xx2,result)

      RETURN
      END

*********************************************************************

      SUBROUTINE itervh2(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,res)
      IMPLICIT DOUBLE PRECISION (a-z)
      DIMENSION y(4),y1(4),y2(4),y12(4),p(4,4),c(4,4),ht(4),hu(4),ph(4)
      INTEGER I,J

      d1 = x1u-x1l
      d2 = x2u-x2l
      d1d2 = d1*d2

****  local variables  ****
      t = (x1-x1l)/d1
      u = (x2-x2l)/d2

****  local derivatives  ****
      p(1,1) = y(1)
      p(2,1) = y(2)
      p(2,2) = y(3)
      p(1,2) = y(4)
      p(3,1) = y1(1)*d1
      p(4,1) = y1(2)*d1
      p(4,2) = y1(3)*d1
      p(3,2) = y1(4)*d1
      p(1,3) = y2(1)*d2
      p(2,3) = y2(2)*d2
      p(2,4) = y2(3)*d2
      p(1,4) = y2(4)*d2
      p(3,3) = y12(1)*d1d2
      p(4,3) = y12(2)*d1d2
      p(4,4) = y12(3)*d1d2
      p(3,4) = y12(4)*d1d2

****  Hermite polynomials  ****
      t1 = t-1.
      ht(1) = (2.*t+1.)*t1*t1
      ht(2) = t*t*(-2.*t+3.)
      ht(3) = t*t1*t1
      ht(4) = t*t*t1

      u1 = u-1
      hu(1) = (2.*u+1.)*u1*u1
      hu(2) = u*u*(-2.*u+3.)
      hu(3) = u*u1*u1
      hu(4) = u*u*u1

      do 10 I=1,4
         ph(I) = 0.d0
         do 10 J=1,4
 10         ph(I) = ph(I) + p(I,J)*hu(J)

      res = 0.d0
      do 20 I=1,4
 20      res = res + ht(I)*ph(I)

      RETURN
      END

*********************************************************************
****                   First derivative: df/dx                   ****
*********************************************************************
      DOUBLE PRECISION FUNCTION d1fvh2(x1,x2,f1,f2)
      IMPLICIT DOUBLE PRECISION (a-z)

      d1fvh2 = (f2-f1)/(x2-x1)

      END

*********************************************************************
****                Second derivative: d2f/dxdy                  ****
*********************************************************************
      DOUBLE PRECISION FUNCTION d2fvh2(x1,x2,y1,y2,f11,f12,f21,f22)
      IMPLICIT DOUBLE PRECISION (a-z)

      d2fvh2 = (f22 - f21 - f12 + f11)/((x2-x1)*(y2-y1))

      END
