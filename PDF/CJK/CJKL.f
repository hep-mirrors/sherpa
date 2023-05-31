*********************************************************************
****   LO parametrization of the parton densities in the real    ****
****             photon for the CJKL fit based on                ****
****                                                             ****
****        F.Cornet, P.Jankowski, M.Krawczyk and A.Lorca        ****
****      "A New 5 Flavour LO Analysis and Parametrization       ****
****         of Parton Distributions in the Real Photon"         ****
****               Phys. Rev. D68: 014010, 2003                  ****
****                      hep-ph/0212160                         ****
*********************************************************************
****                                                             ****
****   valid for 10^(-5) < x < 1 and 1 < Q^2 < 2*10^5 GeV^2      ****
****                                                             ****
****       x      - Bjorken x variable                           ****
****       Q2     - square of momentum scale (in GeV**2)         ****
****   XPDF(-5:5) - matrix containing x*f(x,Q2)/alfa             ****
****                                                             ****
****   PDF =   -5 ,  -4 ,  -3 ,  -2 ,  -1 ,0 ,1,2,3,4,5          ****
****         b_bar,c_bar,s_bar,u_bar,d_bar,gl,d,u,s,c,b          ****
****                                                             ****
****   heavy quarks masses:                                      ****
****             M_charm=1.3, M_beauty=4.3 GeV                   ****
****                                                             ****
****   Lambda_QCD values for active N_f falvours:                ****
****   N_f=   3      4      5                                    ****
****        0.314  0.280  0.221 GeV                              ****
****                                                             ****
****  LIPARTS: only light partons (gl,d,u,s)                     ****
****  PARTONS: (gl,d,u,s,c,b) parton densities                   ****
****  PARTF2:  (gl,d,u,s,c,b) parton densities                   ****
****           & structure function F2 according to Eq. (22)     ****
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
****  Last changes - 09 DEC 2002                                 ****
*********************************************************************

      SUBROUTINE LIPARTS(x,Q2,XPDF)
      IMPLICIT DOUBLE PRECISION (a-z)
      dimension XPDF(-3:3),PART(-5:5)
      integer I

****  x * par / alfa  ****
      call PARAM(1,x,Q2,PART,F2)
      do 1 I=-3,3
 1       XPDF(I)=PART(I)

      END

      SUBROUTINE PARTONS(x,Q2,XPDF)
      IMPLICIT DOUBLE PRECISION (a-z)
      dimension XPDF(-5:5)

****  x * par / alfa  ****
      call PARAM(2,x,Q2,XPDF,F2)

      END

      SUBROUTINE PARTF2(x,Q2,XPDF,F2)
      IMPLICIT DOUBLE PRECISION (a-z)
      dimension XPDF(-5:5)

****  x * par / alfa  ****
****  F2 / alfa       ****
      call PARAM(3,x,Q2,XPDF,F2)

      END

**********************************

      SUBROUTINE PARAM(OPT,x,Q2,XPDF,F2)
      IMPLICIT DOUBLE PRECISION (a-z)
      Parameter (Pi=3.1415926535897932d0, e1=1.d0/3.d0, e2=2.d0/3.d0,
     % alfa=7.29735308D-3)
      dimension XPDF(-5:5),res(2)
      logical t
      integer OPT,I,step,flav
      common /flav/ flav
      common /mass/ mc,mb

      mc = 1.3d0
      mb = 4.3d0
      Lam3 = 0.314d0

      call POINTLIKE(OPT,x,Q2,glpl,uppl,dnpl,chpl,btpl)
      call HADRONLIKE(OPT,x,Q2,glhd,vlhd,sthd,chhd,bthd)

****  x * par / alfa  ****
      XPDF(0) = glhd + glpl
      XPDF(2) = 0.5d0*vlhd + sthd + uppl
      XPDF(1) = 0.5d0*vlhd + sthd + dnpl
      XPDF(3) = sthd + dnpl
      XPDF(4) = chhd + chpl
      XPDF(5) = bthd + btpl

      do 1 I=1,5
 1       XPDF(-I) = XPDF(I)
      do 3 I=-5,5
 3       if (XPDF(I).LT.0.d0) XPDF(I)=0.d0

      F2alfa = 0.d0
      if (OPT.LE.2) goto 2

      glu = alfa*XPDF(0)/x
      up  = alfa*XPDF(2)/x
      dn  = alfa*XPDF(1)/x
      str = alfa*XPDF(3)/x
      chm = alfa*XPDF(4)/x
      bot = alfa*XPDF(5)/x

      mc2 = mc*mc
      ec2 = e2*e2
      ec4 = ec2*ec2

      mb2 = mb*mb
      eb2 = e1*e1
      eb4 = eb2*eb2

      zetac = x*(1.d0+4.d0*mc2/Q2)
      zetab = x*(1.d0+4.d0*mb2/Q2)

      if (zetac.GE.1.d0) chm = 0.d0
      if (zetab.GE.1.d0) bot = 0.d0

**********************************************************************
****   The Bethe-Heitler cross-section for gamma*gamma -> ccbar   ****
**********************************************************************
      beta2 = 1.d0-4.d0*mc2*x/((1.d0-x)*Q2)
      if ((x.LT.1.d0).AND.(beta2.GT.0.d0)) then
         beta = DSQRT(beta2)
         mcQ = 4.d0*mc2/Q2
         F2c = x*3.d0*ec4*alfa/Pi*(
     %       beta*( -1.d0+8.d0*x*(1.d0-x)-x*(1.d0-x)*mcQ )
     %       + ( x*x+(1.-x)*(1.-x)+x*(1.-3.*x)*mcQ
     %       - x*x*mcQ*mcQ/2. )
     %       *DLOG((1.+beta)/(1.-beta)) )
      else
         F2c = 0.d0
      endif
**********************************************************************
****      Substraction of the overlaping Bethe-Heitler term       ****
**********************************************************************
      if (zetac.LT.1.d0) then
         F2cover = zetac*3.d0*ec4*alfa/Pi
     %           *(zetac*zetac+(1.d0-zetac)*(1.d0-zetac))*DLOG(Q2/mc2)
      else
         F2cover = 0.d0
      endif

**********************************************************************
****   The Bethe-Heitler cross-section for gamma*gamma -> bbbar   ****
**********************************************************************
      beta2 = 1.d0-4.d0*mb2*x/((1.d0-x)*Q2)
      if ((x.LT.1.d0).AND.(beta2.GT.0.d0)) then
         beta = DSQRT(beta2)
         mbQ = 4.d0*mb2/Q2
         F2b = x*3.d0*eb4*alfa/Pi*(
     %       beta*( -1.d0+8.d0*x*(1.d0-x)-x*(1.d0-x)*mbQ )
     %       + ( x*x+(1.-x)*(1.-x)+x*(1.-3.*x)*mbQ
     %       - x*x*mbQ*mbQ/2. )
     %       *DLOG((1.+beta)/(1.-beta)) )
      else
         F2b = 0.d0
      endif
**********************************************************************
****      Substraction of the overlaping Bethe-Heitler term       ****
**********************************************************************
      if (zetab.LT.1.d0) then
         F2bover = zetab*3.d0*eb4*alfa/Pi
     %           *(zetab*zetab+(1.d0-zetab)*(1.d0-zetab))*DLOG(Q2/mb2)
      else
         F2bover = 0.d0
      endif

**********************************************************************
****       1:   CHARM - 'resolved' = gamma*G -> ccbar             ****
****       2:   Substraction of the overlapping term              ****
**********************************************************************

      if (zetac.LT.1.d0) then
        step = 5
        eps = 1.d-6
        flav = 1
        call intxr(x,step,zetac,Q2,zetac,1.d0,step,eps,res,t)
        ALS = alfas(Q2,3,Lam3)
        F2cres = ec2 * ALS/(2.d0*Pi) * res(1)
        F2cresover = zetac*ec2 * ALS/Pi * res(2) * DLOG(Q2/mc2)
      else
        F2cres = 0.d0
        F2cresover = 0.d0
      endif

**********************************************************************
****       1:   BEAUTY - 'resolved' = gamma*G -> bbbar            ****
****       2:   Substraction of the overlapping term              ****
**********************************************************************

      if (zetab.LT.1.d0) then
        step = 5
        eps = 1.d-6
        flav = 2
        call intxr(x,step,zetab,Q2,zetab,1.d0,step,eps,res,t)
        ALS = alfas(Q2,3,Lam3)
        F2bres = eb2 * ALS/(2.d0*Pi) * res(1)
        F2bresover = zetab*eb2 * ALS/Pi * res(2) * DLOG(Q2/mb2)
      else
        F2bres = 0.d0
        F2bresover = 0.d0
      endif

**********************************************************************

      F2charm = 2.d0*x*ec2*chm + F2c - F2cover + F2cres - F2cresover
      if (F2charm.LT.0.d0) F2charm = 0.d0

      F2beauty = 2.d0*x*eb2*bot + F2b - F2bover + F2bres - F2bresover
      if (F2beauty.LT.0.d0) F2beauty = 0.d0

      F2 = 2.d0*x*(ec2*up + eb2*dn + eb2*str) + F2charm + F2beauty

      F2 = F2/alfa

 2    continue

      end

*********************************************************************
*********************************************************************

      SUBROUTINE POINTLIKE(OPT,x,Q2,glpl,uppl,dnpl,chpl,btpl)
      IMPLICIT DOUBLE PRECISION (a-z)
      dimension pargl(19),parup(19),pardn(19)
      dimension parch1(21),parch2(21),parbt1(21),parbt2(21)
      integer OPT

****  gluon       ****
      data pargl/ -0.43865d0, 2.7174d0, 0.36752d0, 0.086893d0,
     % 0.010556d0, -0.099005d0, 1.0648d0, 3.6717d0, 2.1944d0,
     % 0.236795d0, -0.19994d0, -0.34992d0, 0.049525d0, 0.34830d0,
     % 0.14342d0, 2.5071d0, 1.9358d0, -0.11849d0, 0.028124d0 /

****  up          ****
      data parup/ -1.0711d0, 3.1320d0, 0.69243d0,
     % -0.058266d0, 0.0097377d0, -0.0068345d0,
     % 0.22297d0, 6.4289d0, 1.7302d0, 0.87940d0,
     % 2.6878d0, 0.20506d0, -0.10617d0, 0.15211d0,
     % 0.013567d0, 2.2802d0, 0.76997d0, -0.110241d0,
     % -0.040252d0 /

****  down = str  ****
      data pardn/ -1.1357d0, 3.1187d0, 0.66290d0,
     % 0.098814d0, -0.092892d0, -0.0066140d0,
     % -0.31385d0, 6.4671d0, 1.6996d0, 11.777d0,
     % -11.124d0, -0.067300d0, 0.049949d0, 0.020427d0,
     % -0.0037558d0, 2.2834d0, 0.84262d0, 0.034760d0,
     % -0.20135d0 /

****  charm  ****
      data parch1/ 2.9808d0, 28.682d0, 2.4863d0,
     % -0.18826d0, 0.18508d0, -0.0014153d0, -0.48961d0,
     % 0.20911d0, 2.7644d0, -7.6307d0, 394.58d0,
     % 0.13565d0, -0.11764d0, -0.011510d0, 0.18810d0,
     % -2.8544d0, 0.93717d0, 5.6807d0, -541.82d0,
     % 14.256d0, 200.82d0 /

      data parch2/ -1.8095d0, 7.9399d0, 0.041563d0,
     % -0.54831d0, 0.19484d0, -0.39046d0, 0.12717d0,
     % 8.7191d0, 4.2616d0, -0.30307d0, 7.2383d0,
     % 0.33412d0, 0.041562d0, 0.37194d0, 0.059280d0,
     % 3.0194d0, 0.73993d0, 0.29430d0, -1.5995d0,
     % 0.d0, 0.d0 /

****  bottom  ****
      data parbt1/ 2.2849d0, 6.0408d0, -0.11577d0,
     % -0.26971d0, 0.27033d0, 0.0022862d0, 0.30807d0,
     % 14.812d0, 1.7148d0, 3.8140d0, 2.2292d0,
     % 0.17942d0, -0.18358d0, -0.0016837d0, -0.10490d0,
     % -1.2977d0, 2.3532d0, -1.0514d0, 20.194d0,
     % 0.0061059d0, 0.053734d0 /

      data parbt2/ -5.0607d0, 16.590d0, 0.87190d0,
     % -0.72790d0, -0.62903d0, -2.4467d0, 0.56575d0,
     % 1.4687d0, 1.1706d0, -0.084651d0, 9.6036d0,
     % 0.36549d0, 0.56817d0, 1.6783d0, -0.19120d0,
     % 9.6071d0, 0.99674d0, -0.083206d0, -3.4864d0,
     % 0.d0,0.d0 /

      glpl = pl(x,Q2,pargl)
      uppl = pl(x,Q2,parup)
      dnpl = pl(x,Q2,pardn)

      if (OPT.EQ.1) goto 1

      if (Q2.LE.10.d0) then
         chpl = cpl(x,Q2,parch1)
      else
         chpl = cpl(x,Q2,parch2)
      endif

      if (Q2.LE.100.d0) then
         btpl = bpl(x,Q2,parbt1)
      else
         btpl = bpl(x,Q2,parbt2)
      endif

 1    continue

      END

**********************
**********************

      SUBROUTINE HADRONLIKE(OPT,x,Q2,glhd,vlhd,sthd,chhd,bthd)
      IMPLICIT DOUBLE PRECISION (a-z)
      dimension pargl(18),parvl(10),parst(14)
      dimension parch1(17),parch2(17),parbt1(16),parbt2(16)
      integer OPT

****  gluon       ****
      data pargl/ 0.59945d0, 1.1285d0, -0.19898d0,
     % 1.9942d0, -1.9848d0, -0.34948d0, 1.0012d0,
     % 1.2287d0, 4.9230d0, 0.21294d0, 0.57414d0,
     % -1.8306d0, 1.4136d0, 0.47058d0, 0.99767d0,
     % 2.4447d0, 0.18526d0, 2.7450d0 /

****  valence     ****
      data parvl/ 1.0898d0, 0.78391d0, 0.42654d0,
     % -1.6576d0, 0.96155d0, 0.38087d0, -0.068720d0,
     % -1.2128d0, 1.7075d0, 1.8441d0 /

****  strange     ****
      data parst/ 0.71660d0, 0.72289d0, 0.60478d0,
     % 4.2106d0, 4.1494d0, 4.5179d0, 5.2812d0,
     % 1.0497d0, -0.21562d0, 0.036160d0, -0.85835d0,
     % 0.34866d0, 1.9219d0, -0.15200d0 /

****  charm  ****
      data parch1/ 5.6729d0, 1.6248d0, -2586.4d0, 2695.0d0,
     % 1.5146d0, -3.9185d0, 3.6126d0, 1.4575d0,
     % -0.70433d0, 1910.1d0, -1688.2d0, 3.1028d0,
     % 11.738d0, -1.0291d0, 0.d0, 0.d0, 0.d0 /

      data parch2/ -1.6470d0, -0.78809d0, -2.0561d0,
     % 2.1266d0, 3.0301d0, 4.1282d0, 0.89599d0,
     % 0.72738d0, 0.90278d0, 0.75576d0, 0.66383d0,
     % -1.7499d0, 1.6929d0, 1.2761d0, -0.15061d0,
     % -0.26292d0, 1.6466d0 /

****  bottom  ****
      data parbt1/ -10.210d0, 0.82278d0, -99.613d0,
     % 492.61d0, 3.3917d0, 5.6829d0, -2.0137d0,
     % -2.2296d0, 0.081818d0, 171.25d0, -420.45d0,
     % 0.084256d0, -0.23571d0, 4.6955d0, 0.d0, 0.d0 /

      data parbt2/ 2.4198d0, -0.98933d0, -2.1109d0,
     % 9.0196d0, 3.6455d0, 4.6196d0, 0.66454d0,
     % 0.40703d0, 0.42366d0, 1.2711d0, -3.6082d0,
     % -4.1353d0, 2.4212d0, 1.1109d0, 0.15817d0, 2.3615d0 /

      glhd = glu(x,Q2,pargl)
      vlhd = val(x,Q2,parvl)
      sthd = str(x,Q2,parst)

      if (OPT.EQ.1) goto 1

      if (Q2.LE.10.d0) then
         chhd = chm(x,Q2,parch1)
      else
         chhd = chm(x,Q2,parch2)
      endif

      if (Q2.LE.100.d0) then
         bthd = bot(x,Q2,parbt1)
      else
         bthd = bot(x,Q2,parbt2)
      endif

 1    continue

      END

*********************************************************************

      SUBROUTINE GLUON(x,Q2,glun)
      IMPLICIT DOUBLE PRECISION (a-z)
      Parameter (alfa=7.29735308D-3)
      dimension plglu(19),hadglu(18)

****  point-like       ****
      data plglu/ -0.4387d0, 2.717d0, 0.3675d0,
     % 0.08689d0, 0.01056d0, -0.09900d0, 1.065d0,
     % 3.672d0, 2.194d0, 0.2368d0, -0.1999d0,
     % -0.3499d0, 0.04953d0, 0.3483d0, 0.1434d0,
     % 2.507d0, 1.936d0, -0.1185d0, 0.02812d0 /

****  hadron-like      ****
      data hadglu/ 0.599449d0, 1.12849d0, -0.198980d0, 1.99417d0,
     % -1.98480d0, -0.349479d0, 1.00120d0, 1.22871d0, 4.92304d0,
     % 0.212938d0, 0.574140d0, -1.83060d0, 1.41359d0, 0.470584d0,
     % 0.997665d0, 2.44470d0, 0.185258d0, 2.74500d0 /

      glpl = pl(x,Q2,plglu)
      glhd = glu(x,Q2,hadglu)

      glun = alfa*(glpl + glhd)/x

      END

********************************************************************
****                      POINT-LIKE                            ****
********************************************************************

      DOUBLE PRECISION FUNCTION pl(x,Q2,PAR)
      IMPLICIT DOUBLE PRECISION (a-z)
      dimension PAR(19)
      Parameter (Pi=3.1415926535897932d0)

      Lam2 = 0.221d0*0.221d0
      s = DLOG(DLOG(Q2/Lam2)/DLOG(0.25d0/Lam2))
      dlg = DLOG(1.d0/x)

      alfa  = PAR(1)
      alfap = PAR(2)
      beta  = PAR(3)
      A =  PAR(4)  + PAR(12)*s
      B =  PAR(5)  + PAR(13)*s
      C =  PAR(6)  + PAR(14)*s
      D =  PAR(7)  + PAR(15)*s
      E =  PAR(8)  + PAR(16)*s
      EP = PAR(9)  + PAR(17)*s
      AS = PAR(10) + PAR(18)*s
      BS = PAR(11) + PAR(19)*s

      pl = s**alfa*x**AS*(A+B*DSQRT(x)+C*x**BS)
      pl = pl + s**alfap*DEXP(-E + DSQRT(EP*s**beta*dlg))
      pl = 9.d0/(4.d0*Pi)*DLOG(Q2/Lam2)*pl*(1.d0-x)**D

      END

**********
**********

      DOUBLE PRECISION FUNCTION cpl(x,Q2,PAR)
      IMPLICIT DOUBLE PRECISION (a-z)
      dimension PAR(21)
      Parameter (Pi=3.1415926535897932d0)

      Lam2 = 0.221d0*0.221d0
      s = DLOG(DLOG(Q2/Lam2)/DLOG(0.25d0/Lam2))
      s2 = s*s

      cpl = 0.d0
      y = x + 1.d0 - Q2/(Q2+6.76d0)
      if (y.GE.1.d0) goto 10

      dlg = DLOG(1.d0/x)

      alfa  = PAR(1)
      alfap = PAR(2)
      beta  = PAR(3)
      A =  PAR(4)  + PAR(12)*s
      B =  PAR(5)  + PAR(13)*s
      C =  PAR(6)  + PAR(14)*s
      D =  PAR(7)  + PAR(15)*s
      E =  PAR(8)  + PAR(16)*s + PAR(20)*s2
      EP = PAR(9)  + PAR(17)*s
      AS = PAR(10) + PAR(18)*s
      BS = PAR(11) + PAR(19)*s + PAR(21)*s2

      cpl = s**alfa*y**AS*(A+B*DSQRT(y)+C*y**BS)
      cpl = cpl + s**alfap*DEXP(-E + DSQRT(EP*s**beta*dlg))
      cpl = 9.d0/(4.d0*Pi)*DLOG(Q2/Lam2)*cpl*(1.d0-y)**D

 10   continue

      END

**********
**********

      DOUBLE PRECISION FUNCTION bpl(x,Q2,PAR)
      IMPLICIT DOUBLE PRECISION (a-z)
      dimension PAR(21)
      Parameter (Pi=3.1415926535897932d0)

      Lam2 = 0.221d0*0.221d0
      s = DLOG(DLOG(Q2/Lam2)/DLOG(0.25d0/Lam2))
      s2 = s*s
      ds = DSQRT(s)

      bpl = 0.d0
      y = x + 1.d0 - Q2/(Q2+73.96d0)
      if (y.GE.1.d0) goto 10

      dlg = DLOG(1.d0/x)

      alfa  = PAR(1)
      alfap = PAR(2)
      beta  = PAR(3)
      A =  PAR(4)  + PAR(12)*s
      B =  PAR(5)  + PAR(13)*s + PAR(20)*s2
      C =  PAR(6)  + PAR(14)*s
      D =  PAR(7)  + PAR(15)*s
      E =  PAR(8)  + PAR(16)*s
      EP = PAR(9)  + PAR(17)*s + PAR(21)*ds
      AS = PAR(10) + PAR(18)*s
      BS = PAR(11) + PAR(19)*s

      bpl = s**alfa*y**AS*(A+B*DSQRT(y)+C*y**BS)
      bpl = bpl + s**alfap*DEXP(-E + DSQRT(EP*s**beta*dlg))
      bpl = 9.d0/(4.d0*Pi)*DLOG(Q2/Lam2)*bpl*(1.d0-y)**D

 10   continue

      END

********************************************************************
****                     HADRON-LIKE                            ****
********************************************************************

      DOUBLE PRECISION FUNCTION val(x,Q2,PAR)
      IMPLICIT DOUBLE PRECISION (a-z)
      dimension PAR(10)
      Parameter (Pi=3.1415926535897932d0)

      Lam2 = 0.221d0*0.221d0
      s = DLOG(DLOG(Q2/Lam2)/DLOG(0.25d0/Lam2))

      AC = PAR(1) + PAR(6)*s
      AS = PAR(2) + PAR(7)*s
      BC = PAR(3) + PAR(8)*s
      C  = PAR(4) + PAR(9)*s
      D  = PAR(5) + PAR(10)*s

      val = AC*x**AS*(1.d0+BC*DSQRT(x)+C*x)
      val = val*(1.d0-x)**D

      END

**********
**********

      DOUBLE PRECISION FUNCTION glu(x,Q2,PAR)
      IMPLICIT DOUBLE PRECISION (a-z)
      dimension PAR(18)
      Parameter (Pi=3.1415926535897932d0)

      Lam2 = 0.221d0*0.221d0
      s = DLOG(DLOG(Q2/Lam2)/DLOG(0.25d0/Lam2))
      dlg = DLOG(1.d0/x)

      alfa = PAR(1)
      beta = PAR(2)
      AC = PAR(3)  + PAR(11)*s
      BC = PAR(4)  + PAR(12)*s
      C  = PAR(5)  + PAR(13)*s
      AS = PAR(6)  + PAR(14)*s
      BS = PAR(7)  + PAR(15)*s
      E  = PAR(8)  + PAR(16)*s
      EP = PAR(9)  + PAR(17)*s
      D  = PAR(10) + PAR(18)*s

      glu = x**AS*(AC+BC*DSQRT(x)+C*x)
      glu = glu + s**alfa*DEXP(-E+DSQRT(EP*s**beta*dlg))
      glu = glu*(1.d0-x)**D

      END

**********
**********

      DOUBLE PRECISION FUNCTION str(x,Q2,PAR)
      IMPLICIT DOUBLE PRECISION (a-z)
      dimension PAR(14)
      Parameter (Pi=3.1415926535897932d0)

      Lam2 = 0.221d0*0.221d0
      s = DLOG(DLOG(Q2/Lam2)/DLOG(0.25d0/Lam2))
      dlg = DLOG(1.d0/x)

      alfa = PAR(1)
      AS = PAR(2)  + PAR(9)*s
      AC = PAR(3)  + PAR(10)*s
      BC = PAR(4)  + PAR(11)*s
      D  = PAR(5)  + PAR(12)*s
      E  = PAR(6)  + PAR(13)*s
      EP = PAR(7)  + PAR(14)*s
      beta = PAR(8)

      str = s**alfa/(dlg**AS)*(1.d0+AC*DSQRT(x)+BC*x)
      str = str*(1.d0-x)**D*DEXP(-E+DSQRT(EP*s**beta*dlg))

      END

**********
**********

      DOUBLE PRECISION FUNCTION chm(x,Q2,PAR)
      IMPLICIT DOUBLE PRECISION (a-z)
      dimension PAR(17)
      Parameter (Pi=3.1415926535897932d0)

      Lam2 = 0.221d0*0.221d0
      s = DLOG(DLOG(Q2/Lam2)/DLOG(0.25d0/Lam2))
      s2 = s*s

      chm = 0.d0
      y = x + 1.d0 - Q2/(Q2+6.76d0)
      if (y.GE.1.d0) goto 10

      dlg = DLOG(1.d0/x)

      alfa = PAR(1)
      AS = PAR(2)  + PAR(9)*s
      AC = PAR(3)  + PAR(10)*s
      BC = PAR(4)  + PAR(11)*s
      D  = PAR(5)  + PAR(12)*s + PAR(17)*s2
      E  = PAR(6)  + PAR(13)*s + PAR(16)*s2
      EP = PAR(7)  + PAR(14)*s + PAR(15)*s2
      beta = PAR(8)

      chm = s**alfa/(dlg**AS)*(1.d0+AC*DSQRT(y)+BC*y)
      chm = chm*(1.d0-y)**D*DEXP(-E+EP*DSQRT(s**beta*dlg))

 10   continue

      END

**********
**********

      DOUBLE PRECISION FUNCTION bot(x,Q2,PAR)
      IMPLICIT DOUBLE PRECISION (a-z)
      dimension PAR(16)
      Parameter (Pi=3.1415926535897932d0)

      Lam2 = 0.221d0*0.221d0
      s = DLOG(DLOG(Q2/Lam2)/DLOG(0.25d0/Lam2))
      s2 = s*s

      bot = 0.d0
      y = x + 1.d0 - Q2/(Q2+73.96d0)
      if (y.GE.1.d0) goto 10

      dlg = DLOG(1.d0/x)

      alfa = PAR(1)
      AS = PAR(2)  + PAR(9)*s  + PAR(15)*s2
      AC = PAR(3)  + PAR(10)*s
      BC = PAR(4)  + PAR(11)*s
      D  = PAR(5)  + PAR(12)*s + PAR(16)*s2
      E  = PAR(6)  + PAR(13)*s
      EP = PAR(7)  + PAR(14)*s
      beta = PAR(8)

      bot = s**alfa/(dlg**AS)*(1.d0+AC*DSQRT(y)+BC*y)
      bot = bot*(1.d0-y)**D*DEXP(-E+EP*DSQRT(s**beta*dlg))

 10   continue

      END

*************************************************************************
****                    Running alpha strong                         ****
*************************************************************************

      DOUBLE PRECISION FUNCTION alfas(Q2,Nf,LAM)
      IMPLICIT DOUBLE PRECISION (a-z)
      PARAMETER(PI=3.1415926535898d0)
      INTEGER Nf

      alfas = 12.0d0*PI/((33.0d0-2.0d0*Nf)*DLOG(Q2/(LAM*LAM)) )

      END

****************************************************************************

      subroutine intxr(x0,step,x,Q2,a,b,n,eps,result,t)
      implicit double precision (a-z)
      integer *4 iadr(60)
      integer n,I,ind,step
      dimension aa(60),bb(60),result(2),reslt2(2),eps1(2),eps2(2),
     %          ra1(2),ra2(2),rb1(2),rb2(2)
      logical t
      t=.true.

      DO 110 I=1,2
         result(I)=0.d0
110      reslt2(I)=0.d0

      ind=1
      iadr(1)=-1
      aa(1)=a
      bb(1)=b
1     c=(aa(ind)+bb(ind))/2.0d0

      call gauscxr(x0,step,x,Q2,aa(ind),c,ra1,ra2)
      call gauscxr(x0,step,x,Q2,c,bb(ind),rb1,rb2)

      DO 200 I=1,2
        eps1(I)=dabs(ra1(I)-ra2(I))/(dabs(ra1(I)+result(I))+1.0d-300)
200     eps2(I)=dabs(rb1(I)-rb2(I))/(dabs(rb1(I)+result(I))+1.0d-300)

      rozn = eps1(1) - eps2(1)
      DO 300 I=2,2
         roz = eps1(I)-eps2(I)
         if (roz.GT.rozn) rozn = roz
300   CONTINUE
      if(rozn) 10,10,20

10    rozn = eps1(1) - eps
      DO 400 I=2,2
         roz = eps1(I) - eps
         if (roz.GT.rozn) rozn = roz
400   CONTINUE
      if(rozn) 12,12,11

11    if(ind-n) 13,15,15
15    t=.false.

12    DO 500 I=1,2
         result(I)=result(I)+ra1(I)
500      reslt2(I)=reslt2(I)+ra2(I)

      iadr(ind)=iadr(ind)+100
      if(iadr(ind)-150) 20,20,30
13    ind=ind+1
      iadr(ind)=0
      aa(ind)=aa(ind-1)
      bb(ind)=(aa(ind-1)+bb(ind-1))/2.
      go to 1
14    iadr(ind)=iadr(ind)+100
      if(iadr(ind)-150) 23,23,30

20    rozn = eps2(1) - eps
      DO 600 I=2,2
         roz = eps2(I) - eps
         if (roz.GT.rozn) rozn = roz
600   CONTINUE
      if(rozn) 22,22,21

21    if(ind-n) 23,25,25
25    t =.false.

22    DO 700 I=1,2
         result(I)=result(I)+rb1(I)
700      reslt2(I)=reslt2(I)+rb2(I)

      iadr(ind)=iadr(ind)+100
      if(iadr(ind)-150) 10,10,30
23    ind=ind+1
      iadr(ind)=1
      aa(ind)=(aa(ind-1)+bb(ind-1))/2.
      bb(ind)=bb(ind-1)
      go to 1
24    iadr(ind)=iadr(ind)+100
      if(iadr(ind)-150) 13,13,30
30    ind=ind-1
      if(iadr(ind+1)-200) 100,14,24

100   eps = dabs(result(1)-reslt2(1))/(dabs(result(1))+1.d-300)
      DO 800 I=2,2
         ep = dabs(result(I)-reslt2(I))/(dabs(result(I)+1.d-300))
         if (ep.GT.eps) eps = ep
800   CONTINUE

      return
      end

****************************************************************************

      subroutine gauscxr(x0,step,x,Q2,a,b,gauskr,gaus)
      implicit double precision(a-z)
      parameter (nvar = 12)
      integer l,I,step
      dimension g(3,8),res1(2),res2(2),res3(2),
     %          gauskr(2),gaus(2),c(2)
      Parameter (Pi=3.1415926535897932d0, alfa=7.29735308D-3)
      data g/
     $9.933798 7588 1716d-1,0.,1.782238 3320 7104d-2,
     $9.602898 5649 7536d-1,1.012285 36290376d-1,4.943939 5002 1394d-2,
     $8.941209 0684 7456d-1,0.  ,8.248229 8931 3584d-2,
     $7.966664 7741 3626d-1,2.223810 3445 3374d-1,1.116463 7082 6840d-1,
     $6.723540 7094 5158d-1,0.  ,1.362631 0925 5172d-1,
     $5.255324 0991 6329d-1,3.137066 4587 7887d-1,1.566526 0616 8188d-1,
     $3.607010 9792 8132d-1,0.  ,1.720706 0855 5211d-1,
     $1.834346 4249 5650d-1,3.626837 8337 8362d-1,1.814000 2506 8035d-1/
      data g39/1.844464 0574 4692d-1/

      do 10 I=1,2
        gaus(I)=0.d0
 10     gauskr(I)=0.d0

      d=(b-a)/2.
      e=(b+a)/2.

      z3 = e

      call GLUON(z3,Q2,glu3)
      res3(1) = glu3*fun(x0/z3,Q2)
      res3(2) = glu3/z3*Pqg(x/z3)

      do 100 l=1,8
      y=d*g(1,l)

      z1 = e+y
      z2 = e-y

      call GLUON(z1,Q2,glu1)
      res1(1) = glu1*fun(x0/z1,Q2)
      res1(2) = glu1/z1*Pqg(x/z1)

      call GLUON(z2,Q2,glu2)
      res2(1) = glu2*fun(x0/z2,Q2)
      res2(2) = glu2/z2*Pqg(x/z2)

      do 20 I=1,2
         c(I) = res1(I) + res2(I)
         gaus(I) = gaus(I)+c(I)*g(2,l)
 20      gauskr(I) = gauskr(I)+c(I)*g(3,l)
100   continue
      do 30 I=1,2
         gaus(I) = d*gaus(I)
 30      gauskr(I) = d*(gauskr(I)+g39*res3(I))
      return
      end

****************************************************************************

      DOUBLE PRECISION FUNCTION fun(X,Q2)
      IMPLICIT DOUBLE PRECISION (a-z)
      integer flav
      common /flav/ flav
      common /mass/ mc,mb

      if (flav.EQ.1) then

        mc2 = mc*mc
        beta2 = 1.d0-4.d0*mc2*X / ((1.d0-X)*Q2)

        if (beta2.GT.0.d0) then
           mcQ = 4.d0*mc2/Q2
           beta = DSQRT(beta2)
           fun = X*(
     %         beta*( -1.d0+8.d0*X*(1.d0-X)-X*(1.d0-X)*mcQ )
     %         + ( X*X+(1.-X)*(1.-X)+X*(1.-3.*X)*mcQ
     %         - X*X*mcQ*mcQ/2. )
     %         *DLOG((1.+beta)/(1.-beta)) )
        else
           fun = 0.d0
        endif

      elseif (flav.EQ.2) then

        mb2 = mb*mb
        beta2 = 1.d0-4.d0*mb2*X / ((1.d0-X)*Q2)

        if (beta2.GT.0.d0) then
           mbQ = 4.d0*mb2/Q2
           beta = DSQRT(beta2)
           fun = X*(
     %         beta*( -1.d0+8.d0*X*(1.d0-X)-X*(1.d0-X)*mbQ )
     %         + ( X*X+(1.-X)*(1.-X)+X*(1.-3.*X)*mbQ
     %         - X*X*mbQ*mbQ/2. )
     %         *DLOG((1.+beta)/(1.-beta)) )
        else
           fun = 0.d0
        endif

      endif


      END

****************************************************************************

      DOUBLE PRECISION FUNCTION Pqg(X)
      IMPLICIT DOUBLE PRECISION (a-z)

      Pqg = 1.d0/2.d0*( X*X + (1.d0-X)*(1.d0-X) )

      END
