CDECK  ID>, SSRSGB. 
C--------------------------------------------------------------------
      FUNCTION SSRSGB(QS)
C--------------------------------------------------------------------
C
C     Calculate b quark self energy
C     according to Pierce et al. formulae adapted to Isajet
C
C
      IMPLICIT NONE
      COMMON/SSLUN/LOUT
      INTEGER LOUT
      SAVE /SSLUN/
C          SUSY parameters
C          AMGLSS               = gluino mass
C          AMULSS               = up-left squark mass
C          AMELSS               = left-selectron mass
C          AMERSS               = right-slepton mass
C          AMNiSS               = sneutrino mass for generation i
C          TWOM1                = Higgsino mass = - mu
C          RV2V1                = ratio v2/v1 of vev's
C          AMTLSS,AMTRSS        = left,right stop masses
C          AMT1SS,AMT2SS        = light,heavy stop masses
C          AMBLSS,AMBRSS        = left,right sbottom masses
C          AMB1SS,AMB2SS        = light,heavy sbottom masses
C          AMLLSS,AMLRSS        = left,right stau masses
C          AML1SS,AML2SS        = light,heavy stau masses
C          AMZiSS               = signed mass of Zi
C          ZMIXSS               = Zi mixing matrix
C          AMWiSS               = signed Wi mass
C          GAMMAL,GAMMAR        = Wi left, right mixing angles
C          AMHL,AMHH,AMHA       = neutral Higgs h0, H0, A0 masses
C          AMHC                 = charged Higgs H+ mass
C          ALFAH                = Higgs mixing angle
C          AAT                  = stop trilinear term
C          THETAT               = stop mixing angle
C          AAB                  = sbottom trilinear term
C          THETAB               = sbottom mixing angle
C          AAL                  = stau trilinear term
C          THETAL               = stau mixing angle
C          AMGVSS               = gravitino mass
C          MTQ                  = top mass at MSUSY
C          MBQ                  = bottom mass at MSUSY
C          MLQ                  = tau mass at MSUSY
C          FBMA                 = b-Yukawa at mA scale
C          VUQ                  = Hu vev at MSUSY
C          VDQ                  = Hd vev at MSUSY
      COMMON/SSPAR/AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS(4,4)
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,
     $VUQ,VDQ
      REAL AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,VUQ,VDQ
      REAL AMZISS(4)
      EQUIVALENCE (AMZISS(1),AMZ1SS)
      SAVE /SSPAR/
C          Standard model parameters
C          AMUP,...,AMTP        = quark masses
C          AME,AMMU,AMTAU       = lepton masses
C          AMW,AMZ              = W,Z masses
C          GAMW,GAMZ            = W,Z widths
C          ALFAEM,SN2THW,ALFA3  = SM couplings
C          ALQCD4               = 4 flavor lambda
      COMMON/SSSM/AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      REAL AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      SAVE /SSSM/
C          Frozen couplings from RG equations:
C     GSS( 1) = g_1        GSS( 2) = g_2        GSS( 3) = g_3
C     GSS( 4) = y_tau      GSS( 5) = y_b        GSS( 6) = y_t
C     GSS( 7) = M_1        GSS( 8) = M_2        GSS( 9) = M_3
C     GSS(10) = A_tau      GSS(11) = A_b        GSS(12) = A_t
C     GSS(13) = M_h1^2     GSS(14) = M_h2^2     GSS(15) = M_er^2
C     GSS(16) = M_el^2     GSS(17) = M_dnr^2    GSS(18) = M_upr^2
C     GSS(19) = M_upl^2    GSS(20) = M_taur^2   GSS(21) = M_taul^2
C     GSS(22) = M_btr^2    GSS(23) = M_tpr^2    GSS(24) = M_tpl^2
C     GSS(25) = mu         GSS(26) = B          GSS(27) = Y_N
C     GSS(28) = M_nr       GSS(29) = A_n        GSS(30) = log(vdq)
C     GSS(31) = log(vuq)
C          Masses:
C     MSS( 1) = glss     MSS( 2) = upl      MSS( 3) = upr
C     MSS( 4) = dnl      MSS( 5) = dnr      MSS( 6) = stl
C     MSS( 7) = str      MSS( 8) = chl      MSS( 9) = chr
C     MSS(10) = b1       MSS(11) = b2       MSS(12) = t1
C     MSS(13) = t2       MSS(14) = nuel     MSS(15) = numl
C     MSS(16) = nutl     MSS(17) = el-      MSS(18) = er-
C     MSS(19) = mul-     MSS(20) = mur-     MSS(21) = tau1
C     MSS(22) = tau2     MSS(23) = z1ss     MSS(24) = z2ss
C     MSS(25) = z3ss     MSS(26) = z4ss     MSS(27) = w1ss
C     MSS(28) = w2ss     MSS(29) = hl0      MSS(30) = hh0
C     MSS(31) = ha0      MSS(32) = h+
C          Unification:
C     MGUTSS  = M_GUT    GGUTSS  = g_GUT    AGUTSS  = alpha_GUT
      COMMON /SUGMG/ MSS(32),GSS(31),MGUTSS,GGUTSS,AGUTSS,FTGUT,
     $FBGUT,FTAGUT,FNGUT
      REAL MSS,GSS,MGUTSS,GGUTSS,AGUTSS,FTGUT,FBGUT,FTAGUT,FNGUT
      SAVE /SUGMG/
      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,ASM3
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3
      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG
      SAVE /SUGPAS/
      COMPLEX ZAZB1(4),ZAZB2(4),ZBZB1(4),ZBZB2(4),ZI
      COMPLEX ZAWD(2)
      COMPLEX ZAWT1(2),ZAWT2(2),ZBWT1(2),ZBWT2(2)
      COMPLEX*16 SSB0,SSB1,ZZZ,SIGWI,SIGZI
      REAL SR2,PI,E,G,TANB,BETA,SINB,COSB,COSA,SINA,COST,SINT,
     $MG,MT1,MT2,MB1,MB2,GLB,GRB,XM,YM,THX,THY,
     $FB1(4),FB2(4),GB1(4),GB2(4),MZI(4),FT,FB,
     $RSIGB,QS,SSRSGB,ZBW(2),BWP(2),FWT1(2),GWT1(2),FWT2(2),
     $GWT2(2),MWI(2),MW,MZ,MB,MT,SUALFS,COS2W,FAC,
     $COSBE,SINBE,XM3
      REAL*8 LP(11),LPTOT
      REAL*8 REAL8,RSIGWI,RSIGZI

      INTEGER THZ(4),THW(2),I
C
      DATA ZI/(0.,1.)/
C     Recompute weak scale Yukawa couplings including SUSY loops
C     Follow formulae of Pierce et al. NPB491, 3 (1997)
C
      REAL8(ZZZ)=DREAL(ZZZ)
      SR2=SQRT(2.)
      PI=4*ATAN(1.)
      E=SQRT(4*PI*ALFAEM)
      COS2W=1.-SN2THW
      G=G2
      MW=AMW
      MZ=AMZ
      MB=MBQ
C      ASMZ=SUALFS(AMZ**2,.36,AMTP,3)
      FAC=16*PI**2
      TANB=VUQ/VDQ
      BETA=ATAN(TANB)
      SINBE=SIN(BETA)
      COSBE=COS(BETA)
      COSA=COS(ALFAH)
      SINA=SIN(ALFAH)
      COST=COS(THETAT)
      SINT=SIN(THETAT)
      COSB=COS(THETAB)
      SINB=SIN(THETAB)
      MG=ABS(GSS(9))
      XM3=SIGN(1.,MSS(1))*MG
      MT=AMTP
      MT1=MSS(12)
      MT2=MSS(13)
      MB1=MSS(10)
      MB2=MSS(11)
      GLB=-.5+XW/3.
      GRB=-XW/3.
      FT=MTQ/VUQ
      FB=MBQ/VDQ
      XM=1./TAN(GAMMAL)
      YM=1./TAN(GAMMAR)
      THX=SIGN(1.,XM)
      THY=SIGN(1.,YM)
      THW(1)=0
      IF(AMW1SS.LT.0.) THW(1)=1
      THW(2)=0
      IF(AMW2SS.LT.0.) THW(2)=1

C     Neutralino-fermion-sfermion couplings
C     a(Pierce)=-beta^* ; b(Pierce)=-alpha^*
C
      SIGZI=(0.,0.)
      DO I=1,4
        THZ(I)=0
        IF (AMZISS(I).LT.0.) THZ(I)=1
        MZI(I)=ABS(AMZISS(I))
        ZAZB1(I)=ZI*(-ZI)**(THZ(I)-1)*(-G*ZMIXSS(3,I)+GP*
     $  ZMIXSS(4,I)/3.)/SR2*COSB+(-ZI)**(THZ(I))*FB*ZMIXSS(2,I)*SINB
        ZAZB2(I)=ZI*(-ZI)**(THZ(I)-1)*(-G*ZMIXSS(3,I)+GP*
     $  ZMIXSS(4,I)/3.)/SR2*SINB-(-ZI)**(THZ(I))*FB*ZMIXSS(2,I)*COSB
        ZBZB1(I)=-ZI*(ZI)**(THZ(I)-1)*(-2)*GP/3./SR2*ZMIXSS(4,I)*SINB
     $  -(ZI)**(THZ(I))*FB*ZMIXSS(2,I)*COSB
        ZBZB2(I)=ZI*(ZI)**(THZ(I)-1)*(-2)*GP/3./SR2*ZMIXSS(4,I)*COSB
     $  -(ZI)**(THZ(I))*FB*ZMIXSS(2,I)*SINB
        FB1(I)=ZAZB1(I)*CONJG(ZAZB1(I))+ZBZB1(I)*CONJG(ZBZB1(I))
        GB1(I)=2*REAL(ZAZB1(I)*CONJG(ZBZB1(I)))
        FB2(I)=ZAZB2(I)*CONJG(ZAZB2(I))+ZBZB2(I)*CONJG(ZBZB2(I))
        GB2(I)=2*REAL(ZAZB2(I)*CONJG(ZBZB2(I)))
        SIGZI=SIGZI+FB1(I)*SSB1(QS,MZI(I),MB1)+GB1(I)*MZI(I)/MBQ*
     $  SSB0(QS,MZI(I),MB1)+FB2(I)*SSB1(QS,MZI(I),MB2)
     $  +GB2(I)*MZI(I)/MBQ*SSB0(QS,MZI(I),MB2)
      END DO
      RSIGZI=REAL8(SIGZI)
C
C     Chargino-fermion-sfermion couplings; I labels chargino
C
      ZAWD(1)=ZI*(-1)**THW(1)*G*SIN(GAMMAR)
      ZAWD(2)=ZI*(-1)**THW(2)*THY*G*COS(GAMMAR)
      ZBW(1)=-(-1)**THW(1)*FT*COS(GAMMAR)
      ZBW(2)=(-1)**THW(2)*THY*FT*SIN(GAMMAR)
      BWP(1)=-FB*COS(GAMMAL)
      BWP(2)=FB*THX*SIN(GAMMAL)
      SIGWI=(0.,0.)
      MWI(1)=ABS(AMW1SS)
      MWI(2)=ABS(AMW2SS)
      DO I=1,2
        ZAWT1(I)=-ZI*ZAWD(I)*COST+ZBW(I)*SINT
        ZAWT2(I)=-ZI*ZAWD(I)*SINT-ZBW(I)*COST
        ZBWT1(I)=-BWP(I)*COST
        ZBWT2(I)=-BWP(I)*SINT
        FWT1(I)=ZAWT1(I)*CONJG(ZAWT1(I))+ZBWT1(I)*CONJG(ZBWT1(I))
        GWT1(I)=2*REAL(CONJG(ZBWT1(I))*ZAWT1(I))
        FWT2(I)=ZAWT2(I)*CONJG(ZAWT2(I))+ZBWT2(I)*CONJG(ZBWT2(I))
        GWT2(I)=2*REAL(CONJG(ZBWT2(I))*ZAWT2(I))
        SIGWI=SIGWI+FWT1(I)*SSB1(QS,MWI(I),MT1)+
     $MWI(I)/MBQ*GWT1(I)*SSB0(QS,MWI(I),MT1)+
     $FWT2(I)*SSB1(QS,MWI(I),MT2)+
     $MWI(I)/MBQ*GWT2(I)*SSB0(QS,MWI(I),MT2)
      END DO
      RSIGWI=REAL8(SIGWI)
      LP(1)=ASM3/3./PI*(REAL8(SSB1(QS,MG,MB1))+
     $  REAL8(SSB1(QS,MG,MB2))-SIN(2*THETAB)*XM3/MBQ*
     $  (REAL8(SSB0(QS,MG,MB1))-REAL8(SSB0(QS,MG,MB2))))
      LP(2)=.5*FB**2*COSA**2*(REAL8(SSB1(QS,MB,AMHH))+
     $REAL8(SSB0(QS,MB,AMHH)))/FAC
      LP(3)=.5*FB**2*SINA**2*(REAL8(SSB1(QS,MB,AMHL))+
     $REAL8(SSB0(QS,MB,AMHL)))/FAC
      LP(4)=.5*FB**2*SINBE**2*(REAL8(SSB1(QS,MB,AMHA))-
     $REAL8(SSB0(QS,MB,AMHA)))/FAC
      LP(5)=.5*FB**2*COSBE**2*(REAL8(SSB1(QS,MB,MZ))-
     $REAL8(SSB0(QS,MB,MZ)))/FAC
      LP(6)=.5*((FT**2*COSBE**2+FB**2*SINBE**2)*
     $REAL8(SSB1(QS,MT,AMHC))+(G2**2+FT**2*SINBE**2+FB**2*COSBE**2)
     $*REAL8(SSB1(QS,MT,MW)))/FAC
      LP(7)=FT**2*SINBE**2*(REAL8(SSB0(QS,MT,AMHC))-
     $REAL8(SSB0(QS,MT,MW)))/FAC
      LP(8)=0.
      LP(9)=+G2**2/COS2W*((GLB**2+GRB**2)*REAL8(SSB1(QS,MB,MZ))
     $+4*GLB*GRB*REAL8(SSB0(QS,MB,MZ)))/FAC
      LP(10)=.5*RSIGZI/FAC
      LP(11)=.5*RSIGWI/FAC
      LPTOT=0.D0
      DO I=1,11
        LPTOT=LPTOT+LP(I)
C        WRITE(6,*) 'LP(',I,')=',LP(I)
      END DO
      SSRSGB=LPTOT
100   RETURN
      END
