CDECK  ID>, SSRSGL. 
C--------------------------------------------------------------------
      FUNCTION SSRSGL(QS)
C--------------------------------------------------------------------
C
C     Calculate tau lepton self energy
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
      COMPLEX ZAZL1(4),ZAZL2(4),ZBZL1(4),ZBZL2(4),ZI
      COMPLEX ZAWL(2)
      COMPLEX ZAWSN(2),ZBWSN(2)
      COMPLEX*16 SSB0,SSB1,ZZZ,SIGWI,SIGZI
      REAL SR2,PI,E,G,TANB,BETA,SINB,COSB,COSA,SINA,COSL,SINL,
     $MG,ML1,ML2,MSN,GLL,GRL,XM,YM,THX,THY,
     $FL1(4),FL2(4),GL1(4),GL2(4),MZI(4),FL,ML,
     $RSIGL,QS,SSRSGL,BWPP(2),FWSN(2),GWSN(2),
     $MWI(2),MW,MZ,COS2W,FAC,
     $COSBE,SINBE
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
      ML=MLQ
      FAC=16*PI**2
      TANB=VUQ/VDQ
      BETA=ATAN(TANB)
      SINBE=SIN(BETA)
      COSBE=COS(BETA)
      COSA=COS(ALFAH)
      SINA=SIN(ALFAH)
      COSL=COS(THETAL)
      SINL=SIN(THETAL)
      MSN=AMN3SS
      ML1=AML1SS
      ML2=AML2SS
      GLL=-.5+XW
      GRL=-XW
      FL=MLQ/VDQ
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
        ZAZL1(I)=ZI*(-1)*(-ZI)**(THZ(I)-1)*(G*ZMIXSS(3,I)+GP*
     $  ZMIXSS(4,I))/SR2*COSL+(-ZI)**(THZ(I))*FL*ZMIXSS(2,I)*SINL
        ZAZL2(I)=ZI*(-1)*(-ZI)**(THZ(I)-1)*(G*ZMIXSS(3,I)+GP*
     $  ZMIXSS(4,I))/SR2*SINL-(-ZI)**(THZ(I))*FL*ZMIXSS(2,I)*COSL
        ZBZL1(I)=-ZI*(-1)*(ZI)**(THZ(I)-1)*SR2*GP*ZMIXSS(4,I)*SINL
     $  -(ZI)**(THZ(I))*FL*ZMIXSS(2,I)*COSL
        ZBZL2(I)=ZI*(-1)*(ZI)**(THZ(I)-1)*SR2*GP*ZMIXSS(4,I)*COSL
     $  -(ZI)**(THZ(I))*FL*ZMIXSS(2,I)*SINL
        FL1(I)=ZAZL1(I)*CONJG(ZAZL1(I))+ZBZL1(I)*CONJG(ZBZL1(I))
        GL1(I)=2*REAL(ZAZL1(I)*CONJG(ZBZL1(I)))
        FL2(I)=ZAZL2(I)*CONJG(ZAZL2(I))+ZBZL2(I)*CONJG(ZBZL2(I))
        GL2(I)=2*REAL(ZAZL2(I)*CONJG(ZBZL2(I)))
        SIGZI=SIGZI+FL1(I)*SSB1(QS,MZI(I),ML1)+GL1(I)*MZI(I)/ML*
     $  SSB0(QS,MZI(I),ML1)+FL2(I)*SSB1(QS,MZI(I),ML2)
     $  +GL2(I)*MZI(I)/ML*SSB0(QS,MZI(I),ML2)
      END DO
      RSIGZI=REAL8(SIGZI)
C
C     Chargino-fermion-sfermion couplings; I labels chargino
C
      ZAWL(1)=ZI*(-1)**THW(1)*G*SIN(GAMMAR)
      ZAWL(2)=ZI*(-1)**THW(2)*THY*G*COS(GAMMAR)
      BWPP(1)=-FL*COS(GAMMAL)
      BWPP(2)=FL*THX*SIN(GAMMAL)
      SIGWI=(0.,0.)
      MWI(1)=ABS(AMW1SS)
      MWI(2)=ABS(AMW2SS)
      DO I=1,2
        ZAWSN(I)=-ZI*ZAWL(I)
        ZBWSN(I)=-BWPP(I)
        FWSN(I)=ZAWSN(I)*CONJG(ZAWSN(I))+ZBWSN(I)*CONJG(ZBWSN(I))
        GWSN(I)=2*REAL(CONJG(ZBWSN(I))*ZAWSN(I))
        SIGWI=SIGWI+FWSN(I)*SSB1(QS,MWI(I),MSN)+
     $  MWI(I)/MLQ*GWSN(I)*SSB0(QS,MWI(I),MSN)
      END DO
      RSIGWI=REAL8(SIGWI)
      LP(1)=0.
      LP(2)=.5*FL**2*COSA**2*(REAL8(SSB1(QS,ML,AMHH))+
     $REAL8(SSB0(QS,ML,AMHH)))/FAC
      LP(3)=.5*FL**2*SINA**2*(REAL8(SSB1(QS,ML,AMHL))+
     $REAL8(SSB0(QS,ML,AMHL)))/FAC
      LP(4)=.5*FL**2*SINBE**2*(REAL8(SSB1(QS,ML,AMHA))-
     $REAL8(SSB0(QS,ML,AMHA)))/FAC
      LP(5)=.5*FL**2*COSBE**2*(REAL8(SSB1(QS,ML,MZ))-
     $REAL8(SSB0(QS,ML,MZ)))/FAC
      LP(6)=.5*(FL**2*SINBE**2*REAL8(SSB1(QS,0.,AMHC))+
     $(G2**2+FL**2*COSBE**2)*REAL8(SSB1(QS,0.,MW)))/FAC
      LP(7)=0.
      LP(8)=0.
      LP(9)=+G2**2/COS2W*((GLL**2+GRL**2)*REAL8(SSB1(QS,ML,MZ))
     $+4*GLL*GRL*REAL8(SSB0(QS,ML,MZ)))/FAC
      LP(10)=.5*RSIGZI/FAC
      LP(11)=.5*RSIGWI/FAC
      LPTOT=0.D0
      DO I=1,11
        LPTOT=LPTOT+LP(I)
C        WRITE(6,*) 'LP(',I,')=',LP(I)
      END DO
      SSRSGL=LPTOT
100   RETURN
      END
