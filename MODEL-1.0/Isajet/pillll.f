CDECK  ID>, PILLLL. 
        REAL FUNCTION PILLLL(P2,G1,G2,G3,CTHW)
C-----------------------------------------------------------------------
C          PILLLL: tau_L slepton self-energy
C     Taken from Damien M. Pierce, Jonathan A. Bagger, Konstantin T. Matchev,
C     Ren-jie Zhang, Nucl.Phys.B491:3-67,1997, hep-ph/9606211
C          P2 = 4-momentum squared
C          CTHW = Cos(theta_W) in DR bar scheme
C-----------------------------------------------------------------------
      IMPLICIT NONE
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
      COMMON/SSINF/XLAM
      DOUBLE PRECISION XLAM
      COMPLEX*16 SSB0,SSG,SSF
      DOUBLE PRECISION SSA0
      COMPLEX TEMP,N(4,4),AC0BBL(4),BC0BBL(4),FLLLL(4),GLLLL(4)
      REAL P2,P,GG,GGP,CTHW,CHWW2,BE,FL,C(4),DD(4),G1,G2,G3
     $,MH0(4),MHP(2),LHLLLI(4,2),LHLLN(2),FNLLL(4),VMAT(2,2)
     $,UMAT(2,2),AP0BBL(4),BP0BBL(4)
     $,APPTBL(2),BPPTBL(2),ACPTBL(2),BCPTBL(2)
     $,LSBLBL(2),LSBLBR(2),LHBLBL(2),LHBLBR(2),LSBLBI(2,2)
     $,LHBLTL(2)
      REAL COST,SINT,COST2,SINT2,COSB,SINB,COSB2,SINB2
     $,COSL,SINL,COSL2,SINL2,THX,THY
      REAL SINA2,COSA2,SINBE2,COSBE2,I3UL,I3DL,I3EL,I3NL,YUR,YUL
     $,YDR,YDL,YER,YEL,YNL,EEL,SWW2,GEL,SINA,COSA,SINBE,COSBE
      INTEGER I,II,THIW1,THIW2
      COMPLEX IMAG
      PARAMETER (IMAG=(0.,1.))
      REAL PI,SR2
      PI=4*ATAN(1.)
      SR2=SQRT(2.)
      P=SQRT(P2)
      CHWW2=CTHW**2
      COST=COS(THETAT)
      SINT=SIN(THETAT)
      COSB=COS(THETAB)
      SINB=SIN(THETAB)
      COSL=COS(THETAL)
      SINL=SIN(THETAL)
      COST2=COST**2
      SINT2=1-COST2
      COSB2=COSB**2
      SINB2=1-COSB2
      COSL2=COSL**2
      SINL2=1-COSL2
      GG=G2
      GGP=SQRT(3./5.)*G1
      BE=ATAN(VUQ/VDQ)
      FL=MLQ/VDQ
      SINA2=SIN(ALFAH)**2
      COSA2=COS(ALFAH)**2
      SINBE2=SIN(BE)**2
      COSBE2=COS(BE)**2
      SINA=SIN(ALFAH)
      COSA=COS(ALFAH)
      SINBE=SIN(BE)
      COSBE=COS(BE)
      C(1)=-(COSA2-SINA2)
      C(2)=COSA2-SINA2
      C(3)=-(COSBE2-SINBE2)
      C(4)=COSBE2-SINBE2
      DD(1)=COSA2
      DD(2)=SINA2
      DD(3)=COSBE2
      DD(4)=SINBE2
      MH0(1)=MSS(30)
      MH0(2)=MSS(29)
      MH0(3)=AMZ
      MH0(4)=MSS(31)
      MHP(1)=AMW
      MHP(2)=MSS(32)
      I3UL=1./2.
      I3DL=-1./2.
      I3EL=-1./2.
      I3NL=1./2.
      YUR=-4./3.
      YUL=1./3.
      YDR=2./3.
      YDL=1./3.
      YER=2.
      YEL=-1.
      YNL=-1.
      EEL=-1.
      SWW2=1-(AMW/AMZ)**2
      GEL=I3EL-EEL*SWW2
      THX=SIGN(1.,TAN(GAMMAL))
      THY=SIGN(1.,TAN(GAMMAR))
      IF (SIGN(1.,AMW1SS).EQ.1.) THEN
         THIW1=0
      ELSE
         THIW1=1
      END IF
      IF (SIGN(1.,AMW2SS).EQ.1.) THEN
         THIW2=0
      ELSE
         THIW2=1
      END IF
      VMAT(1,1)=-(-1)**THIW1*SIN(GAMMAR)
      VMAT(1,2)=(-1)**THIW1*COS(GAMMAR)
      VMAT(2,1)=-(-1)**THIW2*THY*COS(GAMMAR)
      VMAT(2,2)=-(-1)**THIW2*THY*SIN(GAMMAR)
      UMAT(1,1)=-SIN(GAMMAL)
      UMAT(1,2)=COS(GAMMAL)
      UMAT(2,1)=-THX*COS(GAMMAL)
      UMAT(2,2)=-THX*SIN(GAMMAL)
      DO II=1,4
        IF (SIGN(1.,AMZISS(II)).EQ.1.) THEN
          I=0
        ELSE
          I=1
        END IF
        N(II,1)=-IMAG**I*ZMIXSS(4,II)
        N(II,2)=-IMAG**I*ZMIXSS(3,II)
        N(II,3)=IMAG**I*ZMIXSS(2,II)
        N(II,4)=IMAG**I*ZMIXSS(1,II)
      ENDDO
      DO I=1,4
        AP0BBL(I)=0.
        BP0BBL(I)=0.
        AC0BBL(I)=0.
        BC0BBL(I)=0.
      ENDDO
      BP0BBL(1)=GGP/SR2*YEL
      BP0BBL(2)=SR2*GG*I3EL
      AP0BBL(3)=FL
      DO I=1,4
        DO II=1,4
          AC0BBL(I)=AC0BBL(I)+CONJG(N(I,II))*AP0BBL(II)
          BC0BBL(I)=BC0BBL(I)+N(I,II)*BC0BBL(II)
        ENDDO
      ENDDO
      DO I=1,4
        FLLLL(I)=CONJG(AC0BBL(I))*AC0BBL(I)
     $+CONJG(BC0BBL(I))*BC0BBL(I)
        GLLLL(I)=CONJG(BC0BBL(I))*AC0BBL(I)
     $+CONJG(AC0BBL(I))*BC0BBL(I)
      ENDDO
      DO I=1,2
        APPTBL(I)=0.
        BPPTBL(I)=0.
        ACPTBL(I)=0.
        BCPTBL(I)=0.
      ENDDO
      BPPTBL(1)=GG
      DO I=1,2
        DO II=1,2
          ACPTBL(I)=ACPTBL(I)+VMAT(I,II)*APPTBL(II)
          BCPTBL(I)=BCPTBL(I)+UMAT(I,II)*BPPTBL(II)
        ENDDO
      ENDDO
      DO I=1,2
        FNLLL(I)=ACPTBL(I)*ACPTBL(I)+BCPTBL(I)*BCPTBL(I)
      ENDDO
      LSBLBL(1)=GG*AMZ/CTHW*GEL*COSBE+SR2*FL*MLQ
      LSBLBR(1)=-FL/SR2*AAL
      LSBLBL(2)=-GG*AMZ/CTHW*GEL*SINBE
      LSBLBR(2)=-FL/SR2*TWOM1
C     First index is for Higgs sector, second is for sfermion
      LSBLBI(1,1)=LSBLBL(1)*COSL-LSBLBR(1)*SINL
      LSBLBI(1,2)=LSBLBL(1)*SINL+LSBLBR(1)*COSL
      LSBLBI(2,1)=LSBLBL(2)*COSL-LSBLBR(2)*SINL
      LSBLBI(2,2)=LSBLBL(2)*SINL+LSBLBR(2)*COSL
      LHLLLI(1,1)=COSA*LSBLBI(1,1)-SINA*LSBLBI(2,1)
      LHLLLI(1,2)=COSA*LSBLBI(1,2)-SINA*LSBLBI(2,2)
      LHLLLI(2,1)=SINA*LSBLBI(1,1)+COSA*LSBLBI(2,1)
      LHLLLI(2,2)=SINA*LSBLBI(1,2)+COSA*LSBLBI(2,2)
      LHBLBL(1)=0.
      LHBLBR(1)=-FL/SR2*(-TWOM1*SINBE-AAL*COSBE)
      LHBLBL(2)=0.
      LHBLBR(2)=-FL/SR2*(-TWOM1*COSBE+AAL*SINBE)
      LHLLLI(3,1)=LHBLBL(1)*COSL-LHBLBR(1)*SINL
      LHLLLI(3,2)=LHBLBL(1)*SINL+LHBLBR(1)*COSL
      LHLLLI(4,1)=LHBLBL(2)*COSL-LHBLBR(2)*SINL
      LHLLLI(4,2)=LHBLBL(2)*SINL+LHBLBR(2)*COSL
      LHBLTL(1)=-GG*AMW/SR2*(COSBE2-SINBE2)+FL*MLQ*COSBE
      LHBLTL(2)=GG*AMW/SR2*2.*COSBE*SINBE-FL*MLQ*SINBE
      LHLLN(1)=LHBLTL(1)
      LHLLN(2)=LHBLTL(2)
      TEMP=FL**2*(SINL2*SSA0(MSS(21))+COSL2*SSA0(MSS(22)))
     $+4.*GG**2/CHWW2*GEL**2*SSA0(AMZ)+2.*GG**2*SSA0(AMW)
     $+4.*PI/137.036*(COSL2*SSF(P,MSS(21),0.)
     $+SINL2*SSF(P,MSS(22),0.))+GG**2/CHWW2*GEL**2
     $*(COSL2*SSF(P,MSS(21),AMZ)+SINL2*SSF(P,MSS(22),AMZ))
     $+GG**2/2.*SSF(P,MSS(16),AMW)
     $+GG**2/4.*(COSL2*SSA0(MSS(21))+SINL2*SSA0(MSS(22))
     $+2.*SSA0(MSS(16)))
     $+GG**2*I3EL*(3.*(I3UL*(SSA0(MSS(2))+SSA0(MSS(8))
     $+COST2*SSA0(MSS(13))+SINT2*SSA0(MSS(12)))
     $+I3DL*(SSA0(MSS(4))+SSA0(MSS(6))
     $+COSB2*SSA0(MSS(10))+SINB2*SSA0(MSS(11))))
     $+I3EL*(SSA0(MSS(17))+SSA0(MSS(19))
     $+COSL2*SSA0(MSS(21))+SINL2*SSA0(MSS(22)))
     $+I3NL*(SSA0(MSS(14))+SSA0(MSS(15))+SSA0(MSS(16))))
     $+GGP**2/4.*YEL**2*(COSL2*SSA0(MSS(21))+SINL2*SSA0(MSS(22)))
     $+GGP**2/4.*YEL*(3.*(YUL*(SSA0(MSS(2))+SSA0(MSS(8))
     $+COST2*SSA0(MSS(13))+SINT2*SSA0(MSS(12)))
     $+YUR*(SSA0(MSS(3))+SSA0(MSS(9))+SINT2*SSA0(MSS(13))
     $+COST2*SSA0(MSS(12)))+YDL*(SSA0(MSS(4))+SSA0(MSS(6))
     $+COSB2*SSA0(MSS(10))+SINB2*SSA0(MSS(11)))
     $+YDR*(SSA0(MSS(5))+SSA0(MSS(7))
     $+SINB2*SSA0(MSS(10))+COSB2*SSA0(MSS(11))))
     $+YEL*(SSA0(MSS(17))+SSA0(MSS(19))
     $+COSL2*SSA0(MSS(21))+SINL2*SSA0(MSS(22)))
     $+YER*(SSA0(MSS(18))+SSA0(MSS(20))
     $+SINL2*SSA0(MSS(21))+COSL2*SSA0(MSS(22)))
     $+YNL*(SSA0(MSS(14))+SSA0(MSS(15))+SSA0(MSS(16))))
      DO I=1,4
        TEMP=TEMP+(FL**2*DD(I)-GG**2*GEL/2./CHWW2*C(I))
     $*SSA0(MH0(I))/2.
      ENDDO
      DO I=3,4
        TEMP=TEMP+GG*2*(GEL/2./CHWW2-I3EL)*C(I)
     $*SSA0(MHP(I-2))
      ENDDO
      DO I=1,4
          TEMP=TEMP+(LHLLLI(I,1))**2*SSB0(P2,MH0(I),MSS(21))
     $+(LHLLLI(I,2))**2*SSB0(P2,MH0(I),MSS(22))
      ENDDO
      DO I=1,2
          TEMP=TEMP+(LHLLN(I))**2*SSB0(P2,MSS(16),MHP(I))
      ENDDO
      DO I=1,4
        TEMP=TEMP+FLLLL(I)*SSG(P,ABS(MSS(22+I)),AMTAU)
     $-2.*GLLLL(I)*ABS(MSS(22+I))*MLQ*SSB0(P2,ABS(MSS(22+I)),AMTAU)
      ENDDO
      DO I=1,2
        TEMP=TEMP+FNLLL(I)*SSG(P,ABS(MSS(26+I)),0.)
      ENDDO
      PILLLL=REAL(TEMP)/16./PI**2
      RETURN
      END
