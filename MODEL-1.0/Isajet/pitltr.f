CDECK  ID>, PITLTR. 
        REAL FUNCTION PITLTR(P2,G1,G2,G3,CTHW)
C-----------------------------------------------------------------------
C          PITLTR: off-diagonal t squark self-energy
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
      COMPLEX TEMP,N(4,4),AC0TTR(4),AC0TTL(4),BC0TTR(4),BC0TTL(4)
     $,FTTLR(4),GTTLR(4)
      REAL P2,P,GG,GGP,CTHW,CHWW2,BE,FT,FB,G1,G2,G3
     $,MH0(4),MHP(2),LHTLTI(4,2),LHTLBI(2,2),FBTLR(4),VMAT(2,2)
     $,UMAT(2,2),GBTLR(4),AP0TTR(4),AP0TTL(4),BP0TTR(4),BP0TTL(4)
     $,APPBTL(2),APPBTR(2),BPPBTL(2),BPPBTR(2)
     $,ACPBTL(2),ACPBTR(2),BCPBTL(2),BCPBTR(2)
     $,LSTLTL(2),LSTLTR(2),LHTLTL(2),LHTLTR(2),LSTLTI(2,2)
     $,LHTLBL(2),LHTLBR(2)
     $,LHTRTI(4,2),LHTRBI(2,2),LSTRTR(2),LHTRTR(2),LSTRTI(2,2)
     $,LHTRBR(2),LHTRBL(2)
      REAL COST,SINT,COST2,SINT2,COSB,SINB,COSB2,SINB2
     $,COSL,SINL,COSL2,SINL2,THX,THY
      REAL SINA2,COSA2,SINBE2,COSBE2,I3UL,YUR,YUL
     $,EUL,EUR,SWW2,GUL,GUR,SINA,COSA,SINBE,COSBE
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
      FT=MTQ/VUQ
      FB=MBQ/VDQ
      SINA2=SIN(ALFAH)**2
      COSA2=COS(ALFAH)**2
      SINBE2=SIN(BE)**2
      COSBE2=COS(BE)**2
      SINA=SIN(ALFAH)
      COSA=COS(ALFAH)
      SINBE=SIN(BE)
      COSBE=COS(BE)
      MH0(1)=MSS(30)
      MH0(2)=MSS(29)
      MH0(3)=AMZ
      MH0(4)=MSS(31)
      MHP(1)=AMW
      MHP(2)=MSS(32)
      I3UL=1./2.
      YUR=-4./3.
      YUL=1./3.
      EUL=2./3.
      EUR=-2./3.
      SWW2=1-(AMW/AMZ)**2
      EUR=-2./3.
      GUR=-EUR*SWW2
      GUL=I3UL-EUL*SWW2
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
        AP0TTR(I)=0.
        AP0TTL(I)=0.
        BP0TTR(I)=0.
        BP0TTL(I)=0.
        AC0TTR(I)=0.
        AC0TTL(I)=0.
        BC0TTR(I)=0.
        BC0TTL(I)=0.
      ENDDO
      AP0TTR(1)=GGP/SR2*YUR
      BP0TTL(1)=GGP/SR2*YUL
      BP0TTL(2)=SR2*GG*I3UL
      AP0TTL(4)=FT
      BP0TTR(4)=FT
      DO I=1,4
        DO II=1,4
          AC0TTR(I)=AC0TTR(I)+CONJG(N(I,II))*AP0TTR(II)
          AC0TTL(I)=AC0TTL(I)+CONJG(N(I,II))*AP0TTL(II)
          BC0TTR(I)=BC0TTR(I)+N(I,II)*BC0TTR(II)
          BC0TTL(I)=BC0TTL(I)+N(I,II)*BC0TTL(II)
        ENDDO
      ENDDO
      DO I=1,4
        FTTLR(I)=CONJG(AC0TTL(I))*AC0TTR(I)
     $+CONJG(BC0TTL(I))*BC0TTR(I)
        GTTLR(I)=CONJG(BC0TTL(I))*AC0TTR(I)
     $+CONJG(AC0TTL(I))*BC0TTR(I)
      ENDDO
      DO I=1,2
        APPBTL(I)=0.
        APPBTR(I)=0.
        BPPBTL(I)=0.
        BPPBTR(I)=0.
        ACPBTL(I)=0.
        ACPBTR(I)=0.
        BCPBTL(I)=0.
        BCPBTR(I)=0.
      ENDDO
      APPBTL(1)=GG
      BPPBTL(2)=-FB
      APPBTR(2)=-FT
      DO I=1,2
        DO II=1,2
          ACPBTL(I)=ACPBTL(I)+VMAT(I,II)*APPBTL(II)
          ACPBTR(I)=ACPBTR(I)+VMAT(I,II)*APPBTR(II)
          BCPBTL(I)=BCPBTL(I)+UMAT(I,II)*BPPBTL(II)
          BCPBTR(I)=BCPBTR(I)+UMAT(I,II)*BPPBTR(II)
        ENDDO
      ENDDO
      DO I=1,2
        FBTLR(I)=ACPBTL(I)*ACPBTR(I)+BCPBTL(I)*BCPBTR(I)
        GBTLR(I)=BCPBTL(I)*ACPBTR(I)+ACPBTL(I)*BCPBTR(I)
      ENDDO
      LSTLTL(1)=GG*AMZ/CTHW*GUL*COSBE
      LSTLTR(1)=-FT/SR2*TWOM1
      LSTLTL(2)=-GG*AMZ/CTHW*GUL*SINBE+SR2*FT*MTQ
      LSTLTR(2)=-FT/SR2*AAT
C     First index is for Higgs sector, second is for sfermion
      LSTLTI(1,1)=LSTLTL(1)*COST-LSTLTR(1)*SINT
      LSTLTI(1,2)=LSTLTL(1)*SINT+LSTLTR(1)*COST
      LSTLTI(2,1)=LSTLTL(2)*COST-LSTLTR(2)*SINT
      LSTLTI(2,2)=LSTLTL(2)*SINT+LSTLTR(2)*COST
      LHTLTI(1,1)=COSA*LSTLTI(1,1)-SINA*LSTLTI(2,1)
      LHTLTI(1,2)=COSA*LSTLTI(1,2)-SINA*LSTLTI(2,2)
      LHTLTI(2,1)=SINA*LSTLTI(1,1)+COSA*LSTLTI(2,1)
      LHTLTI(2,2)=SINA*LSTLTI(1,2)+COSA*LSTLTI(2,2)
      LHTLTL(1)=0.
      LHTLTR(1)=FT/SR2*(-TWOM1*COSBE-AAT*SINBE)
      LHTLTL(2)=0.
      LHTLTR(2)=-FT/SR2*(-TWOM1*SINBE+AAT*COSBE)
      LHTLTI(3,1)=LHTLTL(1)*COST-LHTLTR(1)*SINT
      LHTLTI(3,2)=LHTLTL(1)*SINT+LHTLTR(1)*COST
      LHTLTI(4,1)=LHTLTL(2)*COST-LHTLTR(2)*SINT
      LHTLTI(4,2)=LHTLTL(2)*SINT+LHTLTR(2)*COST
      LHTLBL(1)=-GG*AMW/SR2*(COSBE2-SINBE2)-FT*MTQ*SINBE
     $+FB*MBQ*COSBE
      LHTLBR(1)=FB*(-TWOM1*SINBE-AAB*COSBE)
      LHTLBL(2)=GG*AMW/SR2*2.*COSBE*SINBE-FT*MTQ*COSBE
     $-FB*MBQ*SINBE
      LHTLBR(2)=FB*(-TWOM1*COSBE+AAB*SINBE)
      LHTLBI(1,1)=LHTLBL(1)*COSB-LHTLBR(1)*SINB
      LHTLBI(1,2)=LHTLBL(1)*SINB+LHTLBR(1)*COSB
      LHTLBI(2,1)=LHTLBL(2)*COSB-LHTLBR(2)*SINB
      LHTLBI(2,2)=LHTLBL(2)*SINB+LHTLBR(2)*COSB
      LSTRTR(1)=GG*AMZ/CTHW*GUR*COSBE
      LSTRTR(2)=-GG*AMZ/CTHW*GUR*SINBE+SR2*FT*MTQ
      LSTRTI(1,1)=LSTLTR(1)*COST-LSTRTR(1)*SINT
      LSTRTI(1,2)=LSTLTR(1)*SINT+LSTRTR(1)*COST
      LSTRTI(2,1)=LSTLTR(2)*COST-LSTRTR(2)*SINT
      LSTRTI(2,2)=LSTLTR(2)*SINT+LSTRTR(2)*COST
      LHTRTI(1,1)=COSA*LSTRTI(1,1)-SINA*LSTRTI(2,1)
      LHTRTI(1,2)=COSA*LSTRTI(1,2)-SINA*LSTRTI(2,2)
      LHTRTI(2,1)=SINA*LSTRTI(1,1)+COSA*LSTRTI(2,1)
      LHTRTI(2,2)=SINA*LSTRTI(1,2)+COSA*LSTRTI(2,2)
      LHTRTR(1)=0.
      LHTRTR(2)=0.
      LHTRTI(3,1)=-LHTLTR(1)*COST-LHTRTR(1)*SINT
      LHTRTI(3,2)=-LHTLTR(1)*SINT+LHTRTR(1)*COST
      LHTRTI(4,1)=-LHTLTR(2)*COST-LHTRTR(2)*SINT
      LHTRTI(4,2)=-LHTLTR(2)*SINT+LHTRTR(2)*COST
      LHTRBR(1)=0.
      LHTRBL(1)=-FT*(-TWOM1*COSBE-AAT*SINBE)
      LHTRBR(2)=-FT*MBQ*COSBE-FB*MTQ*SINBE
      LHTRBL(2)=FB*(-TWOM1*SINBE+AAT*COSBE)
      LHTRBI(1,1)=LHTRBL(1)*COSB-LHTRBR(1)*SINB
      LHTRBI(1,2)=LHTRBL(1)*SINB+LHTRBR(1)*COSB
      LHTRBI(2,1)=LHTRBL(2)*COSB-LHTRBR(2)*SINB
      LHTRBI(2,2)=LHTRBL(2)*SINB+LHTRBR(2)*COSB
      TEMP=4.*G3**2/3.*(4.*MSS(1)*MTQ*SSB0(P2,MSS(1),AMTP)
     $-SINT*COST*(SSF(P,MSS(13),0.)-SSF(P,MSS(12),0.)
     $-SSA0(MSS(13))+SSA0(MSS(12))))
     $-FT/2.*3.*FT*2.*SINT*COST*(SSA0(MSS(13))-SSA0(MSS(12)))
     $-GGP**2/4.*YUL*YUR*SINT*COST*(SSA0(MSS(13))-SSA0(MSS(12)))
     $-(2./3.)**2*4.*PI/137.036*SINT*COST*(SSF(P,MSS(13),0.)
     $-SSF(P,MSS(12),0.))+GG**2/CHWW2*GUL*GUR*SINT*COST
     $*(SSF(P,MSS(13),AMZ)-SSF(P,MSS(12),AMZ))
      DO I=1,4
          TEMP=TEMP+LHTLTI(I,1)*LHTRTI(I,1)
     $*SSB0(P2,MH0(I),MSS(13))+LHTLTI(I,2)*LHTRTI(I,2)
     $*SSB0(P2,MH0(I),MSS(12))
      ENDDO
      DO I=1,2
        DO II=1,2
          TEMP=TEMP+LHTLBI(I,II)*LHTRBI(I,II)*SSB0(P2,MSS(9+II),MHP(I))
        ENDDO
      ENDDO
      DO I=1,4
        TEMP=TEMP+FTTLR(I)*SSG(P,ABS(MSS(22+I)),AMTP)
     $-2.*GTTLR(I)*ABS(MSS(22+I))*MTQ*SSB0(P2,ABS(MSS(22+I)),AMTP)
      ENDDO
      DO I=1,2
        TEMP=TEMP+FBTLR(I)*SSG(P,ABS(MSS(26+I)),AMBT)
     $-2.*GBTLR(I)*ABS(MSS(26+I))*MBQ*SSB0(P2,ABS(MSS(26+I)),AMBT)
      ENDDO
      PITLTR=REAL(TEMP)/16./PI**2
      RETURN
      END
