CDECK  ID>, SUGMAS. 
C---------------------------------------------------------------
      SUBROUTINE SUGMAS(G0,ILOOP,IMODEL)
C---------------------------------------------------------------
C
C     Compute tree level sparticle masses; output to MSS, XISAIN
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
      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,ASM3
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3
      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG
      SAVE /SUGPAS/
C     XSUGIN contains the inputs to SUGRA:
C     XSUGIN(1) = M_0        XSUGIN(2) = M_(1/2)  XSUGIN(3) = A_0
C     XSUGIN(4) = tan(beta)  XSUGIN(5) = sgn(mu)  XSUGIN(6) = M_t
C     XSUGIN(7) = SUG BC scale
C     XGMIN(1) = LAM         XGMIN(2)  = M_MES    XGMIN(3)  = XN5
C     XGMIN(4) = tan(beta)   XGMIN(5)  = sgn(mu)  XGMIN(6) = M_t
C     XGMIN(7) = CGRAV       XGMIN(8)  =RSL       XGMIN(9)  = DEL_HD
C     XGMIN(10)  = DEL_HU    XGMIN(11) = DY       XGMIN(12) = N5_1
C     XGMIN(13)  = N5_2      XGMIN(14) = N5_3
C     XNRIN(1) = M_N3        XNRIN(2) = M_MAJ     XNRIN(3) = ANSS
C     XNRIN(4) = M_N3SS
C     XISAIN contains the MSSMi inputs in natural order.
      COMMON /SUGXIN/ XISAIN(24),XSUGIN(7),XGMIN(14),XNRIN(4)
      REAL XISAIN,XSUGIN,XGMIN,XNRIN
      SAVE /SUGXIN/
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
      REAL MSB1,MSB2,MST1,MST2
      REAL G0(31)
      REAL SUGMFN,SUALFS,SSPOLE,MHP,MGLMGL,MHPS,
     $RDEL,ASMGL,DELHPS,M1S,M2S,FNB,FCN,
     $MB,FNT,MT,MW,TANB,BETA,COSB,COTB,SINB,MZ,COS2B,
     $PI,T2S,G,AL,MSSS,AT,AB,BRKT,B2S,T1S,TERM,B1S,Q,
     $FNL,MSL1,MSL2,COS2W,MAJXS,MAJX,
     $SIGST,SIGT,SIGSB,SIGB,SIGSL,SIGL,RZT,RZB,RZL
      REAL AA,BB,CC,DA,DB,DC,L1,L2,EVAL1,RL1,RL2
      DOUBLE PRECISION SSMQCD
      INTEGER IALLOW,ILOOP,MHLNEG,MHCNEG,IMODEL
C
C          Statement function
C
      SUGMFN(Q)=Q**2*(LOG(Q**2/HIGFRZ**2)-1.)
C
      MHPNEG=0
      PI=4.*ATAN(1.)
      XW=.232
      G=G0(2)
      COS2W=1.-SN2THW
      TANB=VUQ/VDQ
      MT=AMT
      MZ=AMZ
      MW=AMW
      AMTP=MT
      BETA=ATAN(TANB)
      COTB=1./TANB
      SINB=SIN(BETA)
      COSB=COS(BETA)
      SIN2B=SIN(2*BETA)
      COS2B=COS(2*BETA)
      AT=G0(12)
      AB=G0(11)
      AL=G0(10)
      MLQ=G0(4)*VDQ
      MBQ=G0(5)*VDQ
      MTQ=G0(6)*VUQ
C
C          Compute some masses from RGE solution to prepare for SSMASS,
C          which computes the rest.
C
      MSSS=G0(19)+AMUP**2+(.5-2*XW/3.)*MZ**2*COS2B
      IF (MSSS.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
C          Squark and slepton masses
      MSS(2)=SQRT(MSSS)
      MSS(3)=SQRT(G0(18)+AMUP**2+2./3.*XW*MZ**2*COS2B)
      MSS(4)=SQRT(G0(19)+AMDN**2+(-.5+XW/3.)*MZ**2*COS2B)
      MSS(5)=SQRT(G0(17)+AMDN**2-1./3.*XW*MZ**2*COS2B)
      MSS(6)=SQRT(G0(19)+AMST**2+(-.5+XW/3.)*MZ**2*COS2B)
      MSS(7)=SQRT(G0(17)+AMST**2-1./3.*XW*MZ**2*COS2B)
      MSS(8)=SQRT(G0(19)+AMCH**2+(.5-2*XW/3.)*MZ**2*COS2B)
      MSS(9)=SQRT(G0(18)+AMCH**2+2./3.*XW*MZ**2*COS2B)
      BRKT=(.5*(G0(24)-G0(22))-COS2B*(4*MW**2-MZ**2)/12.)**2+
     $       MBQ**2*(AB-MU*TANB)**2
      TERM=.5*(G0(24)+G0(22))+MBQ**2-MZ**2*COS2B/4.
      B1S=TERM-SQRT(BRKT)
      B2S=TERM+SQRT(BRKT)
      MSS(10)=SQRT(MAX(0.,B1S))
      MSS(11)=SQRT(MAX(0.,B2S))
      THETAB=ATAN((B1S-MBQ**2+MZ**2*COS2B*(.5-XW/3.)-
     $G0(24))/MBQ/(AB-MU*TANB))
      BRKT=(.5*(G0(24)-G0(23))+COS2B*(8*MW**2-5*MZ**2)/12.)**2+
     $       MTQ**2*(AT-MU*COTB)**2
      TERM=.5*(G0(24)+G0(23))+MTQ**2+MZ**2*COS2B/4.
      T1S=TERM-SQRT(BRKT)
      IF (T1S.LE.0..OR.B1S.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
      T2S=TERM+SQRT(BRKT)
      MSS(12)=SQRT(MAX(0.,T1S))
      MSS(13)=SQRT(MAX(0.,T2S))
      THETAT=ATAN((T1S-MTQ**2+MZ**2*COS2B*(-.5+2*XW/3.)-
     $G0(24))/MTQ/(AT-MU*COTB))
      MSSS=G0(16)+.5*MZ**2*COS2B
      IF (MSSS.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
      MSS(14)=SQRT(MSSS)
      MSS(15)=MSS(14)
      MSSS=G0(21)+.5*MZ**2*COS2B
      IF (MSSS.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
      MSS(16)=SQRT(MSSS)
      MSSS=G0(16)+AME**2-.5*(2*MW**2-MZ**2)*COS2B
      IF (MSSS.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
      MSS(17)=SQRT(MSSS)
      MSSS=G0(15)+AME**2+(MW**2-MZ**2)*COS2B
      IF (MSSS.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
      MSS(18)=SQRT(MSSS)
      MSS(19)=SQRT(G0(16)+AMMU**2-.5*(2*MW**2-MZ**2)*COS2B)
      MSS(20)=SQRT(G0(15)+AMMU**2+(MW**2-MZ**2)*COS2B)
      BRKT=(.5*(G0(21)-G0(20))-COS2B*(4*MW**2-3*MZ**2)/4.)**2+
     $       MLQ**2*(AL-MU*TANB)**2
      TERM=.5*(G0(21)+G0(20))+MLQ**2-MZ**2*COS2B/4.
      T1S=TERM-SQRT(BRKT)
      IF (T1S.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
      T2S=TERM+SQRT(BRKT)
      MSS(21)=SQRT(MAX(0.,T1S))
      MSS(22)=SQRT(MAX(0.,T2S))
      THETAL=ATAN((T1S-MLQ**2+MZ**2*COS2B*(.5-XW)-
     $G0(21))/MLQ/(AL-MU*TANB))

C          A0 mass
      M1S=MU**2+G0(13)
      M2S=MU**2+G0(14)
      MSB1=MSS(10)
      MSB2=MSS(11)
      MST1=MSS(12)
      MST2=MSS(13)
      MSL1=MSS(21)
      MSL2=MSS(22)
      MB=AMBT
      FNT=(SUGMFN(MST2)-SUGMFN(MST1))/(MST2**2-MST1**2)
     $*AT*MTQ**2/SINB**2
      FNB=(SUGMFN(MSB2)-SUGMFN(MSB1))/(MSB2**2-MSB1**2)
     $*AB*MBQ**2/COSB**2
      FNL=(SUGMFN(MSL2)-SUGMFN(MSL1))/(MSL2**2-MSL1**2)
     $*AL*MLQ**2/COSB**2
      FCN=FNT+FNB+FNL/3.
C      DELHPS=3*G0(2)**2*MU*(COTB+TANB)/32./PI**2/MW**2*FCN
      DELHPS=6*MU*(COTB+TANB)/32./PI**2/VEV**2*FCN
      RDEL=SQRT(ABS(DELHPS))
C     Try Javier/Xerxes improved mA formula
      RZT=.5*(G0(24)-G0(23))+MZ**2*(8*COS2W-5.)*COS2B/12.
      RZB=.5*(G0(24)-G0(22))-MZ**2*(4*COS2W-1.)*COS2B/12.
      RZL=.5*(G0(21)-G0(20))-MZ**2*(4*COS2W-3.)*COS2B/4.
C      ZT=.25*(MST2**2-MST1**2)**2-(G0(6)*VEV)**2*(AT*SINB-MU*COSB)**2
C      ZB=.25*(MSB2**2-MSB1**2)**2-(G0(5)*VEV)**2*(AB*COSB-MU*SINB)**2
C      ZL=.25*(MSL2**2-MSL1**2)**2-(G0(4)*VEV)**2*(AL*COSB-MU*SINB)**2
      SIGT=3*G0(6)**2*SUGMFN(MTQ)/8./PI/PI
      SIGB=-3*G0(5)**2*SUGMFN(MBQ)/8./PI/PI
      SIGL=-1*G0(4)**2*SUGMFN(MLQ)/8./PI/PI
      SIGST=-3.*((SUGMFN(MST1)+SUGMFN(MST2))*(G0(6)**2-.25*
     $(G**2+GP**2))+(SUGMFN(MST2)-SUGMFN(MST1))/(MST2**2-MST1**2)*
     $(G0(6)**2*(AT**2-MU**2)-(G**2+GP**2)*(8*COS2W-5.)*
     $RZT/6.))/16./PI/PI
      SIGSB=3.*((SUGMFN(MSB1)+SUGMFN(MSB2))*(G0(5)**2-.25*
     $(G**2+GP**2))+(SUGMFN(MSB2)-SUGMFN(MSB1))/(MSB2**2-MSB1**2)*
     $(G0(5)**2*(AB**2-MU**2)-(G**2+GP**2)*(4*COS2W-1.)*
     $RZB/6.))/16./PI/PI
      SIGSL=1.*((SUGMFN(MSL1)+SUGMFN(MSL2))*(G0(4)**2-.25*
     $(G**2+GP**2))+(SUGMFN(MSL2)-SUGMFN(MSL1))/(MSL2**2-MSL1**2)*
     $(G0(4)**2*(AL**2-MU**2)-(G**2+GP**2)*(2*COS2W-1.5)*
     $RZL/6.))/16./PI/PI
C          Tree level mhp not needed at this point so fix if negative
      IF (ILOOP.EQ.0) THEN
        MHPS=M1S+M2S
        IF (MHPS.LT.0.) MHPS=0.
      ELSE
C        MHPS=B*MU*(COTB+TANB)+DELHPS
C        Use improved Javier formula
        MHPS=(TANB**2+1.)/(TANB**2-1.)*(G0(13)-G0(14)+SIGT+SIGST+
     $SIGSB+SIGB+SIGSL+SIGL)-AMZ**2
        IF (MHPS.LT.0.) THEN
          MHPNEG=1
          MHPS=1.
        END IF
      END IF
      MHP=SQRT(MHPS)
C     WRITE(6,*) 'MHP=',MHP,' MAJX=',MAJX
      MSS(31)=MHP
C     APPLY XERXES' TEST FOR PROPER POTENTIAL SHAPE AT THE ORIGIN
C     REMOVE THIS CONSTRAINT ON 4/7/00
C      IF (ILOOP.EQ.1) THEN
C      L1=MIN(G0(24),G0(23))
C      L2=MAX(G0(24),G0(23))
C      RL1=SQRT(L1)
C      RL2=SQRT(L2)
C      DA=3*G0(6)**2*AT**2/ABS(G0(24)-G0(23))/16./PI**2*
C     $(-SUGMFN(RL1)+SUGMFN(RL2))
C      DB=3*G0(6)**2/16./PI**2*
C     $(SUGMFN(RL1)*(1.-AT**2/ABS(G0(24)-G0(23)))+SUGMFN(RL2)*
C     $(1.+AT**2/ABS(G0(24)-G0(23))))
C      DC=-3*G0(6)**2*AT*MU/ABS(G0(24)-G0(23))/16./PI**2*
C     $(-SUGMFN(RL1)+SUGMFN(RL2))
C      AA=M1S+DA
C      BB=M2S+DB
C      CC=-B*MU+DC
C      EVAL1=((AA+BB)-SQRT((AA+BB)**2-4*(AA*BB-CC*CC)))/2.
C      IF (EVAL1.GE.0) THEN
C        NOGOOD=7
C      END IF
C      END IF
C
C          Initialize SUSY parameters in /SSPAR/:
C
      AMGLSS=ABS(G0(9))
      AMULSS=MSS(2)
      AMURSS=MSS(3)
      AMDLSS=MSS(4)
      AMDRSS=MSS(5)
      AMSLSS=MSS(6)
      AMSRSS=MSS(7)
      AMCLSS=MSS(8)
      AMCRSS=MSS(9)
      AMN1SS=MSS(14)
      AMN2SS=MSS(15)
      AMN3SS=MSS(16)
      AMELSS=MSS(17)
      AMERSS=MSS(18)
      AMMLSS=MSS(19)
      AMMRSS=MSS(20)
      TWOM1=-MU
      RV2V1=1./XTANB
      AMTLSS=SIGN(1.,G0(24))*SQRT(ABS(G0(24)))
      AMTRSS=SIGN(1.,G0(23))*SQRT(ABS(G0(23)))
      AMBLSS=SQRT(G0(24))
      AMBRSS=SQRT(G0(22))
      AMLLSS=SQRT(G0(21))
      AMLRSS=SQRT(G0(20))
      AMB1SS=MSS(10)
      AMB2SS=MSS(11)
      AMT1SS=MSS(12)
      AMT2SS=MSS(13)
      AML1SS=MSS(21)
      AML2SS=MSS(22)
      AAT=G0(12)
      AAB=G0(11)
      AAL=G0(10)
      AMHA=MHP
C
C          Use SSMASS to diagonalize neutralino and chargino mass
C          matrices and calculate Higgs masses.
C
      MHLNEG=0
      MHCNEG=0
      CALL SSMASS(G0(9),G0(7),G0(8),IALLOW,ILOOP,MHLNEG,MHCNEG,IMODEL)
      IF(MHLNEG.EQ.1.OR.MHCNEG.EQ.1) THEN
        NOGOOD=8
      ENDIF
C      IF(IALLOW.NE.0.AND.ILOOP.NE.0) THEN
C        NOGOOD=5
C        GO TO 100
C      ENDIF
C
C          Save results also in MSS; re-save radiative corrected masses
C
      MSS(2)=AMULSS
      MSS(3)=AMURSS
      MSS(4)=AMDLSS
      MSS(5)=AMDRSS
      MSS(6)=AMSLSS
      MSS(7)=AMSRSS
      MSS(8)=AMCLSS
      MSS(9)=AMCRSS
      MSS(10)=AMB1SS
      MSS(11)=AMB2SS
      MSS(12)=AMT1SS
      MSS(13)=AMT2SS
      MSS(14)=AMN1SS
      MSS(15)=AMN2SS
      MSS(16)=AMN3SS
      MSS(17)=AMELSS
      MSS(18)=AMERSS
      MSS(19)=AMMLSS
      MSS(20)=AMMRSS
      MSS(21)=AML1SS
      MSS(22)=AML2SS
      MSS(23)=AMZ1SS
      MSS(24)=AMZ2SS
      MSS(25)=AMZ3SS
      MSS(26)=AMZ4SS
      MSS(27)=AMW1SS
      MSS(28)=AMW2SS
      MSS(29)=AMHL
      MSS(30)=AMHH
      MSS(31)=AMHA
      MSS(32)=AMHC
C          Gluino pole mass
      MGLMGL=G0(9)
      ASMGL=SUALFS(MGLMGL**2,.36,MT,3)
      MSS(1)=SIGN(1.,G0(9))*SSPOLE(ABS(MGLMGL),MGLMGL**2,ASMGL)
      AMGLSS=ABS(MSS(1))
      GSS(9)=G0(9)
C
100   RETURN
      END
