CDECK  ID>, SSHWW.  
      SUBROUTINE SSHWW
C-----------------------------------------------------------------------
C     Calculate HL, HH -> WW, ZZ, using either the on-shell matrix
C     element if kinematically allowed or the WW* or ZZ* matrix
C     element from Eqn.(6) for Keung and Marciano (PRD. 84: 248).
C     For the latter, save the mode as W(Z) f fbar, and require that
C     MH > MW + 2 * MB.
C
C     Bisset's GBDCY
C-----------------------------------------------------------------------
      IMPLICIT NONE
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
C          Temporary parameters for functions
      COMMON/SSTMP/TMP(10),ITMP(10)
      REAL TMP
      INTEGER ITMP
      SAVE /SSTMP/
C          SM ident code definitions. These are standard ISAJET but
C          can be changed.
      INTEGER IDUP,IDDN,IDST,IDCH,IDBT,IDTP
      INTEGER IDNE,IDE,IDNM,IDMU,IDNT,IDTAU
      INTEGER IDGL,IDGM,IDW,IDZ,IDH
      PARAMETER (IDUP=1,IDDN=2,IDST=3,IDCH=4,IDBT=5,IDTP=6)
      PARAMETER (IDNE=11,IDE=12,IDNM=13,IDMU=14,IDNT=15,IDTAU=16)
      PARAMETER (IDGL=9,IDGM=10,IDW=80,IDZ=90,IDH=81)
C          SUSY ident code definitions. They are chosen to be similar
C          to those in versions < 6.50 but may be changed.
      INTEGER ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1
      INTEGER ISNEL,ISEL,ISNML,ISMUL,ISNTL,ISTAU1
      INTEGER ISUPR,ISDNR,ISSTR,ISCHR,ISBT2,ISTP2
      INTEGER ISNER,ISER,ISNMR,ISMUR,ISNTR,ISTAU2
      INTEGER ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,ISGL
      INTEGER ISHL,ISHH,ISHA,ISHC
      INTEGER ISGRAV
      INTEGER IDTAUL,IDTAUR
      PARAMETER (ISUPL=21,ISDNL=22,ISSTL=23,ISCHL=24,ISBT1=25,ISTP1=26)
      PARAMETER (ISNEL=31,ISEL=32,ISNML=33,ISMUL=34,ISNTL=35,ISTAU1=36)
      PARAMETER (ISUPR=41,ISDNR=42,ISSTR=43,ISCHR=44,ISBT2=45,ISTP2=46)
      PARAMETER (ISNER=51,ISER=52,ISNMR=53,ISMUR=54,ISNTR=55,ISTAU2=56)
      PARAMETER (ISGL=29)
      PARAMETER (ISZ1=30,ISZ2=40,ISZ3=50,ISZ4=60,ISW1=39,ISW2=49)
      PARAMETER (ISHL=82,ISHH=83,ISHA=84,ISHC=86)
      PARAMETER (ISGRAV=91)
      PARAMETER (IDTAUL=10016,IDTAUR=20016)
C
      EXTERNAL SSHWW1,SSHWW2
      DOUBLE PRECISION SSHWW1,SSHWW2
      DOUBLE PRECISION PI,SR2,G2,BETA,ALPHA,SW2,CW2,CAB2,SAB2,MW,MZ
     $,MH,COUPL,LOWER,UPPER,FWW1,FWW2,FWW3,FWW,DWID,FZZ
      DOUBLE PRECISION SSDINT,SSDLAM
      REAL WID
      REAL BRZN,BRZL,BRZU,BRZD,BRWL,BRWQ
      INTEGER IDHHA,IH
C          Hard wired Z branching ratios
      DATA BRZN,BRZL,BRZU,BRZD/.06839,.03442,.11792,.15191/
      DATA BRWL,BRWQ/.11111,.33333/
C
C          Mass matrix parameters
C
      PI=4*ATAN(1.D0)
      SR2=SQRT(2.D0)
      G2=4*PI*ALFAEM/SN2THW
      BETA=ATAN(1.0/RV2V1)
      ALPHA=ALFAH
      SW2=SN2THW
      CW2=1.-SN2THW
      CAB2=(DCOS(ALPHA+BETA))**2
      SAB2=1.0-CAB2
      MW=AMW
      MZ=AMZ
C
C          WW* and ZZ* decays
C
      DO 100 IH=1,2
        IF(IH.EQ.1) THEN
          MH=AMHL
          IDHHA=ISHL
          COUPL=SAB2
        ELSE
          MH=AMHH
          IDHHA=ISHH
          COUPL=CAB2
        ENDIF
C          H -> W + W* -> W + f + fbar
        TMP(1)=MH
        IF(MH.GT.MW+2*AMBT.AND.MH.LE.2*MW) THEN
          LOWER=2*MW/MH
          UPPER=1+MW**2/MH**2
          IF (LOWER.LT.0.998D0) THEN
            IF (UPPER.LE.1.001D0) THEN
              FWW1=SSDINT(LOWER,SSHWW1,0.998D0)
              FWW2=SSDINT(0.998D0,SSHWW1,UPPER)
              FWW=FWW1+FWW2
            ELSEIF(UPPER.GT.1.001D0) THEN
              FWW1=SSDINT(LOWER,SSHWW1,0.998D0)
              FWW2=SSDINT(0.998D0,SSHWW1,1.001D0)
              FWW3=SSDINT(1.001D0,SSHWW1,UPPER)
              FWW=FWW1+FWW2+FWW3
            ENDIF
          ELSE IF (0.998D0.LT.LOWER.AND.LOWER.LT.1.001D0) THEN
            IF (UPPER.LE.1.001D0) THEN
              FWW=SSDINT(LOWER,SSHWW1,UPPER)
            ELSEIF(UPPER.GT.1.001D0) THEN
              FWW1=SSDINT(LOWER,SSHWW1,1.001D0)
              FWW2=SSDINT(1.001D0,SSHWW1,UPPER)
              FWW=FWW1+FWW2
            ENDIF
          ELSE IF (LOWER.GT.1.001D0) THEN
            FWW=SSDINT(LOWER,SSHWW1,UPPER)
          END IF
          DWID=3*(G2**2)*MH*FWW/(512.0*PI**3)
          WID=DWID*COUPL
          CALL SSSAVE(IDHHA,0.5*BRWL*WID,IDW,IDE,-IDNE,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,0.5*BRWL*WID,IDW,IDMU,-IDNM,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,0.5*BRWL*WID,IDW,IDTAU,-IDNT,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,0.5*BRWQ*WID,IDW,-IDUP,IDDN,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,0.5*BRWQ*WID,IDW,-IDCH,IDST,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,0.5*BRWL*WID,-IDW,-IDE,IDNE,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,0.5*BRWL*WID,-IDW,-IDMU,IDNM,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,0.5*BRWL*WID,-IDW,-IDTAU,IDNT,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,0.5*BRWQ*WID,-IDW,IDUP,-IDDN,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,0.5*BRWQ*WID,-IDW,IDCH,-IDST,0,0)
          CALL SSSVME(9)
        ENDIF
C          H -> Z + Z* -> Z + f + fbar
        IF(MH.GT.MZ+2*AMBT.AND.MH.LE.2*MZ) THEN
          LOWER=2*MZ/MH
          UPPER=1+MZ**2/MH**2
          FZZ=SSDINT(LOWER,SSHWW2,UPPER)
          DWID=7.0-40*SW2/3+160*SW2**2/9
          DWID=DWID/CW2**2
          DWID=DWID*G2**2*MH*FZZ/(2048*PI**3)
          WID=DWID*COUPL
          CALL SSSAVE(IDHHA,BRZN*WID,IDZ,IDNE,-IDNE,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,BRZN*WID,IDZ,IDNM,-IDNM,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,BRZN*WID,IDZ,IDNT,-IDNT,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,BRZL*WID,IDZ,IDE,-IDE,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,BRZL*WID,IDZ,IDMU,-IDMU,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,BRZL*WID,IDZ,IDTAU,-IDTAU,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,BRZU*WID,IDZ,IDUP,-IDUP,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,BRZU*WID,IDZ,IDCH,-IDCH,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,BRZD*WID,IDZ,IDDN,-IDDN,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,BRZD*WID,IDZ,IDST,-IDST,0,0)
          CALL SSSVME(9)
          CALL SSSAVE(IDHHA,BRZD*WID,IDZ,IDBT,-IDBT,0,0)
          CALL SSSVME(9)
        ENDIF
100   CONTINUE
C
C          HH -> WW, ZZ
C          If these are allowed, the WW* and ZZ* are not.
C
      MH=AMHH
      IF(MH.GT.2*MW) THEN
        DWID=3+(MH/MW)**4/4-(MH/MW)**2
        DWID=DWID*G2*CAB2*MW**2/(16.0*PI*MH**3)
        WID=DWID*SQRT(SSDLAM(MH**2,MW**2,MW**2))
        CALL SSSAVE(ISHH,WID,IDW,-IDW,0,0,0)
      ENDIF
      IF(MH.GT.2*MZ) THEN
        DWID=3+(MH/MZ)**4/4-(MH/MZ)**2
        DWID=DWID*G2*CAB2*MW**2/(16.0*PI*MH**3)/(2.0*CW2**2)
        WID=DWID*SQRT(SSDLAM(MH**2,MZ**2,MZ**2))
        CALL SSSAVE(ISHH,WID,IDZ,IDZ,0,0,0)
      ENDIF
C
      RETURN
      END
