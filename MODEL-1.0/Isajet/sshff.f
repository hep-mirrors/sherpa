CDECK  ID>, SSHFF.  
      SUBROUTINE SSHFF(IMODEL)
C-----------------------------------------------------------------------
C
C     Calculate all decays higgs -> f fbar, including QCD radiative
C     corrections for quarks.
C
C     Bisset's SETFAC, WDHFFN, QCDRAD
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
      COMMON/SSLUN/LOUT
      INTEGER LOUT
      SAVE /SSLUN/
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
      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,ASM3
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3
      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG
      SAVE /SUGPAS/
C
      DOUBLE PRECISION PI,SR2,DWID,XG2,MHIH,BETA,BEFAC,ALFAC,MH,MF
     $,MFRUN,FACTOR,ALAM,MF1,MF2,SUM,MF1RUN,MF2RUN,COLOR,TEMP1
     $,QCDFAC
      DOUBLE PRECISION MFIFF(9),MFIF1(6),MFIF2(6)
      DOUBLE PRECISION SSDLAM,SSMQCD,SSHFF1
      REAL WID,COSB,MBMA
      INTEGER IH,IDIH,IFF,IDF,ID1,ID2
      INTEGER IDIFF(9),IDIF1(6),IDIF2(6),IMODEL
C
      PI=4.*ATAN(1.)
      SR2=SQRT(2.)
      BETA=ATAN(1./RV2V1)
      COSB=COS(BETA)
      XG2=4.0*PI*ALFAEM/SN2THW
      MBMA=FBMA*COSB*VEV
C
C          Loop over HL, HH, HA and fermions
C
      MFIFF(1)=AME
      IDIFF(1)=IDE
      MFIFF(2)=AMMU
      IDIFF(2)=IDMU
      MFIFF(3)=AMTAU
      IDIFF(3)=IDTAU
      MFIFF(4)=AMDN
      IDIFF(4)=IDDN
      MFIFF(5)=AMST
      IDIFF(5)=IDST
      MFIFF(6)=AMBT
      IDIFF(6)=IDBT
      MFIFF(7)=AMUP
      IDIFF(7)=IDUP
      MFIFF(8)=AMCH
      IDIFF(8)=IDCH
      MFIFF(9)=AMTP
      IDIFF(9)=IDTP
C
      DO 100 IH=1,3
        IF(IH.EQ.1) THEN
          MH=AMHL
          IDIH=ISHL
          BEFAC=COS(BETA)
          ALFAC=SIN(ALFAH)
        ELSEIF(IH.EQ.2) THEN
          MH=AMHH
          IDIH=ISHH
          BEFAC=COS(BETA)
          ALFAC=COS(ALFAH)
        ELSE
          MH=AMHA
          IDIH=ISHA
          BEFAC=1/TAN(BETA)
          ALFAC=1.
        ENDIF
C
C          Down type fermions
C
        DO 110 IFF=1,6
          MF=MFIFF(IFF)
          IDF=IDIFF(IFF)
          FACTOR=1.-4.*MF**2/MH**2
          IF(FACTOR.LE.0) GO TO 110
          FACTOR=SQRT(FACTOR)
          IF(IFF.GE.4) THEN
             COLOR=3.
             IF (IMODEL.NE.0.AND.IFF.EQ.6.AND.IH.GE.2) THEN
               MFRUN=MBMA
             ELSE
               MFRUN=SSMQCD(MF,MH)
             END IF
             QCDFAC=SSHFF1(MH,MF,IH)
          ELSE
             COLOR=1.
             MFRUN=MF
             QCDFAC=1.
          ENDIF
          DWID=XG2*MFRUN**2*MH*ALFAC**2/(32.*PI*AMW**2*BEFAC**2)
          IF(IH.EQ.1.OR.IH.EQ.2) THEN
            DWID=DWID*FACTOR**3
          ELSEIF(IH.EQ.3) THEN
            DWID=DWID*FACTOR
          ENDIF
          DWID=DWID*COLOR*QCDFAC
          WID=DWID
          CALL SSSAVE(IDIH,WID,IDF,-IDF,0,0,0)
110     CONTINUE
C
C          Up type fermions
C
        IF(IH.EQ.1) THEN
          BEFAC=SIN(BETA)
          ALFAC=COS(ALFAH)
        ELSEIF(IH.EQ.2) THEN
          BEFAC=SIN(BETA)
          ALFAC=SIN(ALFAH)
        ELSE
          BEFAC=TAN(BETA)
          ALFAC=1.
        ENDIF
        DO 120 IFF=7,9
          MF=MFIFF(IFF)
          IDF=IDIFF(IFF)
          FACTOR=1.-4.*MF**2/MH**2
          IF(FACTOR.LE.0) GO TO 120
          FACTOR=SQRT(FACTOR)
          MFRUN=SSMQCD(MF,MH)
          QCDFAC=SSHFF1(MH,MF,IH)
          DWID=XG2*MFRUN**2*MH*ALFAC**2/(32.*PI*AMW**2*BEFAC**2)
          IF(IH.EQ.1.OR.IH.EQ.2) THEN
            DWID=DWID*FACTOR**3
          ELSEIF(IH.EQ.3) THEN
            DWID=DWID*FACTOR
          ENDIF
          DWID=3.*DWID*QCDFAC
          WID=DWID
          CALL SSSAVE(IDIH,WID,IDF,-IDF,0,0,0)
120     CONTINUE
100   CONTINUE
C
C           HC decays. F1 has Iz=+1/2, F2 has Iz=-1/2
C
      MFIF1(1)=0
      IDIF1(1)=IDNE
      MFIF2(1)=AME
      IDIF2(1)=IDE
      MFIF1(2)=0
      IDIF1(2)=IDNM
      MFIF2(2)=AMMU
      IDIF2(2)=IDMU
      MFIF1(3)=0
      IDIF1(3)=IDNT
      MFIF2(3)=AMTAU
      IDIF2(3)=IDTAU
      MFIF1(4)=AMUP
      IDIF1(4)=IDUP
      MFIF2(4)=AMDN
      IDIF2(4)=IDDN
      MFIF1(5)=AMCH
      IDIF1(5)=IDCH
      MFIF2(5)=AMST
      IDIF2(5)=IDST
      MFIF1(6)=AMTP
      IDIF1(6)=IDTP
      MFIF2(6)=AMBT
      IDIF2(6)=IDBT
      MH=AMHC
C
      DO 200 IFF=1,6
        MF1=MFIF1(IFF)
        MF2=MFIF2(IFF)
        ID1=IDIF1(IFF)
        ID2=IDIF2(IFF)
        SUM=MF1+MF2
        ALAM=SSDLAM(MH**2,MF1**2,MF2**2)
        IF(ALAM.LE.0.OR.SUM.GE.MH) GO TO 200
        IF(IFF.LE.3) THEN
          MF1RUN=MF1
          MF2RUN=MF2
          COLOR=1
        ELSE
          MF1RUN=SSMQCD(MF1,MH)
          IF (MF2.EQ.AMBT) THEN
            MF2RUN=MBMA
          ELSE
            MF2RUN=SSMQCD(MF2,MH)
          END IF
          COLOR=3
        ENDIF
        TEMP1=MF1RUN**2*1./TAN(BETA)**2+MF2RUN**2*TAN(BETA)**2
        TEMP1=TEMP1*(MH**2-MF1**2-MF2**2)-4.*MF1**2*MF2**2
        IF (TEMP1.LT.0.0) GO TO 200
        DWID=XG2*COLOR*SQRT(ALAM)*TEMP1/MH**3/(32.0*PI*AMW**2)
        WID=DWID
        CALL SSSAVE(ISHC,WID,ID1,-ID2,0,0,0)
200   CONTINUE
C
      RETURN
      END
