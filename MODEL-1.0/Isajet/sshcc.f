CDECK  ID>, SSHCC.  
      SUBROUTINE SSHCC
C-----------------------------------------------------------------------
C     Calculates the decay widths of all neutral Higgses into all
C     possible pairs of charginos.
C
C     Bisset's CHGINO
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
      DOUBLE PRECISION PI,SR2,XM,THETX,YM,THETY,SGL,CGL,SGR,CGR
     $,MW1,MW2,THETM,THETP,G2,GP2,BETA,ALPHA,T1,MH,M1,M2,LAMB
     $,DWID,TEMP,TEMPXY
      DOUBLE PRECISION MHI(3),IDHI(3),SHP(3),SHM(3),SH(3),PH(3)
      DOUBLE PRECISION SSDLAM
      REAL WID
      INTEGER NUMH,IDHHA
C
      PI=4.*ATAN(1.D0)
      SR2=SQRT(2.D0)
      XM=1./TAN(GAMMAL)
      THETX=SIGN(1.D0,XM)
      YM=1./TAN(GAMMAR)
      THETY=SIGN(1.D0,YM)
      SGL=1/(DSQRT(1+XM**2))
      CGL=SGL*XM
      SGR=1/(DSQRT(1+YM**2))
      CGR=SGR*YM
      MW1=DBLE(ABS(AMW1SS))
      MW2=DBLE(ABS(AMW2SS))
      THETM=SIGN(1.,AMW1SS)
      THETP=SIGN(1.,AMW2SS)
      G2=4*PI*ALFAEM/SN2THW
      GP2=4*PI*ALFAEM/(1-SN2THW)
      BETA=ATAN(1.0/RV2V1)
      ALPHA=ALFAH
C
C          Arrays for loops
C
      MHI(1)=AMHL
      MHI(2)=AMHH
      MHI(3)=AMHA
      IDHI(1)=ISHL
      IDHI(2)=ISHH
      IDHI(3)=ISHA
C          The following came from Bisset's MASZIN, but with L,H,P
C          replaced by a generic H and a subscript.
      TEMPXY=0.5*THETX*THETY*(-THETP)
      SHP(1)=SIN(ALPHA)*CGR*SGL+COS(ALPHA)*CGL*SGR
      SHP(1)=SHP(1)*TEMPXY
      SHM(1)=SIN(ALPHA)*SGR*CGL+COS(ALPHA)*SGL*CGR
      SHM(1)=SHM(1)*0.5*THETM
      SH(1)=-THETX*SGR*SGL*SIN(ALPHA)*THETM
      PH(1)=-SH(1)
      T1=THETX*CGL*CGR*COS(ALPHA)*THETM
      SH(1)=SH(1)+T1
      PH(1)=PH(1)-T1
      T1=THETY*SGL*SGR*COS(ALPHA)*THETP
      SH(1)=SH(1)-T1
      PH(1)=PH(1)-T1
      T1=THETY*CGL*CGR*SIN(ALPHA)*THETP
      SH(1)=SH(1)+T1
      PH(1)=PH(1)+T1
      SH(1)=0.5*SH(1)
      PH(1)=0.5*PH(1)
      SHP(2)=COS(ALPHA)*CGR*SGL-SIN(ALPHA)*CGL*SGR
      SHP(2)=SHP(2)*TEMPXY
      SHM(2)=COS(ALPHA)*SGR*CGL-SIN(ALPHA)*SGL*CGR
      SHM(2)=SHM(2)*0.5*THETM
      SH(2)=-THETX*SGR*SGL*COS(ALPHA)*THETM
      PH(2)=-SH(2)
      T1=THETX*CGL*CGR*SIN(ALPHA)*THETM
      SH(2)=SH(2)-T1
      PH(2)=PH(2)+T1
      T1=THETY*SGL*SGR*SIN(ALPHA)*THETP
      SH(2)=SH(2)+T1
      PH(2)=PH(2)+T1
      T1=THETY*CGL*CGR*COS(ALPHA)*THETP
      SH(2)=SH(2)+T1
      PH(2)=PH(2)+T1
      SH(2)=0.5*SH(2)
      PH(2)=0.5*PH(2)
      SHP(3)=SIN(BETA)*CGR*SGL+COS(BETA)*CGL*SGR
      SHP(3)=SHP(3)*0.5*THETX*THETY*(-THETP)
      SHM(3)=SIN(BETA)*SGR*CGL+COS(BETA)*SGL*CGR
      SHM(3)=SHM(3)*0.5*THETM
      SH(3)=-THETX*SGR*SGL*SIN(BETA)*THETM
      PH(3)=-SH(3)
      T1=THETX*CGL*CGR*COS(BETA)*THETM
      SH(3)=SH(3)+T1
      PH(3)=PH(3)-T1
      T1=THETY*SGL*SGR*COS(BETA)*THETP
      SH(3)=SH(3)+T1
      PH(3)=PH(3)+T1
      T1=THETY*CGL*CGR*SIN(BETA)*THETP
      SH(3)=SH(3)-T1
      PH(3)=PH(3)-T1
      SH(3)=0.5*SH(3)
      PH(3)=0.5*PH(3)
C
C          Loop over neutral Higgs
C
      DO 100 NUMH=1,3
        MH=MHI(NUMH)
        IDHHA=IDHI(NUMH)
C          w1 + w1
        M1=ABS(AMW1SS)
        M2=M1
        IF(MH.GT.M1+M2) THEN
          LAMB=SSDLAM(MH**2,M1**2,M2**2)
          TEMP=1-4*M1**2/MH**2
          DWID=G2*MH*SHM(NUMH)**2/(4.0*PI)
          DWID=DWID*SQRT(TEMP**3)
          WID=DWID
          CALL SSSAVE(IDHHA,WID,ISW1,-ISW1,0,0,0)
        ENDIF
C          w2 + w2
        M1=ABS(AMW2SS)
        M2=M1
        IF(MH.GT.M1+M2) THEN
          TEMP=1-4*M1**2/MH**2
          DWID=G2*MH*SHP(NUMH)**2/(4*PI)
          DWID=DWID*SQRT(TEMP**3)
          WID=DWID
          CALL SSSAVE(IDHHA,WID,ISW2,-ISW2,0,0,0)
        ENDIF
C          w1 + w2
        M1=ABS(AMW1SS)
        M2=ABS(AMW2SS)
        IF(MH.GT.M1+M2) THEN
          LAMB=SSDLAM(MH**2,M1**2,M2**2)
          DWID=PH(NUMH)**2*(MH**2-(M1-M2)**2)
          DWID=DWID+SH(NUMH)**2*(MH**2-(M1+M2)**2)
          DWID=DWID*G2*SQRT(LAMB)/(16.0*PI*(MH**3))
          WID=DWID
          CALL SSSAVE(IDHHA,WID,ISW1,-ISW2,0,0,0)
          CALL SSSAVE(IDHHA,WID,-ISW1,ISW2,0,0,0)
        ENDIF
100   CONTINUE
C
      RETURN
      END
