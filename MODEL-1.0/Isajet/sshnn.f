CDECK  ID>, SSHNN.  
      SUBROUTINE SSHNN
C-----------------------------------------------------------------------
C     Calculates the decay widths of all neutral Higgses into all
C     possible pairs of neutralinos, and the decay widths of the
C     charged Higgs into any neutralino and any chargino
C
C     Bisset's NEUINO
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
      DOUBLE PRECISION XIJ,XJI,DIJ,TEMP,DWZN,TEMP2,T2,RWZ,SWZ
     $,PI,SR2,XM,THETX,YM,THETY,SGL,CGL,SGR,CGR,MW1,MW2,THETM,THETP
     $,G2,GP2,BETA,ALPHA,MH,M1,M2
      DOUBLE PRECISION SN1,SN2,DWID,LAMB
      DOUBLE PRECISION A(4,4),MHI(3)
      DOUBLE PRECISION SSDLAM
      REAL WID
      INTEGER II,NUMH,I1,I2,IZ,IW,ID1,ID2,IDHHA
      INTEGER IDHI(3),IDZI(4)
C
C          Mass matrix parameters
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
C          The following was calculated in Bisset's MASZIN
      DO 10 II=1,4
         TEMP=SQRT(G2)*ZMIXSS(3,II)+SQRT(GP2)*ZMIXSS(4,II)
         TEMP=TEMP/SR2
         A(1,II)=-TEMP*SGR-SQRT(G2)*ZMIXSS(1,II)*CGR
         A(2,II)=TEMP*CGR-SQRT(G2)*ZMIXSS(1,II)*SGR
         A(3,II)=-TEMP*SGL+SQRT(G2)*ZMIXSS(2,II)*CGL
         A(4,II)=TEMP*CGL+SQRT(G2)*ZMIXSS(2,II)*SGL
10    CONTINUE
C
C          Arrays for loops
C
      MHI(1)=AMHL
      MHI(2)=AMHH
      MHI(3)=AMHA
      IDHI(1)=ISHL
      IDHI(2)=ISHH
      IDHI(3)=ISHA
      IDZI(1)=ISZ1
      IDZI(2)=ISZ2
      IDZI(3)=ISZ3
      IDZI(4)=ISZ4
C
C          Loop over neutral Higgs decays h(numh) into neutralino
C          pairs zi(i1) and zi(i2)
C
      DO 100 NUMH=1,3
        MH=MHI(NUMH)
        IDHHA=IDHI(NUMH)
        DO 110 I1=1,4
          M1=ABS(AMZISS(I1))
          ID1=IDZI(I1)
          DO 120 I2=I1,4
            M2=ABS(AMZISS(I2))
            ID2=IDZI(I2)
            IF(M1+M2.GE.MH) GO TO 120
            LAMB=SSDLAM(MH**2,M1**2,M2**2)
            IF(I1.EQ.I2) THEN
              DIJ = 0.5
            ELSE
              DIJ = 1.0
            ENDIF
            TEMP=-0.5*SIGN(1.,AMZISS(I1))*SIGN(1.,AMZISS(I2))
            XIJ=TEMP*(SQRT(G2)*ZMIXSS(3,I2)-SQRT(GP2)*ZMIXSS(4,I2))
            XJI=TEMP*(SQRT(G2)*ZMIXSS(3,I1)-SQRT(GP2)*ZMIXSS(4,I1))
            IF(NUMH.EQ.1) THEN
              XIJ=XIJ*(ZMIXSS(2,I1)*SIN(ALPHA)-ZMIXSS(1,I1)*COS(ALPHA))
              XJI=XJI*(ZMIXSS(2,I2)*SIN(ALPHA)-ZMIXSS(1,I2)*COS(ALPHA))
            ELSEIF (NUMH .EQ. 2) THEN
              XIJ=XIJ*(ZMIXSS(2,I1)*COS(ALPHA)+ZMIXSS(1,I1)*SIN(ALPHA))
              XJI=XJI*(ZMIXSS(2,I2)*COS(ALPHA)+ZMIXSS(1,I2)*SIN(ALPHA))
            ELSEIF(NUMH.EQ.3) THEN
              XIJ=-XIJ*(ZMIXSS(2,I1)*SIN(BETA)-ZMIXSS(1,I1)*COS(BETA))
              XJI=-XJI*(ZMIXSS(2,I2)*SIN(BETA)-ZMIXSS(1,I2)*COS(BETA))
            ENDIF
            DWID=DIJ*(XIJ+XJI)**2
            DWID=DWID*SQRT(LAMB)/(8.0*PI*(MH**3))
            IF(NUMH.EQ.1.OR.NUMH.EQ.2) THEN
              TEMP2 = ((MH**2)-(M1-2.0*TEMP*M2)**2)
            ELSEIF(NUMH.EQ.3) THEN
              TEMP2=((MH**2)-(M1+2.0*TEMP*M2)**2)
            ENDIF
            DWID=DWID*TEMP2
            WID=DWID
            CALL SSSAVE(IDHHA,WID,ID1,ID2,0,0,0)
120       CONTINUE
110     CONTINUE
100   CONTINUE
C
C          Loop over h+ decays into wi(iw) + zi(iz)
C
      MH=AMHC
      DO 210 IW=1,2
        IF(IW.EQ.1) THEN
          M1=ABS(AMW1SS)
          ID1=ISW1
          SN1=SIGN(1.,AMW1SS)
        ELSE
          M1=ABS(AMW2SS)
          ID1=ISW2
          SN1=SIGN(1.,AMW2SS)
        ENDIF
        DO 220 IZ=1,4
          M2=ABS(AMZISS(IZ))
          ID2=IDZI(IZ)
          SN2=SIGN(1.,AMZISS(IZ))
          IF(M1+M2.GE.MH) GO TO 220
          LAMB=SSDLAM(MH**2,M1**2,M2**2)
          T2=MH**2-M1**2-M2**2
          IF(IW.EQ.1) THEN
            RWZ=COS(BETA)*A(2,IZ)*SN1
            TEMP=SIN(BETA)*A(4,IZ)*SN2
            SWZ=0.5*(RWZ+TEMP)
            RWZ=0.5*(RWZ-TEMP)
          ELSE
            RWZ=COS(BETA)*A(1,IZ)*THETY*SN1
            TEMP=SIN(BETA)*A(3,IZ)*THETX*SN2
            SWZ=0.5*(RWZ+TEMP)
            RWZ=0.5*(RWZ-TEMP)
          ENDIF
          DWID=RWZ**2+SWZ**2
          DWID=DWID*T2
          TEMP=2*M1*M2*(RWZ**2-SWZ**2)
          DWID=(DWID-TEMP)/(8.0*PI*(MH**3))
          DWID=DWID*SQRT(LAMB)
          WID=DWID
          CALL SSSAVE(ISHC,WID,ID1,ID2,0,0,0)
220     CONTINUE
210   CONTINUE
      RETURN
      END
