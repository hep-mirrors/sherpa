CDECK  ID>, SSTPBF.
        SUBROUTINE SSTPBF
C-----------------------------------------------------------------------
C
C     Calculate the top branching ratios.
C     Source: H. Baer (modified by F. Paige)
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
      COMMON/SSLUN/LOUT
      INTEGER LOUT
      SAVE /SSLUN/
C          MXSS         =  maximum number of modes
C          NSSMOD       = number of modes
C          ISSMOD       = initial particle
C          JSSMOD       = final particles
C          GSSMOD       = width
C          BSSMOD       = branching ratio
C          MSSMOD       = decay matrix element pointer
C          LSSMOD       = logical flag used internally by SSME3
      INTEGER MXSS
      PARAMETER (MXSS=1000)
      COMMON/SSMODE/NSSMOD,ISSMOD(MXSS),JSSMOD(5,MXSS),GSSMOD(MXSS)
     $,BSSMOD(MXSS),MSSMOD(MXSS),LSSMOD
      INTEGER NSSMOD,ISSMOD,JSSMOD,MSSMOD
      REAL GSSMOD,BSSMOD
      LOGICAL LSSMOD
      SAVE /SSMODE/
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
      COMMON/SSPAR/AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS(4,4)
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS
      REAL AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS
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
      PARAMETER (ISUPL=21,ISDNL=22,ISSTL=23,ISCHL=24,ISBT1=25,ISTP1=26)
      PARAMETER (ISNEL=31,ISEL=32,ISNML=33,ISMUL=34,ISNTL=35,ISTAU1=36)
      PARAMETER (ISUPR=41,ISDNR=42,ISSTR=43,ISCHR=44,ISBT2=45,ISTP2=46)
      PARAMETER (ISNER=51,ISER=52,ISNMR=53,ISMUR=54,ISNTR=55,ISTAU2=56)
      PARAMETER (ISGL=29)
      PARAMETER (ISZ1=30,ISZ2=40,ISZ3=50,ISZ4=60,ISW1=39,ISW2=49)
      PARAMETER (ISHL=82,ISHH=83,ISHA=84,ISHC=86)
      PARAMETER (ISGRAV=91)
C
      COMPLEX ZI,ZONE,ZA,ZB,ZPP,ZPM,ZAUIZ,ZBUIZ
      REAL SSXLAM,G,AL2,BE2,TANB,COTB,GTBW,GTBH,BWLL,GF,BWQQ,PI,SR2
      REAL WID,AS,BS,MZIZ,CS2THW,GP,FT,FB,SNZI,THIZ
      REAL SINT,COST,SINB,COSB,AWI,BWI,AMW1,AMW2,SNWI
      REAL THX,THY,XM,YM,BETA,ASMB,MBMB,MBQ,ASMT,MTMT,MTQ,SUALFS
      DOUBLE PRECISION SSMQCD
      INTEGER IZ,ISZIZ(4)
      DATA ZONE/(1.,0.)/,ZI/(0.,1.)/
C
      PI=4.*ATAN(1.)
      SR2=SQRT(2.)
      G=SQRT(4*PI*ALFAEM/SN2THW)
      GP=G*SQRT(SN2THW/(1.-SN2THW))
      TANB=1./RV2V1
      COTB=1./TANB
      BETA=ATAN(TANB)
      CS2THW=1.-SN2THW
      ASMB=SUALFS(AMBT**2,.36,AMTP,3)
      MBMB=AMBT*(1.-4*ASMB/3./PI)
      MBQ=SSMQCD(DBLE(MBMB),DBLE(AMTP))
      ASMT=SUALFS(AMTP**2,.36,AMTP,3)
      MTMT=AMTP/(1.+4*ASMT/3./PI+(16.11-1.04*(5.-6.63/AMTP))*
     $(ASMT/PI)**2)
      MTQ=SSMQCD(DBLE(MTMT),DBLE(AMTP))
      FB=G*MBQ/SR2/AMW/COS(BETA)
      FT=G*MTQ/SR2/AMW/SIN(BETA)
      SINT=SIN(THETAT)
      COST=COS(THETAT)
      SINB=SIN(THETAB)
      COSB=COS(THETAB)
      ISZIZ(1)=ISZ1
      ISZIZ(2)=ISZ2
      ISZIZ(3)=ISZ3
      ISZIZ(4)=ISZ4
      XM=1./TAN(GAMMAL)
      YM=1./TAN(GAMMAR)
      THX=SIGN(1.,XM)
      THY=SIGN(1.,YM)
C
C          W decays
C
      GF=1.16E-5
      GTBW=GF*AMTP**3*SSXLAM(1.,AMW**2/AMTP**2,AMBT**2/AMTP**2)*
     $((1.-AMBT**2/AMTP**2)**2+AMW**2/AMTP**2*(1.+AMBT**2/AMTP**2)
     $-2*AMW**4/AMTP**4)/(8.*PI*SR2)
      BWQQ=3./9.
      BWLL=1./9.
      CALL SSSAVE(IDTP,BWQQ*GTBW,IDUP,-IDDN,IDBT,0,0)
      CALL SSSAVE(IDTP,BWQQ*GTBW,IDCH,-IDST,IDBT,0,0)
      CALL SSSAVE(IDTP,BWLL*GTBW,-IDE,IDNE,IDBT,0,0)
      CALL SSSAVE(IDTP,BWLL*GTBW,-IDMU,IDNM,IDBT,0,0)
      CALL SSSAVE(IDTP,BWLL*GTBW,-IDTAU,IDNT,IDBT,0,0)
C
C          H+ decays
C
      AL2=(G/2/SR2/AMW*(AMBT*TANB+AMTP*COTB))**2
      BE2=(G/2/SR2/AMW*(AMBT*TANB-AMTP*COTB))**2
      IF (AMTP.GT.(AMBT+AMHC)) THEN
        GTBH=AMTP/16./PI*((AL2+BE2)
     $  *(1.+AMBT**2/AMTP**2-AMHC**2/AMTP**2)
     $  +2*(AL2-BE2)*AMBT/AMTP)
     $  *SQRT(SSXLAM(1.,AMHC**2/AMTP**2,AMBT**2/AMTP**2))
        CALL SSSAVE(IDTP,GTBH,ISHC,IDBT,0,0,0)
      END IF
C
C          t->t_1 + z_i decays
      DO 100 IZ=1,4
        MZIZ=ABS(AMZISS(IZ))
        SNZI=SIGN(1.,AMZISS(IZ))
        IF (SNZI.EQ.1.) THEN
           THIZ=0.
        ELSE
           THIZ=1.
        END IF
        ZAUIZ=ZI**(THIZ-1.)*SNZI*
     $(-G/SR2*ZMIXSS(3,IZ)-GP/3./SR2*ZMIXSS(4,IZ))
        ZBUIZ=ZI**(THIZ-1.)*4*GP*ZMIXSS(4,IZ)/3./SR2
        ZPP=ZI**THIZ
        ZPM=(-ZI)**THIZ
        ZA=((ZI*ZAUIZ-ZPP*FT*ZMIXSS(1,IZ))*COST-
     $(ZI*ZBUIZ-ZPM*FT*ZMIXSS(1,IZ))*SINT)/2.
        ZB=((-ZI*ZAUIZ-ZPP*FT*ZMIXSS(1,IZ))*COST-
     $(ZI*ZBUIZ+ZPM*FT*ZMIXSS(1,IZ))*SINT)/2.
        AS=ZA*CONJG(ZA)
        BS=ZB*CONJG(ZB)
        IF (AMTP.GT.(AMT1SS+MZIZ)) THEN
          WID=(AS*((AMTP+MZIZ)**2-AMT1SS**2)+BS*
     $((AMTP-MZIZ)**2-AMT1SS**2))/16./PI/AMTP*
     $SQRT(SSXLAM(1.,AMT1SS**2/AMTP**2,MZIZ**2/AMTP**2))
          CALL SSSAVE(IDTP,WID,ISZIZ(IZ),ISTP1,0,0,0)
        END IF
100   CONTINUE
C
C       t -> sb_1 + sW_i
C
        AMW1=ABS(AMW1SS)
        AMW2=ABS(AMW2SS)
        IF (AMTP.GT.(AMB1SS+AMW1)) THEN
          SNWI=SIGN(1.,AMW1SS)
          AWI=-G*SIN(GAMMAL)*COSB+FB*COS(GAMMAL)*SINB
          BWI=-FT*(-SNWI)*COS(GAMMAR)
          WID=AMTP*((AWI**2+BWI**2*COSB**2)*(1.+AMW1**2/AMTP**2
     $-AMB1SS**2/AMTP**2)+4*AMW1/AMTP*AWI*BWI*COST)/32./PI*
     $SQRT(SSXLAM(1.,AMW1**2/AMTP**2,AMB1SS**2/AMTP**2))
          CALL SSSAVE(IDTP,WID,ISW1,ISBT1,0,0,0)
        END IF
c
        IF (AMTP.GT.(AMB1SS+AMW2)) THEN
          SNWI=SIGN(1.,AMW2SS)
          AWI=-G*THX*COS(GAMMAL)*COSB-FB*THX*SIN(GAMMAL)*SINB
          BWI=FT*(-SNWI)*THY*SIN(GAMMAR)
          WID=AMTP*((AWI**2+BWI**2*COSB**2)*(1.+AMW2**2/AMTP**2
     $-AMB1SS**2/AMTP**2)+4*AMW2/AMTP*AWI*BWI*COST)/32./PI*
     $SQRT(SSXLAM(1.,AMW2**2/AMTP**2,AMB1SS**2/AMTP**2))
          CALL SSSAVE(IDTP,WID,ISW2,ISBT1,0,0,0)
        END IF
C
C       t -> sb_2 + sW_i
C
        IF (AMTP.GT.(AMB2SS+AMW1)) THEN
          SNWI=SIGN(1.,AMW1SS)
          AWI=-G*SIN(GAMMAL)*SINB-FB*COS(GAMMAL)*COSB
          BWI=-FT*(-SNWI)*COS(GAMMAR)
          WID=AMTP*((AWI**2+BWI**2*SINB**2)*(1.+AMW1**2/AMTP**2
     $-AMB2SS**2/AMTP**2)+4*AMW1/AMTP*AWI*BWI*COST)/32./PI*
     $SQRT(SSXLAM(1.,AMW1**2/AMTP**2,AMB2SS**2/AMTP**2))
          CALL SSSAVE(IDTP,WID,ISW1,ISBT2,0,0,0)
        END IF
c
        IF (AMTP.GT.(AMB2SS+AMW2)) THEN
          SNWI=SIGN(1.,AMW2SS)
          AWI=-G*THX*COS(GAMMAL)*SINB+FB*THX*SIN(GAMMAL)*COSB
          BWI=FT*(-SNWI)*THY*SIN(GAMMAR)
          WID=AMTP*((AWI**2+BWI**2*SINB**2)*(1.+AMW2**2/AMTP**2
     $-AMB2SS**2/AMTP**2)+4*AMW2/AMTP*AWI*BWI*COST)/32./PI*
     $SQRT(SSXLAM(1.,AMW2**2/AMTP**2,AMB2SS**2/AMTP**2))
          CALL SSSAVE(IDTP,WID,ISW2,ISBT2,0,0,0)
        END IF
C
C
C          Normalize branching ratios
C
      CALL SSNORM(IDTP)
C
      RETURN
      END
