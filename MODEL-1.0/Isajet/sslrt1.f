CDECK  ID>, SSLRT1. 
        REAL FUNCTION SSLRT1(SS)
C-----------------------------------------------------------------------
C          SSLRT1: l_R -> l+tau+stau_1
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
C          Temporary parameters for functions
      COMMON/SSTMP/TMP(10),ITMP(10)
      REAL TMP
      INTEGER ITMP
      SAVE /SSTMP/
        REAL SS
        DOUBLE PRECISION ETMX,ETMN,S,MT1,MT,MLR,BEZI,BEZJ,TM,
     ,AI,AJ,BI,BJ,MZI,MZJ,SNZI,SNZJ,XL,BK1,BK2,BK3,BK,WID,SSDLAM
        S=SS
        MT1=AML1SS
        MT=AMTAU
        MLR=TMP(1)
        BEZI=TMP(2)
        BEZJ=TMP(3)
        AI=TMP(4)
        AJ=TMP(5)
        BI=TMP(6)
        BJ=TMP(7)
        MZI=ABS(TMP(8))
        MZJ=ABS(TMP(9))
        SNZI=SIGN(1.0,TMP(8))
        SNZJ=SIGN(1.0,TMP(9))
        TM=SSDLAM(S,MT**2,MT1**2)
        XL=DSQRT(MAX(0.D0,TM))
        ETMN=(S+MT**2-MT1**2-XL*(MLR**2-S)/(MLR**2+S))*(MLR**2+S)/
     ,       (2*S)/(2*MLR)
        ETMX=(S+MT**2-MT1**2+XL*(MLR**2-S)/(MLR**2+S))*(MLR**2+S)/
     ,       (2*S)/(2*MLR)
        BK1=-(ETMX-ETMN)*((ETMX+ETMN)*MLR*S-
     ,       (S+MT**2-MT1**2)*MLR**2)/2.D0
        BK2=(ETMX-ETMN)*((ETMX+ETMN)*MLR-S-MT**2+MT1**2)/2.D0
        BK3=SNZJ*BI*AJ*MZJ+SNZI*BJ*AI*MZI
        BK=BI*BJ*BK1+AI*AJ*MZI*MZJ*SNZI*SNZJ*BK2+BK3*MT*(MLR**2-S)*
     ,     (ETMX-ETMN)/2.D0
        WID=BEZI*BEZJ*BK/(S-MZI**2)/(S-MZJ**2)
        SSLRT1=WID
        RETURN
        END
