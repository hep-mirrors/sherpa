CDECK  ID>, SSGZT.
        REAL FUNCTION SSGZT(E)
C-----------------------------------------------------------------------
C          SSGLBF: glss -> ziss + tp + tb
C          Baer's TOPINT
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
C          Temporary parameters for functions
      COMMON/SSTMP/TMP(10),ITMP(10)
      REAL TMP
      INTEGER ITMP
      SAVE /SSTMP/
C
      REAL E
      DOUBLE PRECISION MG,MT,MS,MZ,SSDLAM,P,PSI,XLOG,PHI,SN,C1,C2
      DOUBLE PRECISION XT,MUT,MUZ,XMIN,XMAX,EMIN,EMAX
C          Convert to double precision
      MG=AMGLSS
      C1=TMP(1)
      C2=TMP(2)
      MS=TMP(3)
      MT=TMP(4)
      MZ=TMP(5)
      SN=TMP(6)
C
      XT=2*E/MG
      MUT=(MT/MG)**2
      MUZ=(MZ/MG)**2
      XMIN=((2.D0-XT)*(1.D0+2*MUT-MUZ-XT)-DSQRT((XT**2-4*MUT)*
     $SSDLAM((1.D0+MUT-XT),MUT,MUZ)))/2.D0/(1.D0-XT+MUT)
      XMAX=((2.D0-XT)*(1.D0+2*MUT-MUZ-XT)+DSQRT((XT**2-4*MUT)*
     $SSDLAM((1.D0+MUT-XT),MUT,MUZ)))/2.D0/(1.D0-XT+MUT)
      EMIN=XMIN*MG/2.D0
      EMAX=XMAX*MG/2.D0
      P=SQRT(E**2-MT**2)
      PSI=P*E*(MG**2-MZ**2-2*MG*E)*
     $DSQRT(SSDLAM((MG**2+MT**2-2*MG*E),MZ**2,MT**2))/MG/
     $(MG**2+MT**2-2*MG*E)/(MG**2+MT**2-2*MG*E-MS**2)**2
      XLOG=DLOG((MG**2+MT**2-2*MG*EMAX-MS**2)/
     $(MG**2+MT**2-2*MG*EMIN-MS**2))
      PHI=MZ*(-(EMAX-EMIN)-(2*E*MG+MZ**2-MT**2-MS**2)*
     $XLOG/2.D0/MG)/
     $2.D0/MG/(MG**2+MT**2-MS**2-2*MG*E)
      SSGZT=C1*PSI+SN*C2*PHI
      RETURN
      END
