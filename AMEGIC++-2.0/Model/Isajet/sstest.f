CDECK  ID>, SSTEST.
      SUBROUTINE SSTEST(IALLOW)
C
C          Test MSSM parameters against existing bounds on SUSY from
C          LEP and SLC:
C          IALLOW = 1    Z1 is not LSP
C          IALLOW = 2    Gamma(Z -> Z1SS Z1SS) < GAMINV
C          IALLOW = 4    Z -> charginos allowed
C          IALLOW = 8    BF(Z -> Z1SS Z2SS)>10^5
C          IALLOW = 16   Z -> squarks, sleptons
C          IALLOW = 32   BR(Z -> Z* HL0) < B(Z -> Z* H(M=MHSM))
C          IALLOW = 64   BR(Z -> HL0 HA0) > 0
C          IALLOW = 128  M(H+) > M(Z)/2
C          where GAMINV is the present bound on the invisible width,
C          and MHSM is the lower bound on the standard Higgs mass.
C
C          Bounds on the other modes are only approximate, but the
C          error in the allowed region of masses must be tiny.
C          Updated by H. Baer on 5/25/95
C
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
      PARAMETER (ISUPL=21,ISDNL=22,ISSTL=23,ISCHL=24,ISBT1=25,ISTP1=26)
      PARAMETER (ISNEL=31,ISEL=32,ISNML=33,ISMUL=34,ISNTL=35,ISTAU1=36)
      PARAMETER (ISUPR=41,ISDNR=42,ISSTR=43,ISCHR=44,ISBT2=45,ISTP2=46)
      PARAMETER (ISNER=51,ISER=52,ISNMR=53,ISMUR=54,ISNTR=55,ISTAU2=56)
      PARAMETER (ISGL=29)
      PARAMETER (ISZ1=30,ISZ2=40,ISZ3=50,ISZ4=60,ISW1=39,ISW2=49)
      PARAMETER (ISHL=82,ISHH=83,ISHA=84,ISHC=86)
      PARAMETER (ISGRAV=91)
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
C     GSS(28) = M_nr       GSS(29) = A_n
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
      COMMON /SUGMG/ MSS(32),GSS(29),MGUTSS,GGUTSS,AGUTSS,FTGUT,
     $FBGUT,FTAGUT,FNGUT
      REAL MSS,GSS,MGUTSS,GGUTSS,AGUTSS,FTGUT,FBGUT,FTAGUT,FNGUT
      SAVE /SUGMG/
      COMMON/XMSSM/GOMSSM,GOSUG,GOGMSB,GOAMSB,AL3UNI
     $,XGLSS,XMUSS,XHASS,XTBSS
     $,XQ1SS,XDRSS,XURSS,XL1SS,XERSS
     $,XQ2SS,XSRSS,XCRSS,XL2SS,XMRSS
     $,XQ3SS,XBRSS,XTRSS,XL3SS,XTARSS,XATSS,XABSS,XATASS
     $,XM1SS,XM2SS,XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU
     $,XLAMGM,XMESGM,XN5GM,XCMGV,XMGVTO
     $,XRSLGM,XDHDGM,XDHUGM,XDYGM,XN51GM,XN52GM,XN53GM
     $,XMN3NR,XMAJNR,XANSS,XNRSS,XSBCS
      SAVE /XMSSM/
      REAL XGLSS,XMUSS,XHASS,XTBSS
     $,XQ1SS,XDRSS,XURSS,XL1SS,XERSS
     $,XQ2SS,XSRSS,XCRSS,XL2SS,XMRSS
     $,XQ3SS,XBRSS,XTRSS,XL3SS,XTARSS,XATSS,XABSS,XATASS
     $,XM1SS,XM2SS
     $,XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU
     $,XLAMGM,XMESGM,XN5GM,XCMGV,XMGVTO
     $,XRSLGM,XDHDGM,XDHUGM,XDYGM,XN51GM,XN52GM,XN53GM
     $,XMN3NR,XMAJNR,XANSS,XNRSS,XSBCS
      LOGICAL GOMSSM,GOSUG,GOGMSB,GOAMSB,AL3UNI

C
      INTEGER IALLOW
      EXTERNAL SSZHX
      REAL MHSM,GAMINV,PI,SR2,G,GP,MZ,MZ1,MZ2,MZ3,MZ4,MW1,MW2,
     $TANB,BETA,COS2B,SIN2B,VS,V,VP,FT,MHL,ALPHA,SUSYCC,
     $GAMSS,W11,GZ1Z1,GAMSM,SSXINT,SSXLAM,COS2W,
     $W12,GZ1Z2,DGAMZ,BFZ,BFZ1Z2
C
C          Current bounds
      DATA MHSM/64./,GAMINV/.0043/,DGAMZ/.0115/,BFZ/1.E-5/
C
C          Initialize
C
      IALLOW=0
      PI=4.*ATAN(1.)
      SR2=SQRT(2.)
      G=SQRT(4*PI*ALFAEM/SN2THW)
      GP=G*SQRT(SN2THW/(1.-SN2THW))
      COS2W=SQRT(1.-SN2THW)
      MZ=AMZ
      MZ1=ABS(AMZ1SS)
      MZ2=ABS(AMZ2SS)
      MZ3=ABS(AMZ3SS)
      MZ4=ABS(AMZ4SS)
      MW1=ABS(AMW1SS)
      MW2=ABS(AMW2SS)
C
C          Check that Z1SS is LSP
C
      IF(MZ1.GT.MW1.OR.MZ1.GT.AMGLSS.OR.MZ1.GT.AMULSS
     $.OR.MZ1.GT.AMERSS.OR.MZ1.GT.AMELSS.OR.MZ1.GT.AMN1SS
     $.OR.MZ1.GT.AMB1SS.OR.MZ1.GT.AMT1SS.OR.MZ1.GT.AML1SS) THEN
        IALLOW=IALLOW+1
      ENDIF
C
C          Z -> Z1SS + Z1SS
C
      IF (AMZ.GT.2*MZ1) THEN
        W11=SQRT(G*G+GP*GP)
     $  *(ZMIXSS(1,1)*ZMIXSS(1,1)-ZMIXSS(2,1)*ZMIXSS(2,1))/4.
        GZ1Z1=SQRT(SSXLAM(MZ**2,MZ1**2,MZ1**2))/12./PI/MZ**3*W11**2
     $  *(2*MZ**2-MZ1**2-MZ1**2-(MZ1**2-MZ1**2)**2/MZ**2
     $  -6*MZ1*MZ1*SIGN(1.,AMZ1SS*AMZ1SS))
        IF(GZ1Z1.GT.GAMINV) THEN
          IALLOW=IALLOW+2
        ENDIF
      ENDIF
C
C          Check for other allowed visible modes modes
C
      IF(AMZ.GT.2*MW1) THEN
        IALLOW=IALLOW+4
      ENDIF
C
C          Check funny Z branching fractions
C
      BFZ1Z2=0.
      IF (AMZ.GT.MZ1+MZ2) THEN
        W12=SQRT(G*G+GP*GP)
     $  *(ZMIXSS(1,1)*ZMIXSS(1,2)-ZMIXSS(2,1)*ZMIXSS(2,2))/4.
        GZ1Z2=SQRT(SSXLAM(MZ**2,MZ1**2,MZ2**2))/6./PI/MZ**3*W12**2
     $  *(2*MZ**2-MZ1**2-MZ2**2-(MZ1**2-MZ2**2)**2/MZ**2
     $  -6*MZ1*MZ2*SIGN(1.,AMZ1SS*AMZ2SS))
        BFZ1Z2=GZ1Z2/GAMZ
      END IF
      IF(BFZ1Z2.GT.BFZ) THEN
        IALLOW=IALLOW+8
      ENDIF
C
      IF(AMZ.GT.2*AMULSS.OR.AMZ.GT.2*AMELSS.OR.AMZ.GT.2*AMERSS
     $.OR.AMZ.GT.2*AMN1SS.OR.AMZ.GT.2*AMB1SS.OR.AMZ.GT.2*AMT1SS)THEN
        IALLOW=IALLOW+16
      ENDIF
C
C          Z -> Higgs modes
C
      TMP(1)=MHSM
      GAMSM=SSXINT(2*MHSM/MZ,SSZHX,(1.+MHSM**2/MZ**2))
C          Z -> hl0 x
      IF(AMZ.GT.AMHL) THEN
        TANB=1./RV2V1
        BETA=ATAN(TANB)
        COS2B=COS(2*BETA)
        SIN2B=SIN(2*BETA)
        VS=2*AMW**2/G**2/(1.+RV2V1**2)
        V=SQRT(VS)
        VP=RV2V1*V
        FT=G*AMTP/SR2/AMW/V*SQRT(V**2+VP**2)
        MHL=AMHL
        ALPHA=ALFAH
        SUSYCC=SIN(ALPHA+BETA)
        TMP(1)=MHL
        GAMSS=SSXINT(2*MHL/AMZ,SSZHX,(1.+MHL**2/AMZ**2))*SUSYCC**2
        IF(GAMSS.GE.GAMSM) IALLOW=IALLOW+32
      ENDIF
C          Z -> hl0 ha0
      IF (AMZ.GT.(AMHL+AMHA)) THEN
        IALLOW=IALLOW+64
      ENDIF
C          Z -> H+ H-
      IF(AMZ.GT.2*AMHC) THEN
        IALLOW=IALLOW+128
      ENDIF
C
      RETURN
      END
