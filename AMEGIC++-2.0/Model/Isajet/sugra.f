CDECK  ID>, SUGRA.
C--------------------------------------------------------------------
      SUBROUTINE SUGRA(M0,MHF,A0,TANB,SGNMU,MT,IMODEL)
C--------------------------------------------------------------------
C
C     Calculate supergravity spectra for ISAJET using as inputs
C     M0    = M_0       = common scalar mass at GUT scale
C     MHF   = M_(1/2)   = common gaugino mass at GUT scale
C     A0    = A_0       = trilinear soft breaking parameter at GUT scale
C     TANB  = tan(beta) = ratio of vacuum expectation values v_1/v_2
C     SGNMU = sgn(mu)   = +-1 = sign of Higgsino mass term
C     MT    = M_t       = mass of t quark
C     M0    = Lambda    = ratio of vevs <F>/<S>
C     MHF   = M_Mes     = messenger scale
C     A0    = n_5       = number of messenger fields
C     IMODEL            = 1 for SUGRA model
C                       = 2 for GMSB model
C                       = 7 for AMSB model
C
C     Uses Runge-Kutta method to integrate RGE's from M_Z to M_GUT
C     and back, putting in correct thresholds. For the first iteration
C     only the first 6 couplings are included and a common threshold
C     is used.
C
C     See /SUGMG/ for definitions of couplings and masses.
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
      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ
      INTEGER NOGOOD,IAL3UN,ITACHY
      SAVE /SUGPAS/
C     XNUSUG contains non-universal GUT scale soft terms for SUGRA:
C     XNUSUG(1)=M1 XNUSUG(2)=M2 XNUSUG(3)=M3
C     XNUSUG(4)=A_tau XNUSUG(5)=A_b XNUSUG(6)=A_t
C     XNUSUG(7)=m_Hd XNUSUG(8)=m_Hu XNUSUG(9)=m_eR XNUSUG(10)=m_eL
C     XNUSUG(11)=m_dR XNUSUG(12)=m_uR XNUSUG(13)=m_uL XNUSUG(14)=m_lR
C     XNUSUG(15)=m_lL XNUSUG(16)=m_bR XNUSUG(17)=m_tR XNUSUG(18)=m_tL
C
      COMMON /SUGNU/ XNUSUG(18)
      REAL XNUSUG
      SAVE /SUGNU/
      COMMON/SSINF/XLAM
      DOUBLE PRECISION XLAM
      REAL GY(7),W1(21),G(29),W2(87)
      REAL G0(29)
      COMPLEX*16 SSB0,SSB1
      DOUBLE PRECISION DDILOG,XLM
      INTEGER IG(29)
      EXTERNAL SURG06,SURG26
      REAL M0,MHF,A0,TANB,SGNMU,MT,XLAMGM,XMESGM,XN5GM
      INTEGER NSTEP
      REAL M2,SUALFE,SUALFS,Q,T,A1I,AGUT,A3I,A2I,MTMT,ASMT,DT,
     $TGUT,TZ,GGUT,SIG2,SIG1,MH1S,MH2S,AGUTI,
     $MUS,MBMZ,MB,MTAU,MZ,MW,SR2,PI,ALEM,MTAMZ,
     $MTAMB,MTAMTA,MBMB,ASMB,BETA,COTB,SINB,COS2B,COSB,XC,
     $MSN,MG,MT1,MT2,MB1,MB2,MW1,MW2,AMU,BTHAT,BBHAT,BLHAT,AM2
      INTEGER II,I,J,IMODEL
      REAL G0SAVE(26),DELG0,DELLIM,THRF,THRG,DY,QOLD
      INTEGER MXITER,NSTEP0
      COMPLEX*16 ZZZ
C LESTER Added the following line
      REAL*8 DREAL
      REAL*8 REAL8
C
      DATA MZ/91.187/,MTAU/1.777/,MB/4.9/,ALEM/.0078186/
C          This choice is a compromise between precision and speed:
      DATA MXITER/20/,NSTEP0/200/,DELLIM/2.E-2/
C
C          Define REAL(COMPLEX*16) for g77. This might need to be
C          changed for 64-bit machines?
C
      REAL8(ZZZ)=DREAL(ZZZ)
C
C          Save input parameters
C
      XSUGIN(1)=M0
      XSUGIN(2)=MHF
      XSUGIN(3)=A0
      XSUGIN(4)=TANB
      XSUGIN(5)=SGNMU
      XSUGIN(6)=MT
      XLAMGM=M0
      XMESGM=MHF
      XN5GM=A0
      XGMIN(1)=XLAMGM
      XGMIN(2)=XMESGM
      XGMIN(3)=XN5GM
      XGMIN(4)=TANB
      XGMIN(5)=SGNMU
      XGMIN(6)=MT
      IF (XGMIN(12).EQ.0.) XGMIN(12)=XN5GM
      IF (XGMIN(13).EQ.0.) XGMIN(13)=XN5GM
      IF (XGMIN(14).EQ.0.) XGMIN(14)=XN5GM
C
C          Compute gauge mediated threshold functions
C
      IF (IMODEL.EQ.2) THEN
        XLM=XLAMGM/XMESGM
        THRF=((1.D0+XLM)*(LOG(1.D0+XLM)-2*DDILOG(XLM/(1.D0+XLM))+
     ,        .5*DDILOG(2*XLM/(1.D0+XLM)))+
     ,       (1.D0-XLM)*(LOG(1.D0-XLM)-2*DDILOG(-XLM/(1.D0-XLM))+
     ,        .5*DDILOG(-2*XLM/(1.D0-XLM))))/XLM**2
        THRG=((1.D0+XLM)*LOG(1.D0+XLM)+(1.D0-XLM)*LOG(1.D0-XLM))/XLM**2
      END IF
C
C          Initialize standard model parameters in /SSSM/:
C
      AMUP=0.0056
      AMDN=0.0099
      AMST=0.199
      AMCH=1.35
      AMBT=5.0
      AMTP=MT
      AMT=MT
      AME=0.511E-3
      AMMU=0.105
      AMTAU=1.777
      AMZ=91.17
      GAMW=2.12
      GAMZ=2.487
      ALFAEM=1./128.
      SN2THW=0.232
      ALFA2=ALFAEM/SN2THW
      ALQCD4=0.177
      ALFA3=0.118
C
      NOGOOD=0
      ITACHY=0
      PI=4.*ATAN(1.)
      SR2=SQRT(2.)
      XW=.2324-1.03E-7*(MT**2-138.**2)
      MW=MZ*SQRT(1.-XW)
      AMW=MW
      A1MZ=5*ALEM/3./(1.-XW)
      A2MZ=ALEM/XW
      G2=SQRT(4*PI*A2MZ)
      GP=SQRT(3./5.*A1MZ*4.*PI)
      XTANB=TANB
      COTB=1./TANB
      BETA=ATAN(TANB)
      SINB=SIN(BETA)
      COSB=COS(BETA)
      SIN2B=SIN(2*BETA)
      COS2B=COS(2*BETA)
      IF (IMODEL.EQ.1) THEN
        MSUSY=SQRT(M0**2+4*MHF**2)
      ELSE IF (IMODEL.EQ.2) THEN
        MSUSY=XLAMGM/100.
      ELSE IF (IMODEL.EQ.7) THEN
        MSUSY=SQRT(M0**2+(.01*MHF)**2)
      END IF
C     USE PIERCE PRESCRIPTION FOR MAGNITUDE OF VEV
C      VEV=SR2*(248.6+0.9*LOG(MSUSY/AMZ)
C      V=SQRT(VEV**2/(1.+COTB))
C     PREVIOUS PRESCRIPTION
      V=SQRT(2*MW**2/G2**2/(1.+COTB**2))
      VP=V/TANB
      VEV=SQRT(V**2+VP**2)
C
C          Compute m(tau), m(b) at z scale using qcd, qed
C
      MTAMTA=MTAU*(1.-SUALFE(MTAU**2)/PI)
      MTAMB=MTAMTA*(SUALFE(MB**2)/SUALFE(MTAU**2))**(-27./76.)
      MTAMZ=MTAMB*(SUALFE(MZ**2)/SUALFE(MB**2))**(-27./80.)
      FTAMZ=MTAMZ/COSB/VEV
      ASMB=SUALFS(MB**2,.36,MT,3)
      MBMB=MB*(1.-4*ASMB/3./PI)
      ASMZ=SUALFS(MZ**2,.36,MT,3)
      MBMZ=MBMB*(ASMZ/ASMB)**(12./23.)*
     $      (SUALFE(MZ**2)/SUALFE(MB**2))**(-3./80.)
      FBMZ=MBMZ/COSB/VEV
      ASMT=SUALFS(MT**2,.36,MT,3)
      MTMT=MT/(1.+4*ASMT/3./PI+(16.11-1.04*(5.-6.63/MT))*(ASMT/PI)**2)
      FTMT=MTMT/SINB/VEV
      FNMZ=SQRT(XNRIN(2)*XNRIN(1)/(SINB*VEV)**2)
      AMNRMJ=XNRIN(2)
C
C          Run the 3 gauge and 3 Yukawa's up to find M_GUT ,A_GUT and
C          Yukawa_GUT
C
C
      NSTEP=NSTEP0
      GY(1)=SQRT(4*PI*A1MZ)
      GY(2)=SQRT(4*PI*A2MZ)
      GY(3)=SQRT(4*PI*ALFA3)
      GY(4)=FTAMZ
      GY(5)=FBMZ
      GY(6)=0.
      GY(7)=0.
      IF (IMODEL.EQ.1.OR.IMODEL.EQ.7) THEN
        IF (XSUGIN(7).EQ.0.) THEN
          MGUT=1.E19
        ELSE
          MGUT=XSUGIN(7)
        END IF
      ELSE IF (IMODEL.EQ.2) THEN
        MGUT=XMESGM
      END IF
      TZ=LOG(MZ/MGUT)
      TGUT=0.
      DT=(TGUT-TZ)/FLOAT(NSTEP)
      DO 200 II=1,NSTEP
        T=TZ+(TGUT-TZ)*FLOAT(II-1)/FLOAT(NSTEP)
        Q=MGUT*EXP(T)
        IF (Q.GT.MT.AND.GY(6).EQ.0.) GY(6)=FTMT
        IF (Q.GT.XNRIN(2).AND.GY(7).EQ.0.) GY(7)=FNMZ
        CALL RKSTP(7,DT,T,GY,SURG06,W1)
        A1I=4*PI/GY(1)**2
        A2I=4*PI/GY(2)**2
        A3I=4*PI/GY(3)**2
        IF (GY(5).GT.10..OR.GY(6).GT.10..OR.GY(7).GT.10.) THEN
          NOGOOD=4
          GO TO 100
        END IF
        IF (A1I.LT.A2I.AND.XSUGIN(7).EQ.0.) GO TO 10
200   CONTINUE
      IF (MGUT.EQ.1.E19) THEN
      WRITE(LOUT,*) 'SUGRA: NO UNIFICATION FOUND'
      GO TO 100
      END IF
10    IF (XSUGIN(7).EQ.0.) THEN
        MGUT=Q
      ELSE
        MGUT=XSUGIN(7)
      END IF
      AGUT=(GY(1)**2/4./PI+GY(2)**2/4./PI)/2.
      GGUT=SQRT(4*PI*AGUT)
      AGUTI=1./AGUT
      FTAGUT=GY(4)
      FBGUT=GY(5)
      FTGUT=GY(6)
      IF (XNRIN(1).EQ.0..AND.XNRIN(2).LT.1.E19) THEN
C       UNIFY FN-FT
        FNGUT=GY(6)
      ELSE
        FNGUT=GY(7)
      END IF
C
C          Define parameters at GUT scale
C
      DO 210 J=1,3
        IF (IMODEL.EQ.1) THEN
          G(J)=GY(J)
          G(J+6)=MHF
          G(J+9)=A0
        ELSE IF (IMODEL.EQ.2) THEN
          G(J)=GY(J)
          G(J+6)=XGMIN(11+J)*XGMIN(8)*THRG*(GY(J)/4./PI)**2*XLAMGM
          G(J+9)=0.
        END IF
210   CONTINUE
C     OVERWRITE ALFA_3 UNIFICATION TO GET ALFA_3(MZ) RIGHT
      IF (IMODEL.EQ.1.AND.IAL3UN.NE.0) G(3)=GGUT
      G(4)=FTAGUT
      G(5)=FBGUT
      G(6)=FTGUT
C     IF NR MAJORANA MASS EXISTS, SET EXTRA NR RGE PARAMETERS
      IF (XNRIN(2).LT.1.E19) THEN
        G(27)=FNGUT
        G(28)=XNRIN(4)**2
        G(29)=XNRIN(3)
      ELSE
        G(27)=0.
        G(28)=0.
        G(29)=0.
      END IF
      IF (IMODEL.EQ.1) THEN
        DO 220 J=13,24
          G(J)=M0**2
220     CONTINUE
C       Set possible non-universal boundary conditions
      DO 230 J=1,6
        IF (XNUSUG(J).LT.1.E19) THEN
          G(J+6)=XNUSUG(J)
        END IF
230   CONTINUE
      DO 231 J=7,18
        IF (XNUSUG(J).LT.1.E19) THEN
          G(J+6)=XNUSUG(J)**2
        END IF
231   CONTINUE
      ELSE IF (IMODEL.EQ.2) THEN
       XC=2*THRF*XLAMGM**2
       DY=SQRT(3./5.)*GY(1)*XGMIN(11)
       G(13)=XC*(.75*XGMIN(13)*(GY(2)/4./PI)**4+.6*.25*
     ,XGMIN(12)*(GY(1)/4./PI)**4)+XGMIN(9)-DY
       G(14)=XC*(.75*XGMIN(13)*(GY(2)/4./PI)**4+.6*.25*
     ,XGMIN(12)*(GY(1)/4./PI)**4)+XGMIN(10)+DY
       G(15)=XC*(.6*XGMIN(12)*(GY(1)/4./PI)**4)+2*DY
       G(16)=XC*(.75*XGMIN(13)*(GY(2)/4./PI)**4+.6*.25*
     ,XGMIN(12)*(GY(1)/4./PI)**4)-DY
       G(17)=XC*(4*XGMIN(14)*(GY(3)/4./PI)**4/3.+.6*XGMIN(12)*
     ,(GY(1)/4./PI)**4/9.)+2*DY/3.
       G(18)=XC*(4*XGMIN(14)*(GY(3)/4./PI)**4/3.+.6*4*XGMIN(12)*
     ,(GY(1)/4./PI)**4/9.)-4*DY/3.
       G(19)=XC*(4*XGMIN(14)*(GY(3)/4./PI)**4/3.+.75*XGMIN(13)*
     ,(GY(2)/4./PI)**4+.6*XGMIN(12)*(GY(1)/4./PI)**4/36.)+DY/3.
       G(20)=G(15)
       G(21)=G(16)
       G(22)=G(17)
       G(23)=G(18)
       G(24)=G(19)
      ELSE IF (IMODEL.EQ.7) THEN
       G(1)=GY(1)
       G(2)=GY(2)
       G(3)=GY(3)
       BLHAT=G(4)*(-9*G(1)**2/5.-3*G(2)**2+3*G(5)**2+4*G(4)**2)
       BBHAT=G(5)*(-7*G(1)**2/15.-3*G(2)**2-16*G(3)**2/3.+
     ,             G(6)**2+6*G(5)**2+G(4)**2)
       BTHAT=G(6)*(-13*G(1)**2/15.-3*G(2)**2-16*G(3)**2/3.+
     ,             6*G(6)**2+G(5)**2)
       G(7)=-33*MHF*G(1)**2/5./16./PI**2
       G(8)=-MHF*G(2)**2/16./PI**2
       G(9)=3*MHF*G(3)**2/16./PI**2
       G(10)=BLHAT*MHF/G(4)/16./PI**2
       G(11)=BBHAT*MHF/G(5)/16./PI**2
       G(12)=BTHAT*MHF/G(6)/16./PI**2
       G(13)=(-99*G(1)**4/50.-3*G(2)**4/2.+3*G(5)*BBHAT+G(4)*BLHAT)*
     ,        MHF**2/(16*PI**2)**2
       G(14)=(-99*G(1)**4/50.-3*G(2)**4/2.+3*G(6)*BTHAT)*
     ,        MHF**2/(16*PI**2)**2
       G(15)=(-198*G(1)**4/25.)*MHF**2/(16*PI**2)**2
       G(16)=(-99*G(1)**4/50.-3*G(2)**4/2.)*MHF**2/(16*PI**2)**2
       G(17)=(-22*G(1)**4/25.+8*G(3)**4)*MHF**2/(16*PI**2)**2
       G(18)=(-88*G(1)**4/25.+8*G(3)**4)*MHF**2/(16*PI**2)**2
       G(19)=(-11*G(1)**4/50.-3*G(2)**4/2.+8*G(3)**4)*
     ,        MHF**2/(16*PI**2)**2
       G(20)=(-198*G(1)**4/25.+2*G(4)*BLHAT)*MHF**2/(16*PI**2)**2
       G(21)=(-99*G(1)**4/50.-3*G(2)**4/2.+G(4)*BLHAT)*
     ,        MHF**2/(16*PI**2)**2
       G(22)=(-22*G(1)**4/25.+8*G(3)**4+2*G(5)*BBHAT)*
     ,MHF**2/(16*PI**2)**2
       G(23)=(-88*G(1)**4/25.+8*G(3)**4+2*G(6)*BTHAT)*
     ,MHF**2/(16*PI**2)**2
       G(24)=(-11*G(1)**4/50.-3*G(2)**4/2.+8*G(3)**4+G(5)*BBHAT+
     ,        G(6)*BTHAT)*MHF**2/(16*PI**2)**2
       DO 234 I=13,24
234      G(I)=G(I)+M0**2
      END IF
      G(25)=0.
      G(26)=0.
      DO 235 I=1,29
        IG(I)=0
235   CONTINUE
C          Check for tachyonic sleptons at GUT scale
      IF (G(15).LT.0..OR.G(16).LT.0.) THEN
        ITACHY=1
      END IF
C
C          Initialize thresholds
C
      MSS(1)=MSUSY
      MSS(2)=MSUSY
      MSS(17)=MSUSY
      MSS(27)=MSUSY
      MSS(31)=MSUSY
      MU=MSUSY
C
C          Evolve parameters from mgut to mz
C
      TZ=LOG(MZ/MGUT)
      TGUT=0.
      DT=(TZ-TGUT)/FLOAT(NSTEP)
C          Freeze Higgs parameters at HIGFRZ = Drees' value
C          AMTLSS, AMTRSS initialized to 0 for later use in HIGFRZ
      IF (IMODEL.EQ.1) THEN
        HIGFRZ=SQRT(M0**2+3*MHF**2)
      ELSE IF (IMODEL.EQ.2) THEN
        HIGFRZ=MSUSY
      ELSE IF (IMODEL.EQ.7) THEN
        HIGFRZ=SQRT(M0**2+(.01*MHF)**2)
      END IF
      AMTLSS=0
      AMTRSS=0
      DO 240 II=1,NSTEP+2
        T=TGUT+(TZ-TGUT)*FLOAT(II-1)/FLOAT(NSTEP)
        QOLD=Q
        Q=MGUT*EXP(T)
        CALL RKSTP(29,DT,T,G,SURG26,W2)
        IF (Q.LT.AMNRMJ.AND.QOLD.GE.AMNRMJ.AND.FNMZ.EQ.0.) THEN
          FNMZ=G(27)
        END IF
        IF (Q.LT.AMNRMJ) THEN
          G(27)=0.
          G(28)=0.
          G(29)=0.
        END IF
        CALL SUGFRZ(Q,G,G0,IG)
        IF (NOGOOD.NE.0) GO TO 100
        IF (Q.LT.MZ) GO TO 20
240   CONTINUE
20    CONTINUE
      ASMZ=G0(3)**2/4./PI
C          Electroweak breaking constraints; tree level
      MUS=(G0(13)-G0(14)*TANB**2)/(TANB**2-1.)-MZ**2/2.
      IF (MUS.LT.0.) THEN
        NOGOOD=2
        GO TO 100
      END IF
      MU=SQRT(MUS)*SIGN(1.,SGNMU)
      B=(G0(13)+G0(14)+2*MUS)*SIN2B/MU/2.
C          Compute tree level masses
      CALL SUGMAS(G0,0,IMODEL)
      IF (NOGOOD.NE.0) GO TO 100
C          Compute effective potential corrections
      CALL SUGEFF(G0,SIG1,SIG2)
      MH1S=G0(13)+SIG1
      MH2S=G0(14)+SIG2
      MUS=(MH1S-MH2S*TANB**2)/(TANB**2-1.)-MZ**2/2.
      IF (MUS.LT.0.) THEN
        NOGOOD=2
        GO TO 100
      END IF
      MU=SQRT(MUS)*SIGN(1.,SGNMU)
      B=(MH1S+MH2S+2*MUS)*SIN2B/MU/2.
C
C     Recompute weak scale Yukawa couplings including SUSY loops
C     Follow formulae of Pierce et al. NPB491, 3 (1997)
C
      M2=G0(8)
      AM2=ABS(M2)
      MSN=MSS(16)
      MG=MSS(1)
      MT1=MSS(12)
      MT2=MSS(13)
      MB1=MSS(10)
      MB2=MSS(11)
      MW1=ABS(MSS(27))
      MW2=ABS(MSS(28))
      AMU=ABS(MU)
      XLAM=LOG(MT**2)
C          Be careful in using our convention vs Pierce et al.
      IF (ABS(COS(THETAT)).LT..707107) THEN
        MTMT=MT/(1.+4*ASMT/3./PI+(16.11-1.04*(5.-6.63/MT))*(ASMT/PI)**2
     $  -ASMT/3./PI*(REAL8(SSB1(MT**2,MG,MT1))+
     $  REAL8(SSB1(MT**2,MG,MT2))-SIN(2*THETAT)*MG/MT*
     $  (REAL8(SSB0(MT**2,MG,MT1))-REAL8(SSB0(MT**2,MG,MT2)))))
      ELSE
        MTMT=MT/(1.+4*ASMT/3./PI+(16.11-1.04*(5.-6.63/MT))*(ASMT/PI)**2
     $  -ASMT/3./PI*(REAL8(SSB1(MT**2,MG,MT2))+
     $  REAL8(SSB1(MT**2,MG,MT1))+SIN(2*THETAT)*MG/MT*
     $  (REAL8(SSB0(MT**2,MG,MT2))-REAL8(SSB0(MT**2,MG,MT1)))))
      END IF
      FTMT=MTMT/SINB/VEV
      XLAM=LOG(MZ**2)
      IF (ABS(COS(THETAB)).LT..707107) THEN
        MBMZ=MBMZ*(1.+ASMZ/3./PI*(REAL8(SSB1(MZ**2,MG,MB1))+
     $  REAL8(SSB1(MZ**2,MG,MB2))-SIN(2*THETAB)*MG/MB*
     $  (REAL8(SSB0(MZ**2,MG,MB1))-REAL8(SSB0(MZ**2,MG,MB2))))
     $  -FTMT**2*MU*(-AAT*TANB+MU)/16./PI**2/(MT1**2-MT2**2)*
     $  (REAL8(SSB0(MZ**2,AMU,MT1))-REAL8(SSB0(MZ**2,AMU,MT2)))+
     $  G2**2*MU*M2*TANB/16./PI**2/(AMU**2-M2**2)*
     $  (SIN(THETAT)**2*(REAL8(SSB0(MZ**2,AM2,MT1))-
     $  REAL8(SSB0(MZ**2,AMU,MT1)))+
     $  COS(THETAT)**2*(REAL8(SSB0(MZ**2,AM2,MT2))-
     $  REAL8(SSB0(MZ**2,AMU,MT2)))))
      ELSE
        MBMZ=MBMZ*(1.+ASMZ/3./PI*(REAL8(SSB1(MZ**2,MG,MB2))+
     $  REAL8(SSB1(MZ**2,MG,MB1))+SIN(2*THETAB)*MG/MB*
     $  (REAL8(SSB0(MZ**2,MG,MB2))-REAL8(SSB0(MZ**2,MG,MB1))))
     $  -FTMT**2*MU*(-AAT*TANB+MU)/16./PI**2/(MT2**2-MT1**2)*
     $  (REAL8(SSB0(MZ**2,AMU,MT2))-REAL8(SSB0(MZ**2,AMU,MT1)))+
     $  G2**2*MU*M2*TANB/16./PI**2/(AMU**2-M2**2)*
     $  (COS(THETAT)**2*(REAL8(SSB0(MZ**2,AM2,MT2))-
     $  REAL8(SSB0(MZ**2,AMU,MT2)))+
     $  SIN(THETAT)**2*(REAL8(SSB0(MZ**2,AM2,MT1))-
     $  REAL8(SSB0(MZ**2,AMU,MT1)))))
      END IF
      FBMZ=MBMZ/COSB/VEV
      MTAMZ=MTAMZ*(1.+G2**2*MU*M2*TANB/16./PI**2/(MUS-M2**2)*
     $(REAL8(SSB0(MZ**2,AM2,MSN))-REAL8(SSB0(MZ**2,AMU,MSN))))
      FTAMZ=MTAMZ/COSB/VEV
C
C          Iterate entire process, increasing NSTEP each time
C          This time, freeze out parameters at sqrt(t_l t_r)
C
      HIGFRZ=MAX(AMZ,(G0(23)*G0(24))**0.25)
      DO 300 I=1,MXITER
        DO 310 J=1,26
310     G0SAVE(J)=G0(J)
        NSTEP=1.2*NSTEP
        CALL SUGRGE(M0,MHF,A0,TANB,SGNMU,MT,G,G0,IG,W2,NSTEP,IMODEL)
        IF(NOGOOD.NE.0) GO TO 100
        DELG0=0.
        DO 320 J=1,24
320     DELG0=MAX(DELG0,ABS((G0(J)-G0SAVE(J))/G0(J)))
        IF(DELG0.LT.DELLIM) GO TO 400
300   CONTINUE
      WRITE(LOUT,1000) MXITER
1000  FORMAT(/' SUGRA WARNING: NO CONVERGENCE IN',I4,' ITERATIONS')
C
C          Save results
C
400   DO 410 I=1,26
        GSS(I)=G0(I)
410   CONTINUE
      MGUTSS=MGUT
      AGUTSS=AGUT
      GGUTSS=GGUT
C
C          Fill XISAIN common block
C
      XISAIN(1)=MSS(1)
      XISAIN(2)=MU
      XISAIN(3)=MSS(31)
      XISAIN(4)=TANB
      XISAIN(5)=SQRT(G0(19))
      XISAIN(6)=SQRT(G0(17))
      XISAIN(7)=SQRT(G0(18))
      XISAIN(8)=SQRT(G0(16))
      XISAIN(9)=SQRT(G0(15))
      XISAIN(10)=XISAIN(5)
      XISAIN(11)=XISAIN(6)
      XISAIN(12)=XISAIN(7)
      XISAIN(13)=XISAIN(8)
      XISAIN(14)=XISAIN(9)
      XISAIN(15)=SQRT(G0(24))
      XISAIN(16)=SQRT(G0(22))
      XISAIN(17)=SQRT(G0(23))
      XISAIN(18)=SQRT(G0(21))
      XISAIN(19)=SQRT(G0(20))
      XISAIN(20)=G0(12)
      XISAIN(21)=G0(11)
      XISAIN(22)=G0(10)
      XISAIN(23)=G0(7)
      XISAIN(24)=G0(8)
      M2=G0(8)
100   RETURN
      END











