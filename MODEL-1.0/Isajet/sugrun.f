CDECK  ID>, SSID.   
      CHARACTER*5 FUNCTION SSID(ID)
C-----------------------------------------------------------------------
C
C     Return character name for ID, assuming the default IDENT codes
C     are used in /SSTYPE/.
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
      COMMON/SSLUN/LOUT
      INTEGER LOUT
      SAVE /SSLUN/
      CHARACTER*5 LABEL(-120:120)
      SAVE LABEL
      INTEGER ID,J
C
      DATA LABEL(0)/'     '/
C
      DATA (LABEL(J),J=1,10)
     $/'UP   ','DN   ','ST   ','CH   ','BT   ','TP   '
     $,'ERROR','ERROR','GL   ','GM   '/
      DATA (LABEL(J),J=-1,-10,-1)
     $/'UB   ','DB   ','SB   ','CB   ','BB   ','TB   '
     $,'ERROR','ERROR','ERROR','ERROR'/
C
      DATA (LABEL(J),J=11,20)
     $/'NUE  ','E-   ','NUM  ','MU-  ','NUT  ','TAU- '
     $,'ERROR','ERROR','ERROR','ERROR'/
      DATA (LABEL(J),J=-11,-20,-1)
     $/'ANUE ','E+   ','ANUM ','MU+  ','ANUT ','TAU+ '
     $,'ERROR','ERROR','ERROR','ERROR'/
C
      DATA (LABEL(J),J=21,30)
     $/'UPL  ','DNL  ','STL  ','CHL  ','BT1  ','TP1  '
     $,'ERROR','ERROR','GLSS ','Z1SS '/
      DATA (LABEL(J),J=-21,-30,-1)
     $/'UBL  ','DBL  ','SBL  ','CBL  ','BB1  ','TB1  '
     $,'ERROR','ERROR','ERROR','ERROR'/
C
      DATA (LABEL(J),J=31,40)
     $/'NUEL ','EL-  ','NUML ','MUL- ','NUTL ','TAU1-'
     $,'ERROR','ERROR','W1SS+','Z2SS '/
      DATA (LABEL(J),J=-31,-40,-1)
     $/'ANUEL','EL+  ','ANUML','MUL+ ','ANUTL','TAU1+'
     $,'ERROR','ERROR','W1SS-','ERROR'/
C
      DATA (LABEL(J),J=41,50)
     $/'UPR  ','DNR  ','STR  ','CHR  ','BT2  ','TP2  '
     $,'ERROR','ERROR','W2SS+','Z3SS '/
      DATA (LABEL(J),J=-41,-50,-1)
     $/'UBR  ','DBR  ','SBR  ','CBR  ','BB2  ','TB2  '
     $,'ERROR','ERROR','W2SS-','ERROR'/
C
      DATA (LABEL(J),J=51,60)
     $/'NUER ','ER-  ','NUMR ','MUR- ','NUTR ','TAU2-'
     $,'ERROR','ERROR','ERROR','Z4SS '/
      DATA (LABEL(J),J=-51,-60,-1)
     $/'ANUEL','ER+  ','ANUMR','MUR+ ','ANUTR','TAU2+'
     $,'ERROR','ERROR','ERROR','ERROR'/
C
      DATA (LABEL(J),J=82,86)
     $/'HL0  ','HH0  ','HA0  ','ERROR','H+   '/
      DATA LABEL(-86)/'H-   '/
C
      DATA LABEL(80)/'W+   '/,LABEL(-80)/'W-   '/,LABEL(90)/'Z0   '/
      DATA LABEL(91)/'GVSS '/
      DATA LABEL(120)/'PI+  '/,LABEL(-120)/'PI-  '/
C
      IF(IABS(ID).GT.120) THEN
        WRITE(LOUT,*) 'SSID: ID = ',ID
        STOP99
      ENDIF
      SSID=LABEL(ID)
      RETURN
      END
CDECK  ID>, SSPRT.  
      SUBROUTINE SSPRT(ID)
C-----------------------------------------------------------------------
C
C     Print decay modes for ID. Note these need not be contiguous,
C     so the loop is over all modes in /SSMODE/.
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
C
      INTEGER ID,I,K,NOUT
      CHARACTER*5 SSID,LBLIN,LBLOUT(3)
C
      NOUT=0
      DO 100 I=1,NSSMOD
        IF(ISSMOD(I).NE.ID) GO TO 100
        NOUT=NOUT+1
        LBLIN=SSID(ISSMOD(I))
        DO 110 K=1,3
110     LBLOUT(K)=SSID(JSSMOD(K,I))
        WRITE(LOUT,1000) LBLIN,(LBLOUT(K),K=1,3),GSSMOD(I),BSSMOD(I)
1000    FORMAT(1X,A5,'  -->  ',3(A5,2X),2E15.5)
100   CONTINUE
C
      IF(NOUT.GT.0) WRITE(LOUT,*) ' '
C
      RETURN
      END
CDECK  ID>, SUGPRT. 
C--------------------------------------------------------------------
      SUBROUTINE SUGPRT(IMODEL,IMODIN)
C--------------------------------------------------------------------
C
C     Print SUGRA parameters and results
C     IMODEL = model type for SUGRA
C     IMODIN = input model type to control formatting
C
      IMPLICIT NONE
      COMMON/SSLUN/LOUT
      INTEGER LOUT
      SAVE /SSLUN/
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
      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,ASM3
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3
      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG
      SAVE /SUGPAS/
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
      REAL PI,GPX,SIN2W,ALEMI,AS,TANBQ
      INTEGER IMODEL,J,K,IMODIN
C
C          Entry
C
      PI=4.*ATAN(1.)
      GPX=SQRT(.6)*GSS(1)
      SIN2W=GPX**2/(GSS(2)**2+GPX**2)
      ALEMI=4*PI/GSS(2)**2/SIN2W
      AS=GSS(3)**2/4./PI
      TANBQ=VUQ/VDQ
C
C          Print inputs and GUT couplings for SUGRA/AMSB models
C
      IF(IMODEL.EQ.1.OR.IMODEL.EQ.7) THEN
        IF(IMODEL.EQ.1) THEN
          WRITE(LOUT,1000) XSUGIN(1),XSUGIN(2),XSUGIN(3),XSUGIN(4),
     $    XSUGIN(5),XSUGIN(6)
1000      FORMAT(
     $    ' M_0,  M_(1/2),  A_0,  tan(beta),  sgn(mu),  M_t ='
     $    /4F10.3,2X,F6.1,F10.3)
        ELSE IF (IMODEL.EQ.7) THEN
          WRITE(LOUT,1018) XSUGIN(1),XSUGIN(2),XSUGIN(4),XSUGIN(5),
     $    XSUGIN(6)
1018      FORMAT(
     $    ' M_0,  M_(3/2),  tan(beta),  sgn(mu),  M_t ='
     $    /3F10.3,2X,F6.1,2F10.3)
        END IF
C
C          Write out non-universal GUT scale parameters
        IF(XNUSUG(1).LT.1.E19.OR.XNUSUG(2).LT.1.E19.OR.XNUSUG(3)
     $  .LT.1.E19) THEN
          WRITE(LOUT,1010) XNUSUG(1),XNUSUG(2),XNUSUG(3)
1010      FORMAT(/' M_1(GUT)= ',F8.2,'    M_2(GUT)= ',F8.2,
     $    '    M_3(GUT)= ',F8.2)
        END IF
        IF(XNUSUG(4).LT.1.E19.OR.XNUSUG(5).LT.1.E19.OR.XNUSUG(6)
     $  .LT.1.E19) THEN
          WRITE(LOUT,1011) XNUSUG(4),XNUSUG(5),XNUSUG(6)
1011      FORMAT(/' A_tau(GUT)= ',F8.2,'    A_b(GUT)= ',F8.2,
     $    '    A_t(GUT)= ',F8.2)
        END IF
        IF(XNUSUG(7).LT.1.E19.OR.XNUSUG(8).LT.1.E19) THEN
          WRITE(LOUT,1012) XNUSUG(7),XNUSUG(8)
1012      FORMAT(/' M_Hd(GUT)= ',F8.2,'    M_Hu(GUT)= ',F8.2)
        END IF
        IF (XNUSUG(9).LT.1.E19.OR.XNUSUG(10).LT.1.E19) THEN
          WRITE(LOUT,1013) XNUSUG(9),XNUSUG(10)
1013      FORMAT(/' M_eR(GUT)= ',F8.2,'    M_eL(GUT)= ',F8.2)
        END IF
        IF(XNUSUG(11).LT.1.E19.OR.XNUSUG(12).LT.1.E19.OR.XNUSUG(13)
     $  .LT.1.E19) THEN
          WRITE(LOUT,1014) XNUSUG(11),XNUSUG(12),XNUSUG(13)
1014      FORMAT(' M_dR(GUT)= ',F8.2,'    M_uR(GUT)= ',F8.2,
     $    '    M_uL(GUT)=',F8.2)
        END IF
        IF(XNUSUG(14).LT.1.E19.OR.XNUSUG(15).LT.1.E19) THEN
          WRITE(LOUT,1015) XNUSUG(14),XNUSUG(15)
1015      FORMAT(/' M_tauR(GUT)= ',F8.2,'    M_tauL(GUT)= ',F8.2)
        END IF
        IF(XNUSUG(16).LT.1.E19.OR.XNUSUG(17).LT.1.E19.OR.XNUSUG(18)
     $  .LT.1.E19) THEN
          WRITE(LOUT,1016) XNUSUG(16),XNUSUG(17),XNUSUG(18)
1016      FORMAT(' M_bR(GUT)= ',F8.2,'    M_tR(GUT)= ',F8.2,
     $    '    M_tL(GUT)=',F8.2)
        END IF
        IF(XSUGIN(7).NE.0) THEN
          WRITE(LOUT,1026) XSUGIN(7)
1026      FORMAT(' Q_max= ',E12.4)
        ENDIF
C
C          Right-handed neutrino parameters
        IF (XNRIN(2).LT.1.E19) THEN
          WRITE(LOUT,1017) XNRIN(1),XNRIN(2),XNRIN(3),XNRIN(4),
     $    FNMZ,FNGUT
1017      FORMAT(' Right-handed neutrino parameters:'/
     $    ' M(nu_tau)=',E10.3,'   M(N_R) =',E10.3,
     $    '   A_N=',F8.2,'   M(NRSS)=',F8.2/
     $    ' FN(M_Z)  =',F8.4, '   FN(M_GUT) =',F8.4)
        END IF
C
C          Unification results
        WRITE(LOUT,1001) MGUTSS,GGUTSS,AGUTSS
1001    FORMAT(/' ISASUGRA unification:'/' M_GUT      =',E10.3,
     $  '   g_GUT          =',F5.3,3X,'   alpha_GUT =',F5.3)
        WRITE(LOUT,999) FTGUT,FBGUT,FTAGUT
999     FORMAT(' FT_GUT     =',F6.3,
     $  '       FB_GUT         =',F6.3,3X,'  FL_GUT =',F6.3)
C
C          Print inputs for GMSB models
C
      ELSE IF (IMODEL.EQ.2) THEN
        WRITE(LOUT,1002) (XGMIN(J),J=1,7)
1002    FORMAT(
     $  ' Lambda,  M_mes,  N_5,  tan(beta),  sgn(mu),  M_t,  C_grav='
     $  /2E10.3,2F10.3,2X,F6.1,F10.3,1X,E10.3)
        WRITE(LOUT,1020) (XGMIN(J),J=8,14)
1020    FORMAT(/' GMSB2 model input:'/
     $  ' Rsl,    dmH_d^2,   dmH_u^2,     d_Y,     N5_1,  N5_2,  N5_3='
     $  /F7.3,1X,E10.3,1X,E10.3,1X,E10.3,2X,3F7.3)
        WRITE(LOUT,1003) AMGVSS
1003    FORMAT(/' M(gravitino)=',E10.3)
      END IF
C
C          Weak scale couplings
C
      WRITE(LOUT,1004) ALEMI,SIN2W,AS
1004  FORMAT(/' 1/alpha_em =',F8.2,2X,
     $'   sin**2(thetaw) =',F6.4,2X,'   a_s^DRB   =  ',F5.3)
      WRITE(LOUT,1005) GSS(7),GSS(8),GSS(9)
1005  FORMAT(' M_1        =',F8.2,2X,
     $'   M_2            =',F8.2,'   M_3       =',F8.2)
      WRITE(LOUT,1006) MU,B,HIGFRZ
1006  FORMAT(' mu(Q)      =',F8.2,2X,
     $'   B(Q)           =',F8.2,'   Q         =',F8.2)
      WRITE(LOUT,1007) GSS(13),GSS(14),TANBQ
1007  FORMAT(' M_H1^2     =',E10.3,'   M_H2^2         =',E10.3,
     $' TANBQ     =   ',F6.3)
C
C          Print mass spectrum from ISASUGRA
C
      WRITE(LOUT,2000) MSS(1),MSS(2),MSS(3),MSS(4),MSS(5),MSS(10),
     $MSS(11),MSS(12),MSS(13),MSS(14),MSS(17),MSS(18),MSS(16),
     $MSS(21),MSS(22),MSS(23),MSS(24),MSS(25),MSS(26),MSS(27),
     $MSS(28),MSS(29),MSS(30),MSS(31),MSS(32)
2000  FORMAT(/' ISAJET masses (with signs):'/
     $' M(GL)  =',F9.2/
     $' M(UL)  =',F9.2,'   M(UR)  =',F9.2,'   M(DL)  =',F9.2,
     $'   M(DR) =',F9.2/
     $' M(B1)  =',F9.2,'   M(B2)  =',F9.2,'   M(T1)  =',F9.2,
     $'   M(T2) =',F9.2/
     $' M(SN)  =',F9.2,'   M(EL)  =',F9.2,'   M(ER)  =',F9.2/
     $' M(NTAU)=',F9.2,'   M(TAU1)=',F9.2,'   M(TAU2)=',F9.2/
     $' M(Z1)  =',F9.2,'   M(Z2)  =',F9.2,'   M(Z3)  =',F9.2,
     $'   M(Z4) =',F9.2/
     $' M(W1)  =',F9.2,'   M(W2)  =',F9.2/
     $' M(HL)  =',F9.2,'   M(HH)  =',F9.2,'   M(HA)  =',F9.2,
     $'   M(H+) =',F9.2)
      WRITE(LOUT,2001) THETAT,THETAB,THETAL,ALFAH
2001  FORMAT(/,' theta_t=',F9.4,'   theta_b=',F9.4,
     $'   theta_l=',F9.4,'   alpha_h=',F9.4)
C
C     Write out chargino /neutralino masses/eigenvectors
C
      WRITE(LOUT,3100) AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS
3100  FORMAT(/' NEUTRALINO MASSES (SIGNED) =',4F10.3)
      DO 100 J=1,4
        WRITE(LOUT,3200) J,(ZMIXSS(K,J),K=1,4)
3200    FORMAT(' EIGENVECTOR ',I1,'       =',4F10.5)
100   CONTINUE
      WRITE(LOUT,3300) AMW1SS,AMW2SS
3300  FORMAT(/' CHARGINO MASSES (SIGNED)  =',2F10.3)
      WRITE(LOUT,3400) GAMMAL,GAMMAR
3400  FORMAT(' GAMMAL, GAMMAR             =',2F10.5/)

C
C          Print ISAJET MSSMi equivalent input
C
      WRITE(LOUT,3000)
3000  FORMAT(/' ISAJET equivalent input:')
      WRITE(LOUT,3001) MSS(1),MU,MSS(31),XSUGIN(4)
3001  FORMAT(' MSSMA: ',4F8.2)
      WRITE(LOUT,3002) SQRT(GSS(19)),SQRT(GSS(17)),SQRT(GSS(18)),
     $SQRT(GSS(16)),SQRT(GSS(15))
3002  FORMAT(' MSSMB: ',5F8.2)
      WRITE(LOUT,3003) SIGN(1.,GSS(24))*SQRT(ABS(GSS(24))),
     $SQRT(GSS(22)),SIGN(1.,GSS(23))*SQRT(ABS(GSS(23))),
     $SQRT(GSS(21)),SQRT(GSS(20)),GSS(12),GSS(11),GSS(10)
3003  FORMAT(' MSSMC: ',8F8.2)
      WRITE(LOUT,3004)
3004  FORMAT(' MSSMD: SAME AS MSSMB (DEFAULT)')
      WRITE(LOUT,3005) GSS(7),GSS(8)
3005  FORMAT(' MSSME: ',2F8.2)
      RETURN
      END
CDECK  ID>, SUGRUN. 
      PROGRAM SUGRUN
C
C     Main program to calculate MSSM input parameters for ISAJET
C     from renormalization group equations and supergravity.
C     All external names are of the form SUxxxx.
C     Must link with block data ALDATA.
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
      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,ASM3
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3
      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG
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
C          ISAPW1 is used to check whether ALDATA is loaded
      COMMON/ISAPW/ISAPW1
      CHARACTER*30 ISAPW1
      SAVE /ISAPW/
      CHARACTER*80 FNAME
      REAL M0,MHF,A0,TANB,SGNMU,MT,XLAMGM,XMESGM,XN5GM,AMPL,XCMGV
      INTEGER NSTEP,IMODEL,INUSUG,IMODIN
      INTEGER K,NOUT,IALLOW,IITEST,J
      CHARACTER*40 VERSN,VISAJE
      PARAMETER (NOUT=33)
      INTEGER IDOUT(NOUT)
      CHARACTER*30 ISAPW2
      SAVE ISAPW2
C
      DATA IDOUT/
     $IDTP,ISGL,ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1,ISUPR,ISDNR,
     $ISSTR,ISCHR,ISBT2,ISTP2,ISEL,ISMUL,ISTAU1,ISNEL,ISNML,ISNTL,
     $ISER,ISMUR,ISTAU2,ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,
     $ISHL,ISHH,ISHA,ISHC/
      DATA AMPL/2.4E18/
C          ISAPW2 is used to check whether ALDATA is loaded
      DATA ISAPW2/'ALDATA REQUIRED BY FORTRAN G,H'/
C
C          Initialize
C
      IF(ISAPW1.NE.ISAPW2) THEN
        PRINT*, ' ERROR: BLOCK DATA ALDATA HAS NOT BEEN LOADED.'
        PRINT*, ' ISAJET CANNOT RUN WITHOUT IT.'
        PRINT*, ' PLEASE READ THE FINE MANUAL FOR ISAJET.'
        STOP99
      ENDIF
C
      LOUT=1
      NSTEP=1000
      XNRIN(2)=1.E20
C
      PRINT*,'ENTER output file name (in single quotes):'
      READ*,FNAME
      OPEN(1,FILE=FNAME,STATUS='NEW',FORM='FORMATTED')
      PRINT*,'ENTER 1 for mSUGRA:'
      PRINT*,'ENTER 2 for mGMSB:'
      PRINT*,'ENTER 3 for non-universal SUGRA:'
      PRINT*,'ENTER 4 for SUGRA with truly unified gauge couplings:'
      PRINT*,'ENTER 5 for non-minimal GMSB:'
      PRINT*,'ENTER 6 for SUGRA+right-handed neutrino:'
      PRINT*,'ENTER 7 for anomaly-mediated SUSY breaking:'
      READ*,IMODIN
      IMODEL=IMODIN
      IF (IMODEL.EQ.4) THEN
        IAL3UN=1
        IMODEL=1
      END IF
      IF (IMODEL.EQ.1.OR.IMODEL.EQ.3.OR.IMODEL.EQ.6) THEN
        PRINT*,'ENTER M_0, M_(1/2), A_0, tan(beta), sgn(mu), M_t:'
        READ*,M0,MHF,A0,TANB,SGNMU,MT
        IF (IMODEL.EQ.6) THEN
          IMODEL=1
          PRINT*,' ENTER M(nu_3)[=0], M_Majorana, A_N, M(NRSS)'
          READ*,XNRIN(1),XNRIN(2),XNRIN(3),XNRIN(4)
          GO TO 15
        END IF
        IF (IMODEL.EQ.3) THEN
          IMODEL=1
10        PRINT*,' ENTER 1,...,5 for NUSUGx keyword; 0 to continue:'
          PRINT*,' NUSUG1 = GUT scale gaugino masses'
          PRINT*,' NUSUG2 = GUT scale A terms'
          PRINT*,' NUSUG3 = GUT scale Higgs masses'
          PRINT*,' NUSUG4 = GUT scale 1st/2nd generation masses'
          PRINT*,' NUSUG5 = GUT scale 3rd generation masses'
          PRINT*,' ENTER 6 to activate right-hand neutrino'
          PRINT*,' ENTER 7 to enter alternate high scale Q_max.ne.M_GUT'
          READ*,INUSUG
          IF (INUSUG.EQ.0) THEN
            GO TO 15
          ELSE IF (INUSUG.EQ.1) THEN
            PRINT*,'Enter GUT scale M_1, M_2, M_3:'
            READ*,XNUSUG(1),XNUSUG(2),XNUSUG(3)
C            IF (XNUSUG(3).LE.0.) THEN
C            PRINT*, ' NEGATIVE M_3 IS NOT ALLOWED'
C            STOP 99
C            END IF
          ELSE IF (INUSUG.EQ.2) THEN
            PRINT*,'Enter GUT scale A_t, A_b, A_tau:'
            READ*,XNUSUG(6),XNUSUG(5),XNUSUG(4)
          ELSE IF (INUSUG.EQ.3) THEN
            PRINT*,'Enter GUT scale m_Hd, m_Hu:'
            READ*,XNUSUG(7),XNUSUG(8)
          ELSE IF (INUSUG.EQ.4) THEN
            PRINT*,'Enter GUT scale M(ul), M(dr), M(ur), M(el), M(er):'
            READ*,XNUSUG(13),XNUSUG(11),XNUSUG(12),XNUSUG(10),XNUSUG(9)
          ELSE IF (INUSUG.EQ.5) THEN
            PRINT*,'Enter GUT scale M(tl), M(br), M(tr), M(Ll), M(Lr):'
            READ*,XNUSUG(18),XNUSUG(16),XNUSUG(17),XNUSUG(15),XNUSUG(14)
          ELSE IF (INUSUG.EQ.6) THEN
            PRINT*,' ENTER M(nu_3), M_Majorana, A_N, M(NRSS)'
            READ*,XNRIN(1),XNRIN(2),XNRIN(3),XNRIN(4)
          ELSE IF (INUSUG.EQ.7) THEN
            PRINT*,' ENTER Q_max high scale for SUSY BCs'
            READ*,XSUGIN(7)
          END IF
          GO TO 10
        END IF
      ELSE IF (IMODEL.EQ.2.OR.IMODEL.EQ.5) THEN
          PRINT*,'ENTER Lambda, M_mes, N_5, tan(beta), sgn(mu), ',
     $    'M_t, C_gv:'
          READ*,M0,MHF,A0,TANB,SGNMU,MT,XCMGV
          XGMIN(7)=XCMGV
          XGMIN(8)=1.
          AMGVSS=M0*MHF*XCMGV/SQRT(3.)/AMPL
          IF (IMODEL.EQ.5) THEN
            IMODEL=2
            PRINT*,'Rsl = factor multiplying gaugino masses at M_mes'
            PRINT*,'dmH_d^2, dmH_u^2 = Higgs mass**2 shifts at M_mes'
            PRINT*,'d_Y = mass**2 shifts proportional to Y at M_mes'
            PRINT*,'n5_1,n5_2,n5_3 = n5 values for U(1),SU(2),SU(3)'
            PRINT*,'ENTER Rsl, dmH_d^2, dmH_u^2, d_Y, n5_1, n5_2, n5_3'
            READ*,XGMIN(8),XGMIN(9),XGMIN(10),XGMIN(11),XGMIN(12),
     $      XGMIN(13),XGMIN(14)
            END IF
      ELSE IF (IMODEL.EQ.7) THEN
        PRINT*,'ENTER M_0, M_(3/2), tan(beta), sgn(mu), M_t:'
        READ*,M0,MHF,TANB,SGNMU,MT
        A0=0.
      ELSE
        PRINT*,'Invalid model choice.'
        STOP99
      END IF
C
C          Solve RG equations
C
15    CALL SUGRA(M0,MHF,A0,TANB,SGNMU,MT,IMODEL)
C
C          Print results
C
      VERSN=VISAJE()
      WRITE(LOUT,20) VERSN
20    FORMAT(' ',44('*')/' *',42X,'*'/
     $  ' * ',A40,' *'/
     $  ' *',42X,'*'/' ',44('*')/)
      IF (NOGOOD.EQ.1) THEN
        PRINT*, 'BAD POINT: TACHYONIC PARTICLES!'
        WRITE(LOUT,*) 'BAD POINT: TACHYONIC PARTICLES!'
      ELSE IF (NOGOOD.EQ.2) THEN
        PRINT*, 'BAD POINT: NO EW SYMMETRY BREAKING!'
        WRITE(LOUT,*) 'BAD POINT: NO EW SYMMETRY BREAKING!'
      ELSE IF (NOGOOD.EQ.3) THEN
        PRINT*, 'BAD POINT: M(H_P)^2<0!'
        WRITE(LOUT,*) 'BAD POINT: M(H_P)^2<0!'
      ELSE IF (NOGOOD.EQ.4) THEN
        PRINT*, 'BAD POINT: YUKAWA>10!'
        WRITE(LOUT,*) 'BAD POINT: YUKAWA>10!'
      ELSE IF (NOGOOD.EQ.5.AND.IMODEL.EQ.1) THEN
        PRINT*, 'SUGRA BAD POINT: Z1SS NOT LSP!'
        WRITE(LOUT,*) 'SUGRA BAD POINT: Z1SS NOT LSP!'
      ELSE IF (NOGOOD.EQ.7) THEN
        PRINT*, 'BAD POINT: XT EWSB BAD!'
        WRITE(LOUT,*) 'BAD POINT: XT EWSB BAD!'
      ELSE IF (NOGOOD.EQ.8) THEN
        PRINT*, 'BAD POINT: MHL^2<0!'
        WRITE(LOUT,*) 'BAD POINT: MHL^2<0!'
      ELSE IF (NOGOOD.EQ.-1) THEN
        PRINT*, 'BAD POINT: NO RGE SOLUTION FOUND'
        WRITE(LOUT,*) 'BAD POINT: NO RGE SOLUTION FOUND'
      END IF
      IF (MHPNEG.EQ.1) THEN
        PRINT*, 'BAD POINT: M(H_P)^2<0!!'
        WRITE(LOUT,*) 'BAD POINT: M(H_P)^2<0!!'
        NOGOOD=3
      END IF
      IF(NOGOOD.NE.0) STOP99
      IF(ITACHY.NE.0) THEN
        WRITE(LOUT,*) 'WARNING: TACHYONIC SLEPTONS AT GUT SCALE'
        WRITE(LOUT,*) '         POINT MAY BE INVALID'
      ENDIF
C
C          Print selected model and results
C
      IF(IMODIN.EQ.1) WRITE(LOUT,1001)
1001  FORMAT(//' Minimal supergravity (mSUGRA) model:'/)
      IF(IMODIN.EQ.2) WRITE(LOUT,1002)
1002  FORMAT(//' Minimal gauge mediated (GMSB) model:'/)
      IF(IMODIN.EQ.3) WRITE(LOUT,1003)
1003  FORMAT(//' Non-universal supergravity model:'/)
      IF(IMODIN.EQ.4) WRITE(LOUT,1004)
1004  FORMAT(//' Supergravity model with truly unified couplings:'/)
      IF(IMODIN.EQ.5) WRITE(LOUT,1005)
1005  FORMAT(//' Non-minimal gauge mediated (GMSB) model:'/)
      IF(IMODIN.EQ.6) WRITE(LOUT,1006)
1006  FORMAT(//' Supergravity model with right-handed neutrinos:'/)
      IF(IMODIN.EQ.7) WRITE(LOUT,1007)
1007  FORMAT(//' Anomaly-mediated SUSY breaking model:'/)
      CALL SUGPRT(IMODEL,IMODIN)
C
C          Calculate all masses and decay modes
C
        CALL SSMSSM(XISAIN(1),XISAIN(2),XISAIN(3),
     $ XISAIN(4),XISAIN(5),XISAIN(6),XISAIN(7),XISAIN(8),XISAIN(9),
     $ XISAIN(10),XISAIN(11),XISAIN(12),XISAIN(13),XISAIN(14),
     $ XISAIN(15),XISAIN(16),XISAIN(17),XISAIN(18),XISAIN(19),
     $ XISAIN(20),XISAIN(21),XISAIN(22),XISAIN(23),XISAIN(24),
     $ MT,IALLOW,IMODEL)
C
C          Test parameters
C
      IF(IALLOW.NE.0) THEN
        WRITE(LOUT,2001)
2001    FORMAT(//' MSSM WARNING: Z1SS IS NOT LSP')
      ENDIF
C
      CALL SSTEST(IALLOW)
      IITEST=IALLOW/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,2002)
2002    FORMAT(' MSSM WARNING: Z -> Z1SS Z1SS EXCEEDS BOUND')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,2004)
2004    FORMAT(' MSSM WARNING: Z -> CHARGINOS ALLOWED')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,2008)
2008    FORMAT(' MSSM WARNING: Z -> Z1SS Z2SS TOO BIG')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,2016)
2016    FORMAT(' MSSM WARNING: Z -> SQUARKS, SLEPTONS ALLOWED')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,2032)
2032    FORMAT(' MSSM WARNING: Z -> Z* HL0 EXCEEDS BOUND')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,2064)
2064    FORMAT(' MSSM WARNING: Z -> HL0 HA0 ALLOWED')
      ENDIF
C
      WRITE(LOUT,3600)
3600  FORMAT(//' ISASUSY decay modes:'/
     $' Parent --> daughters',18X,'Width',10X,'Branching ratio'/)
C          Write all modes
      DO 200 J=1,NOUT
        CALL SSPRT(IDOUT(J))
200   CONTINUE
C
      STOP
      END
