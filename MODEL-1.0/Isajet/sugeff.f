CDECK  ID>, SUGEFF. 
C-----------------------------------------------------------------
      SUBROUTINE SUGEFF(G0,SIG1,SIG2)
C-----------------------------------------------------------------
C
C     Compute Higgs mass shift due to 1-loop effective potential
C
      IMPLICIT NONE
      COMMON/SSLUN/LOUT
      INTEGER LOUT
      SAVE /SSLUN/
      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,ASM3
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3
      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG
      SAVE /SUGPAS/
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
      REAL G0(31),SIG1,SIG2
      REAL DT1,DELT1S,SIG1T1,DMSDV2,FB,FT,MST2,MSB2,MSB1,SIG2B1,
     $SIG1B1,SIG2B2,SIG1B2,DB1,SIG1T2,SIG2T1,DELB1S,SIG2T2,MST1,COS2W,
     $MT,PI,COTB,TANB,MB,SIG2B,SIG1B,G,GGP,SIG2T,E,FAC,SIG1T,QS,
     $BETA,SINB,COSB,FAC4,FL,ML,SIG1L,SIG2L,MSL1,MSL2,
     $DELL1S,DL1,SIG1L1,SIG2L1,SIG1L2,SIG2L2,VUQ,VDQ
C
      G=G2
      COS2W=1.-XW
      VUQ=EXP(G0(31))
      VDQ=EXP(G0(30))
      TANB=VUQ/VDQ
      COTB=1./TANB
      BETA=ATAN(TANB)
      SINB=SIN(BETA)
      COSB=COS(BETA)
      PI=4.*ATAN(1.)
      FAC=3./8./PI**2
      FAC4=FAC/3.
      E=EXP(1.)
      QS=HIGFRZ**2
C-----CALCULATE TOP AND BOTTOM CONTRIBUTIONS; USE RUNNING MASSES--
      FL=G0(4)
      FB=G0(5)
      FT=G0(6)
      ML=FL*VDQ
      MB=FB*VDQ
      MT=FT*VUQ
      SIG1T=0.
      SIG2T=-FAC*MT**2*G0(6)**2*LOG(MT**2/E/QS)
      SIG1B=-FAC*MB**2*G0(5)**2*LOG(MB**2/E/QS)
      SIG2B=0.
      SIG1L=-FAC4*ML**2*G0(4)**2*LOG(ML**2/E/QS)
      SIG2L=0.
      GGP=(G**2+GP**2)/2.
      MST1=MSS(12)
      MST2=MSS(13)
      MSB1=MSS(10)
      MSB2=MSS(11)
      MSL1=MSS(21)
      MSL2=MSS(22)
C-----CALCULATE STOP_1 CONTRIBUTION -------------------------------
      DELT1S=(.5*(G0(24)-G0(23))+(8*COS2W-5.)*GGP*
     $      (VDQ**2-VUQ**2)/12.)**2+G0(6)**2*VUQ**2*(G0(12)-MU*COTB)**2
      DT1=.5*(G0(24)-G0(23))+(8*COS2W-5.)*GGP*
     $       (VDQ**2-VUQ**2)/12.
      DMSDV2=GGP/4.-(2*DT1*(8*COS2W-5.)*GGP/12.-
     $        FT**2*MU*(G0(12)*TANB-MU))/2./SQRT(DELT1S)
      SIG1T1=FAC/2.*MST1**2*LOG(MST1**2/E/QS)*DMSDV2
      DMSDV2=-GGP/4.+FT**2-(-2*DT1*(8*COS2W-5.)*GGP/12.+
     $        FT**2*G0(12)*(G0(12)-MU*COTB))/2./SQRT(DELT1S)
      SIG2T1=FAC/2.*MST1**2*LOG(MST1**2/E/QS)*DMSDV2
C-----CALCULATE STOP_2 CONTRIBUTION -------------------------------
      DMSDV2=GGP/4.+(2*DT1*(8*COS2W-5.)*GGP/12.-
     $        FT**2*MU*(G0(12)*TANB-MU))/2./SQRT(DELT1S)
      SIG1T2=FAC/2.*MST2**2*LOG(MST2**2/E/QS)*DMSDV2
      DMSDV2=-GGP/4.+FT**2+(-2*DT1*(8*COS2W-5.)*GGP/12.+
     $        FT**2*G0(12)*(G0(12)-MU*COTB))/2./SQRT(DELT1S)
      SIG2T2=FAC/2.*MST2**2*LOG(MST2**2/E/QS)*DMSDV2
C-----CALCULATE SBOT_1 CONTRIBUTION -------------------------------
      DELB1S=(.5*(G0(24)-G0(22))-(4*COS2W-1.)*GGP*
     $  (VDQ**2-VUQ**2)/12.)**2+G0(5)**2*VDQ**2*(G0(11)-MU*TANB)**2
      DB1=.5*(G0(24)-G0(22))-(4*COS2W-1.)*GGP*
     $       (VDQ**2-VUQ**2)/12.
      DMSDV2=-GGP/4.+FB**2-(-2*DB1*(4*COS2W-1.)*GGP/12.+
     $        FB**2*G0(11)*(G0(11)-MU*TANB))/2./SQRT(DELB1S)
      SIG1B1=FAC/2.*MSB1**2*LOG(MSB1**2/E/QS)*DMSDV2
      DMSDV2=GGP/4.-(2*DB1*(4*COS2W-1.)*GGP/12.-
     $        FB**2*MU*(G0(11)*COTB-MU))/2./SQRT(DELB1S)
      SIG2B1=FAC/2.*MSB1**2*LOG(MSB1**2/E/QS)*DMSDV2
C-----CALCULATE SBOT_2 CONTRIBUTION -------------------------------
      DMSDV2=-GGP/4.+FB**2+(-2*DB1*(4*COS2W-1.)*GGP/12.+
     $        FB**2*G0(11)*(G0(11)-MU*TANB))/2./SQRT(DELB1S)
      SIG1B2=FAC/2.*MSB2**2*LOG(MSB2**2/E/QS)*DMSDV2
      DMSDV2=GGP/4.+(2*DB1*(4*COS2W-1.)*GGP/12.-
     $        FB**2*MU*(G0(11)*COTB-MU))/2./SQRT(DELB1S)
      SIG2B2=FAC/2.*MSB2**2*LOG(MSB2**2/E/QS)*DMSDV2
C-----CALCULATE STAU_1 CONTRIBUTION -------------------------------
      DELL1S=(.5*(G0(21)-G0(20))-(4*COS2W-3.)*GGP*
     $  (VDQ**2-VUQ**2)/4.)**2+G0(4)**2*VDQ**2*(G0(10)-MU*TANB)**2
      DL1=.5*(G0(21)-G0(20))-(4*COS2W-3.)*GGP*
     $       (VDQ**2-VUQ**2)/4.
      DMSDV2=-GGP/4.+FL**2-(-2*DL1*(4*COS2W-3.)*GGP/4.+
     $        FL**2*G0(10)*(G0(10)-MU*TANB))/2./SQRT(DELL1S)
      SIG1L1=FAC4/2.*MSL1**2*LOG(MSL1**2/E/QS)*DMSDV2
      DMSDV2=GGP/4.-(2*DL1*(4*COS2W-3.)*GGP/4.-
     $        FL**2*MU*(G0(10)*COTB-MU))/2./SQRT(DELL1S)
      SIG2L1=FAC4/2.*MSL1**2*LOG(MSL1**2/E/QS)*DMSDV2
C-----CALCULATE STAU_2 CONTRIBUTION -------------------------------
      DMSDV2=-GGP/4.+FL**2+(-2*DL1*(4*COS2W-3.)*GGP/4.+
     $        FL**2*G0(10)*(G0(10)-MU*TANB))/2./SQRT(DELL1S)
      SIG1L2=FAC4/2.*MSL2**2*LOG(MSL2**2/E/QS)*DMSDV2
      DMSDV2=GGP/4.+(2*DL1*(4*COS2W-3.)*GGP/4.-
     $        FL**2*MU*(G0(10)*COTB-MU))/2./SQRT(DELL1S)
      SIG2L2=FAC4/2.*MSL2**2*LOG(MSL2**2/E/QS)*DMSDV2
C-----ADD ALL TERMS ------------------------------------------------
      SIG1=SIG1B+SIG1B1+SIG1B2+SIG1T+SIG1T1+SIG1T2+
     $SIG1L+SIG1L1+SIG1L2
      SIG2=SIG2B+SIG2B1+SIG2B2+SIG2T+SIG2T1+SIG2T2+
     $SIG2L+SIG2L1+SIG2L2
      RETURN
      END
