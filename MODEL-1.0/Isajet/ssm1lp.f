CDECK  ID>, SSM1LP. 
      SUBROUTINE SSM1LP(M1,M2,IALLOW)
C-----------------------------------------------------------------------
C
C          Recalculate sparticle masses including
C          radiative corrections
C          from T. Krupovnickas and H. Baer
C
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
      COMMON/SSINF/XLAM
      DOUBLE PRECISION XLAM
C
      REAL PI,SR2,GP,G1,G2,G3,CS2THW,CTHW,AMTLSQ,AMTRSQ,XX,
     $APD,ADMBC,BETA,SINB,COSB,HIGFRZ,
     $PITLTL,PITRTR,PITLTR,PIBLBL,PIBLBR,PIBRBR,PILLLL,PILLLR,PILRLR,
     $PIELEL,PIERER,PINENE
      INTEGER IALLOW
      REAL SIG0L,SIG0R,SIG0S,SIGPL,SIGPR,SIGPS
     $,SSIG0L(4,4),SSIG0R(4,4),SSIG0S(4,4),SSIGPL(2,2),SSIGPR(2,2)
     $,SSIGPS(2,2)
      REAL AR(4,4),NEWAR(4,4),WR(4),TEMP,WORK(4),V,VP,M1,M2
     $,MPPTRE(2,2),MPP(2,2),MPP2(2,2),ZMIX(4,4)
      INTEGER I,J,K,IERR

      PI=4.*ATAN(1.)
      SR2=SQRT(2.)
C     We will use msbar couplings at M_Z here for now
C     to give consistency between MSSM and SUGRA solutions
      G2=SQRT(4.*PI*ALFAEM/SN2THW)
      GP=G2*SQRT(SN2THW/(1.-SN2THW))
      G1=SQRT(5./3.)*GP
      G3=SQRT(4*PI*.118)
      CS2THW=1.-SN2THW
      CTHW=SQRT(CS2THW)
      V=VUQ
      VP=VDQ
      BETA=ATAN(VUQ/VDQ)
      SINB=SIN(BETA)
      COSB=COS(BETA)
      HIGFRZ=SQRT(MAX(AMZ**2,AMTLSS*AMTRSS*SIGN(1.,AMTLSS*AMTRSS)))
      XLAM=DLOG(DBLE(HIGFRZ**2))
C
C     Refill MSS() for input to self energy routines
C
      MSS(1)=AMGLSS
      MSS(2)=AMULSS
      MSS(3)=AMURSS
      MSS(4)=AMDLSS
      MSS(5)=AMDRSS
      MSS(6)=AMSLSS
      MSS(7)=AMSRSS
      MSS(8)=AMCLSS
      MSS(9)=AMCRSS
      MSS(10)=AMB1SS
      MSS(11)=AMB2SS
      MSS(12)=AMT1SS
      MSS(13)=AMT2SS
      MSS(14)=AMN1SS
      MSS(15)=AMN2SS
      MSS(16)=AMN3SS
      MSS(17)=AMELSS
      MSS(18)=AMERSS
      MSS(19)=AMMLSS
      MSS(20)=AMMRSS
      MSS(21)=AML1SS
      MSS(22)=AML2SS
      MSS(23)=AMZ1SS
      MSS(24)=AMZ2SS
      MSS(25)=AMZ3SS
      MSS(26)=AMZ4SS
      MSS(27)=AMW1SS
      MSS(28)=AMW2SS
      MSS(29)=AMHL
      MSS(30)=AMHH
      MSS(31)=AMHA
      MSS(32)=AMHC
C
C     Neutralino masses
C
      AR(1,1)=0.
      AR(1,2)=-TWOM1
      AR(1,3)=-G2*V/SR2
      AR(1,4)=SQRT(3./5.)*G1*V/SR2
      AR(2,1)=-TWOM1
      AR(2,2)=0.
      AR(2,3)=G2*VP/SR2
      AR(2,4)=-SQRT(3./5.)*G1*VP/SR2
      AR(3,1)=-G2*V/SR2
      AR(3,2)=G2*VP/SR2
      AR(3,3)=M2
      AR(3,4)=0.
      AR(4,1)=SQRT(3./5.)*G1*V/SR2
      AR(4,2)=-SQRT(3./5.)*G1*VP/SR2
      AR(4,3)=0.
      AR(4,4)=M1
      DO I=1,4
        DO J=1,4
          SSIG0L(I,J)=SIG0L(AMZ1SS**2,5-I,5-J,G1,G2,CTHW)
          SSIG0R(I,J)=SIG0R(AMZ1SS**2,5-I,5-J,G1,G2,CTHW)
          SSIG0S(I,J)=SIG0S(AMZ1SS**2,5-I,5-J,G1,G2,CTHW)
          NEWAR(I,J)=AR(I,J)
        ENDDO
      ENDDO
      DO I=1,4
        DO J=1,4
          DO K=1,4
            NEWAR(I,J)=NEWAR(I,J)-(SSIG0R(I,K)*AR(K,J)
     $+AR(I,K)*SSIG0L(K,J)
     $+SSIG0R(J,K)*AR(K,I)+AR(J,K)*SSIG0L(K,I))/2.
          ENDDO
          NEWAR(I,J)=NEWAR(I,J)-(SSIG0S(I,J)+SSIG0S(J,I))/2.
        ENDDO
      ENDDO
      CALL EISRS1(4,4,NEWAR,WR,ZMIX,IERR,WORK)
C      IF (IERR.NE.0) THEN
C        WRITE(LOUT,*) 'EISRS1 ERROR IN SSM1LP, IERR=',IERR
C        STOP99
C      END IF
C       Sort eigenvectors and eigenvalues according to masses
      DO I=1,3
        DO J=I+1,4
          IF (ABS(WR(I)).GT.ABS(WR(J))) THEN
            TEMP=WR(J)
            WR(J)=WR(I)
            WR(I)=TEMP
            DO K=1,4
              TEMP=ZMIX(K,J)
              ZMIX(K,J)=ZMIX(K,I)
              ZMIX(K,I)=TEMP
            ENDDO
          END IF
        ENDDO
      ENDDO
      AMZ1SS=WR(1)
      DO I=1,4
        DO J=1,4
          SSIG0L(I,J)=SIG0L(AMZ2SS**2,5-I,5-J,G1,G2,CTHW)
          SSIG0R(I,J)=SIG0R(AMZ2SS**2,5-I,5-J,G1,G2,CTHW)
          SSIG0S(I,J)=SIG0S(AMZ2SS**2,5-I,5-J,G1,G2,CTHW)
          NEWAR(I,J)=AR(I,J)
        ENDDO
      ENDDO
      DO I=1,4
        DO J=1,4
          DO K=1,4
            NEWAR(I,J)=NEWAR(I,J)-(SSIG0R(I,K)*AR(K,J)
     $+AR(I,K)*SSIG0L(K,J)
     $+SSIG0R(J,K)*AR(K,I)+AR(J,K)*SSIG0L(K,I))/2.
          ENDDO
          NEWAR(I,J)=NEWAR(I,J)-(SSIG0S(I,J)+SSIG0S(J,I))/2.
        ENDDO
      ENDDO
      CALL EISRS1(4,4,NEWAR,WR,ZMIX,IERR,WORK)
C      IF (IERR.NE.0) THEN
C        WRITE(LOUT,*) 'EISRS1 ERROR IN SSM1LP, IERR=',IERR
C        STOP99
C      END IF
C       Sort eigenvectors and eigenvalues according to masses
      DO I=1,3
        DO J=I+1,4
          IF (ABS(WR(I)).GT.ABS(WR(J))) THEN
            TEMP=WR(J)
            WR(J)=WR(I)
            WR(I)=TEMP
            DO K=1,4
              TEMP=ZMIX(K,J)
              ZMIX(K,J)=ZMIX(K,I)
              ZMIX(K,I)=TEMP
            ENDDO
          END IF
        ENDDO
      ENDDO
      AMZ2SS=WR(2)
      DO I=1,4
        DO J=1,4
          SSIG0L(I,J)=SIG0L(AMZ3SS**2,5-I,5-J,G1,G2,CTHW)
          SSIG0R(I,J)=SIG0R(AMZ3SS**2,5-I,5-J,G1,G2,CTHW)
          SSIG0S(I,J)=SIG0S(AMZ3SS**2,5-I,5-J,G1,G2,CTHW)
          NEWAR(I,J)=AR(I,J)
        ENDDO
      ENDDO
      DO I=1,4
        DO J=1,4
          DO K=1,4
            NEWAR(I,J)=NEWAR(I,J)-(SSIG0R(I,K)*AR(K,J)
     $+AR(I,K)*SSIG0L(K,J)
     $+SSIG0R(J,K)*AR(K,I)+AR(J,K)*SSIG0L(K,I))/2.
          ENDDO
          NEWAR(I,J)=NEWAR(I,J)-(SSIG0S(I,J)+SSIG0S(J,I))/2.
        ENDDO
      ENDDO
      CALL EISRS1(4,4,NEWAR,WR,ZMIX,IERR,WORK)
C      IF (IERR.NE.0) THEN
C        WRITE(LOUT,*) 'EISRS1 ERROR IN SSM1LP, IERR=',IERR
C        STOP99
C      END IF
C       Sort eigenvectors and eigenvalues according to masses
      DO I=1,3
        DO J=I+1,4
          IF (ABS(WR(I)).GT.ABS(WR(J))) THEN
            TEMP=WR(J)
            WR(J)=WR(I)
            WR(I)=TEMP
            DO K=1,4
              TEMP=ZMIXSS(K,J)
              ZMIX(K,J)=ZMIX(K,I)
              ZMIX(K,I)=TEMP
            ENDDO
          END IF
        ENDDO
      ENDDO
      AMZ3SS=WR(3)
      DO I=1,4
        DO J=1,4
          SSIG0L(I,J)=SIG0L(AMZ4SS**2,5-I,5-J,G1,G2,CTHW)
          SSIG0R(I,J)=SIG0R(AMZ4SS**2,5-I,5-J,G1,G2,CTHW)
          SSIG0S(I,J)=SIG0S(AMZ4SS**2,5-I,5-J,G1,G2,CTHW)
          NEWAR(I,J)=AR(I,J)
        ENDDO
      ENDDO
      DO I=1,4
        DO J=1,4
          DO K=1,4
            NEWAR(I,J)=NEWAR(I,J)-(SSIG0R(I,K)*AR(K,J)
     $+AR(I,K)*SSIG0L(K,J)
     $+SSIG0R(J,K)*AR(K,I)+AR(J,K)*SSIG0L(K,I))/2.
          ENDDO
          NEWAR(I,J)=NEWAR(I,J)-(SSIG0S(I,J)+SSIG0S(J,I))/2.
        ENDDO
      ENDDO
      CALL EISRS1(4,4,NEWAR,WR,ZMIX,IERR,WORK)
C      IF (IERR.NE.0) THEN
C        WRITE(LOUT,*) 'EISRS1 ERROR IN SSM1LP, IERR=',IERR
C        STOP99
C      END IF
C       Sort eigenvectors and eigenvalues according to masses
      DO I=1,3
        DO J=I+1,4
          IF (ABS(WR(I)).GT.ABS(WR(J))) THEN
            TEMP=WR(J)
            WR(J)=WR(I)
            WR(I)=TEMP
            DO K=1,4
              TEMP=ZMIXSS(K,J)
              ZMIX(K,J)=ZMIX(K,I)
              ZMIX(K,I)=TEMP
            ENDDO
          END IF
        ENDDO
      ENDDO
      AMZ4SS=WR(4)
C
C     Charginos
C
      DO I=1,2
        DO J=1,2
          IF(I.NE.J) THEN
            SSIGPL(I,J)=-SIGPL(AMW1SS**2,J,I,G1,G2,CTHW)
            SSIGPR(I,J)=-SIGPR(AMW1SS**2,J,I,G1,G2,CTHW)
            SSIGPS(I,J)=-SIGPS(AMW1SS**2,J,I,G1,G2,CTHW)
          ELSE
            SSIGPL(I,J)=SIGPL(AMW1SS**2,I,J,G1,G2,CTHW)
            SSIGPR(I,J)=SIGPR(AMW1SS**2,I,J,G1,G2,CTHW)
            SSIGPS(I,J)=SIGPS(AMW1SS**2,I,J,G1,G2,CTHW)
          ENDIF
          MPP(I,J)=0.
        ENDDO
      ENDDO
      MPPTRE(1,1)=M2
      MPPTRE(1,2)=SR2*AMW*SIN(BETA)
      MPPTRE(2,1)=SR2*AMW*COS(BETA)
      MPPTRE(2,2)=TWOM1
      DO I=1,2
        DO J=1,2
          MPP(I,J)=MPPTRE(I,J)
        ENDDO
      ENDDO
      DO I=1,2
        DO J=1,2
          DO K=1,2
            MPP(I,J)=MPP(I,J)-SSIGPR(I,K)*MPPTRE(K,J)
     $-MPPTRE(I,K)*SSIGPL(K,J)
          ENDDO
          MPP(I,J)=MPP(I,J)-SSIGPS(I,J)
        ENDDO
      ENDDO
      DO I=1,2
        DO J=1,2
          MPP2(I,J)=MPP(I,1)*MPP(J,1)+MPP(I,2)*MPP(J,2)
        ENDDO
      ENDDO
      AMW1SS=SIGN(1.,AMW1SS)*SQRT(ABS(MIN((MPP2(1,1)+MPP2(2,2)
     $-SQRT((MPP2(1,1)-MPP2(2,2))**2+4.*MPP2(1,2)*MPP2(2,1)))/2.
     $,(MPP2(1,1)+MPP2(2,2)
     $+SQRT((MPP2(1,1)-MPP2(2,2))**2+4.*MPP2(1,2)*MPP2(2,1)))/2.)))
      DO I=1,2
        DO J=1,2
          IF(I.NE.J) THEN
            SSIGPL(I,J)=-SIGPL(AMW2SS**2,J,I,G1,G2,CTHW)
            SSIGPR(I,J)=-SIGPR(AMW2SS**2,J,I,G1,G2,CTHW)
            SSIGPS(I,J)=-SIGPS(AMW2SS**2,J,I,G1,G2,CTHW)
          ELSE
            SSIGPL(I,J)=SIGPL(AMW2SS**2,I,J,G1,G2,CTHW)
            SSIGPR(I,J)=SIGPR(AMW2SS**2,I,J,G1,G2,CTHW)
            SSIGPS(I,J)=SIGPS(AMW2SS**2,I,J,G1,G2,CTHW)
          ENDIF
          MPP(I,J)=0.
        ENDDO
      ENDDO
      DO I=1,2
        DO J=1,2
          MPP(I,J)=MPPTRE(I,J)
        ENDDO
      ENDDO
      DO I=1,2
        DO J=1,2
          DO K=1,2
            MPP(I,J)=MPP(I,J)-SSIGPR(I,K)*MPPTRE(K,J)
     $-MPPTRE(I,K)*SSIGPL(K,J)
          ENDDO
          MPP(I,J)=MPP(I,J)-SSIGPS(I,J)
        ENDDO
      ENDDO
      DO I=1,2
        DO J=1,2
          MPP2(I,J)=MPP(I,1)*MPP(J,1)+MPP(I,2)*MPP(J,2)
        ENDDO
      ENDDO
      AMW2SS=SIGN(1.,AMW2SS)*SQRT(ABS(MAX((MPP2(1,1)+MPP2(2,2)
     $-SQRT((MPP2(1,1)-MPP2(2,2))**2+4.*MPP2(1,2)*MPP2(2,1)))/2.
     $,(MPP2(1,1)+MPP2(2,2)
     $+SQRT((MPP2(1,1)-MPP2(2,2))**2+4.*MPP2(1,2)*MPP2(2,1)))/2.)))
C
C     Do third generation squarks and sleptons
C
C     Top squarks
      AMTLSQ=SIGN(1.,AMTLSS)*AMTLSS**2
      AMTRSQ=SIGN(1.,AMTRSS)*AMTRSS**2
      XX=(AMGLSS/AMULSS)**2
      APD=AMTLSQ+AMTRSQ+2.*MTQ**2+AMZ**2*(COSB**2-SINB**2)/2.
     $-PITLTL(AMTLSS**2,G1,G2,G3,CTHW)-
     $PITRTR(AMTLSS**2,G1,G2,G3,CTHW)
      ADMBC=(AMTLSQ+MTQ**2+AMZ**2*(COSB**2-SINB**2)*(.5-2./3.*SN2THW)
     $-PITLTL(AMTLSS**2,G1,G2,G3,CTHW))*(AMTRSQ+MTQ**2+AMZ**2
     $*(COSB**2-SINB**2)*2./3.*SN2THW
     $-PITRTR(AMTLSS**2,G1,G2,G3,CTHW))-(MTQ*(-AAT-TWOM1*COSB/SINB)
     $-PITLTR(AMTLSS**2,G1,G2,G3,CTHW))**2
      AMT1SS=SQRT(ABS((APD-SQRT(ABS(APD**2-4.*ADMBC)))/2.))
      AMT2SS=SQRT(ABS((APD+SQRT(ABS(APD**2-4.*ADMBC)))/2.))
      THETAT=ATAN((SQRT(ABS((APD-SQRT(ABS(APD**2-4.*ADMBC)))/2.))**2
     $-MTQ**2+AMZ**2*(COSB**2-SINB**2)
     $*(-.5+2*SN2THW/3.)-AMTLSQ+PITLTL(AMTLSS**2,G1,G2,G3,CTHW))
     $/(MTQ*(TWOM1*COSB/SINB+AAT)+PITLTR(AMTLSS**2,G1,G2,G3,CTHW)))
C
       APD=AMBLSS**2+MBQ**2+AMZ**2*(COSB**2-SINB**2)*(-.5+SN2THW/3.)
     $+AMBRSS**2+MBQ**2+AMZ**2*(COSB**2-SINB**2)*SN2THW/3.
     $-PIBLBL(AMBLSS**2,G1,G2,G3,CTHW)-
     $PIBRBR(AMBLSS**2,G1,G2,G3,CTHW)
      ADMBC=(AMBLSS**2+MBQ**2+AMZ**2*(COSB**2-SINB**2)
     $*(-.5+SN2THW/3.)-PIBLBL(AMBLSS**2,G1,G2,G3,CTHW))*
     $(AMBRSS**2+MBQ**2+AMZ**2*(COSB**2-SINB**2)*SN2THW/3.-
     $PIBRBR(AMBLSS**2,G1,G2,G3,CTHW))-(MBQ*(-AAB-TWOM1*SINB/COSB)-
     $PIBLBR(AMBLSS**2,G1,G2,G3,CTHW))**2
      AMB1SS=SQRT(ABS((APD-SQRT(ABS(APD**2-4.*ADMBC)))/2.))
      AMB2SS=SQRT(ABS((APD+SQRT(ABS(APD**2-4.*ADMBC)))/2.))
      THETAB=ATAN((ABS((APD-SQRT(ABS(APD**2-4.*ADMBC)))/2.)-MBQ**2
     $+AMZ**2*(COSB**2-SINB**2)*(.5-SN2THW/3.)
     $-AMBLSS**2+PIBLBL(AMBLSS**2,G1,G2,G3,CTHW))
     $/(MBQ*(TWOM1*SINB/COSB+AAB)+PIBLBR(AMBLSS**2,G1,G2,G3,CTHW)))
C
      APD=AMLLSS**2+MLQ**2+AMZ**2*(COSB**2-SINB**2)*(-.5+SN2THW)
     $+AMLRSS**2+MLQ**2+AMZ**2*(COSB**2-SINB**2)*(-SN2THW)
     $-PILLLL(AMLLSS**2,G1,G2,G3,CTHW)-PILRLR(AMLLSS**2,G1,G2,G3,CTHW)
      ADMBC=(AMLLSS**2+MLQ**2+AMZ**2*(COSB**2-SINB**2)*(-.5+SN2THW)
     $-PILLLL(AMLLSS**2,G1,G2,G3,CTHW))*(AMLRSS**2+MLQ**2
     $+AMZ**2*(COSB**2-SINB**2)*(-SN2THW)-
     $PILRLR(AMLLSS**2,G1,G2,G3,CTHW))-(MLQ*(-AAL-TWOM1*SINB/COSB)
     $-PILLLR(AMLLSS**2,G1,G2,G3,CTHW))**2
      AML1SS=SQRT(ABS((APD-SQRT(ABS(APD**2-4.*ADMBC)))/2.))
      AML2SS=SQRT(ABS((APD+SQRT(ABS(APD**2-4.*ADMBC)))/2.))
      THETAL=ATAN((ABS((APD-SQRT(ABS(APD**2-4.*ADMBC)))/2.)
     $-MLQ**2+AMZ**2*(COSB**2-SINB**2)*(.5-SN2THW)-
     $AMLLSS**2+PILLLL(AMLLSS**2,G1,G2,G3,CTHW))
     $/(MLQ*(TWOM1*SINB/COSB+AAL)+PILLLR(AMLLSS**2,G1,G2,G3,CTHW)))
C
      AMULSS=SQRT(AMULSS**2*(1.+G3**2/6./PI**2*(1.+3.*XX
     $+(XX-1.)**2*LOG(ABS(XX-1.))-XX**2*LOG(XX))))
      XX=(AMGLSS/AMURSS)**2
      AMURSS=SQRT(AMURSS**2*(1.+G3**2/6./PI**2*(1.+3.*XX
     $+(XX-1.)**2*LOG(ABS(XX-1.))-XX**2*LOG(XX))))
      XX=(AMGLSS/AMDLSS)**2
      AMDLSS=SQRT(AMDLSS**2*(1.+G3**2/6./PI**2*(1.+3.*XX
     $+(XX-1.)**2*LOG(ABS(XX-1.))-XX**2*LOG(XX))))
      XX=(AMGLSS/AMDRSS)**2
      AMDRSS=SQRT(AMDRSS**2*(1.+G3**2/6./PI**2*(1.+3.*XX
     $+(XX-1.)**2*LOG(ABS(XX-1.))-XX**2*LOG(XX))))
      XX=(AMGLSS/AMCLSS)**2
      AMCLSS=SQRT(AMCLSS**2*(1.+G3**2/6./PI**2*(1.+3.*XX
     $+(XX-1.)**2*LOG(ABS(XX-1.))-XX**2*LOG(XX))))
      XX=(AMGLSS/AMCRSS)**2
      AMCRSS=SQRT(AMCRSS**2*(1.+G3**2/6./PI**2*(1.+3.*XX
     $+(XX-1.)**2*LOG(ABS(XX-1.))-XX**2*LOG(XX))))
      XX=(AMGLSS/AMSLSS)**2
      AMSLSS=SQRT(AMSLSS**2*(1.+G3**2/6./PI**2*(1.+3.*XX
     $+(XX-1.)**2*LOG(ABS(XX-1.))-XX**2*LOG(XX))))
      XX=(AMGLSS/AMSRSS)**2
      AMSRSS=SQRT(AMSRSS**2*(1.+G3**2/6./PI**2*(1.+3.*XX
     $+(XX-1.)**2*LOG(ABS(XX-1.))-XX**2*LOG(XX))))
      AMELSS=SQRT(AMELSS**2-PIELEL(AMELSS**2,G1,G2,G3,CTHW))
      AMERSS=SQRT(AMERSS**2-PIERER(AMERSS**2,G1,G2,G3,CTHW))
      AMN1SS=SQRT(AMN1SS**2-PINENE(AMN1SS**2,G1,G2,G3,CTHW))
      AMMLSS=SQRT(AMMLSS**2-PIELEL(AMMLSS**2,G1,G2,G3,CTHW))
      AMMRSS=SQRT(AMMRSS**2-PIERER(AMMRSS**2,G1,G2,G3,CTHW))
      AMN2SS=SQRT(AMN2SS**2-PINENE(AMN2SS**2,G1,G2,G3,CTHW))
      AMN3SS=SQRT(AMN3SS**2-PINENE(AMN3SS**2,G1,G2,G3,CTHW))
C
      RETURN
      END
