C HEPEVT with double precision

C******************************************************************

      SUBROUTINE F2HEPEVT(NHEPW,ISTHEPW,IDHEPW,JMOHEPW,JDAHEPW,
     &                   PHEPW,VHEPW)
CC*****************************************************************
CC
CC
CC*****************************************************************

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C...HEPEVT commonblock.
      PARAMETER (NMXHEP=4000)
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      DOUBLE PRECISION PHEP,VHEP
      SAVE /HEPEVT/

      DOUBLE PRECISION PHEPW,VHEPW

      DIMENSION ISTHEPW(NMXHEP),IDHEPW(NMXHEP)
      DIMENSION JMOHEPW(2,NMXHEP),JDAHEPW(2,NMXHEP)
      DIMENSION PHEPW(5,NMXHEP),VHEPW(4,NMXHEP)

      NEVHEP=0
      NHEP=NHEPW
      DO IHEP=1,NHEP
        ISTHEP(IHEP)     = ISTHEPW(IHEP)
        IDHEP(IHEP)      = IDHEPW(IHEP)
        DO J=1,2
          JMOHEP(J,IHEP) = JMOHEPW(J,IHEP)
          JDAHEP(J,IHEP) = JDAHEPW(J,IHEP)
        ENDDO
        DO J=1,4
          VHEP(J,IHEP)   = VHEPW(J,IHEP)
        ENDDO
        DO J=1,5
          PHEP(J,IHEP)   = PHEPW(J,IHEP)
        ENDDO
      ENDDO

      MSTU(70) = 1
      MSTU(71) = NHEP

cc      CALL PYHEPC(1)

      RETURN
      END

      SUBROUTINE F2Parton(NK,KFJET,PJET)
CC*** *******************************
CC***
CC*** Translating HEPEVT common block to interal kf::code
CC***
CC*** *******************************
cc    to be checked

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
C...HEPEVT commonblock.
      PARAMETER (NMXHEP=4000)
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      DOUBLE PRECISION PHEP,VHEP
      SAVE /HEPEVT/

      DOUBLE PRECISION PJET
      DIMENSION KFJET(2000),PJET(2000,4)
      INTEGER PYCHGE

      NK=0
      DO I=1,NHEP
        IF(ISTHEP(I).EQ.1.OR.ISTHEP(I).EQ.149) THEN
            NK=NK+1
            KFJET(NK) = IDHEP(I)
C            IF(ABS(IDHEP(I)).GT.10.AND.ABS(IDHEP(I)).LT.20) 
C     &           KFJET(NK)=5
C            IF(ABS(IDHEP(I)).GT.100.AND.ABS(IDHEP(I)).LT.1000)
C     &           KFJET(NK)=1
C            IF(ABS(IDHEP(I)).GT.1000.AND.ABS(IDHEP(I)).LT.10000)
C     &           KFJET(NK)=3
C            KQ=PYCHGE(IDHEP(I))
C            IF(KQ.EQ.0.AND.ABS(IDHEP(I)).GT.10) KFJET(NK)=KFJET(NK)+1
C            IF(IDHEP(I).EQ.21) KFJET(NK)=21
C            IF(IDHEP(I).EQ.22) KFJET(NK)=7
            PJET(NK,1)=PHEP(4,I)
            DO J=1,3
              PJET(NK,J+1)=PHEP(J,I)
           ENDDO
        ENDIF
      ENDDO

      RETURN
      END



      SUBROUTINE ALTEST(q2,alam,alp)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)

c     save lambda and set shower lambda
      alams=paru(112)
      paru(112)=parj(81)

c     call alphas routine and store results
      alp=pyalps(q2)
      alam=paru(117)

c     do some test outputs
      print  *, 'alpha(',q2,')=',alp,' ; ',alam

c     restore saved lambda
      paru(112)=alams
 
      RETURN
      END






