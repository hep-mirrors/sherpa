C HEPEVT with double precision

C******************************************************************
C P Y T H I A     6.137

      SUBROUTINE APYINIT(ECM,DA,DB,DS,IHADR,IISR)

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N) 

      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT4/CHAF(500,2)
      CHARACTER CHAF*16
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)

C...Main parameters of run: c.m. energy and number of events.
      NEV=100

C...Select gamma*/Z0 production process.
      MSEL=0
C...QCD
      MSUB(1)=1
C...WW
C      MSUB(25) = 1

C...Only allow Z0 decay to quarks (i.e. no leptonic final states).
      DO 100 IDC=MDCY(23,2),MDCY(23,2)+MDCY(23,3)-1
        IF(IABS(KFDP(IDC,1)).GE.2) MDME(IDC,1)=MIN(0,MDME(IDC,1))
  100 CONTINUE
C...Only allow W+,W- decay to quarks (i.e. no leptonic final states).
      MSTP(48)=1
      MSTJ(101)=5
C...Parton-Shower
c *as*     MSTP(48) = 0

C...QCD Shower
C      CALL sxp_jt74_0395
c *as*      MSTJ(41) = 1
 
c *as*      MSTJ(101) = 1
C...Hadronization off
c *as*      MSTJ(105) = IHADR
c no hadronisation
      MSTJ(105) = 0; 
      MSTP(111) = 0;
c no photon radiation
      MSTJ(41)  = 1;
c no matrix element corrections in shower
c      MSTJ(47)  = 0;
  
c no isr
      MSTP(11) = 0
      MSTP(61)  = 0



      PARJ(41) = DA
      PARJ(42) = DB
      PARJ(21) = DS

      ECMS = ECM


      CALL PYINIT('CMS','e+','e-',ECM)

      END


      SUBROUTINE FINTERF(NHEPW,ISTHEPW,IDHEPW,JMOHEPW,JDAHEPW,
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

      CALL PYHEPC(2)

      MSTU(70) = 1
      MSTU(71) = NHEP

      CALL PYEXEC
      CALL PYHEPC(1)	

      MSTU(70) = 2
      MSTU(72) = NHEP
cc      CALL PYLIST(1)      

      RETURN
      END


      subroutine fpyshower(iev)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)

c no photon radiation
      MSTJ(41)  = 1;  

      call pyevnt
      call pyhepc(1)	

      return
      end


      SUBROUTINE FHAWFACE(NK,KFJET,MOTHER,INVMOM,PJET,XJET)
CC*** *******************************
CC***
CC*** Translating HEPEVT common block to interal kf::code
CC***
CC*** *******************************

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
C...HEPEVT commonblock.
      PARAMETER (NMXHEP=4000)
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      DOUBLE PRECISION PHEP,VHEP
      SAVE /HEPEVT/

      DOUBLE PRECISION PJET,XJET
      DIMENSION KFJET(2000),PJET(2000,4),XJET(2000,4)
      DIMENSION MOTHER(2000),INVMOM(2,2000)
      INTEGER PYCHGE

      NK=0
      DO I=1,NHEP
        IF(ISTHEP(I).EQ.1.OR.ISTHEP(I).EQ.2.OR.ISTHEP(I).EQ.149) THEN
            NK              = NK+1
            KFJET(NK)       = IDHEP(I)
            MOTHER(NK)      = JMOHEP(1,I)
            INVMOM(1,NK)    = JDAHEP(1,I)
            INVMOM(2,NK)    = JDAHEP(2,I)
            PJET(NK,1)=PHEP(4,I)
            DO J=1,3
              PJET(NK,J+1)=PHEP(J,I)
            ENDDO

            XJET(NK,1)=VHEP(4,I)
            DO J=1,3
              XJET(NK,J+1)=VHEP(J,I)
           ENDDO
        ENDIF
      ENDDO

      RETURN
      END








