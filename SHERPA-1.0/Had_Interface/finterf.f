C HEPEVT with double precision

C******************************************************************
C P Y T H I A     6.137

      SUBROUTINE APYINIT(ECM,DA,DB,DS,IHADR,IISR)

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N) 
c *as* new common blocks:
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
c *as* old common blocks:
c      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
c      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
c      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
c      COMMON/PYDAT4/CHAF(500,2)
c      CHARACTER CHAF*16
c      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
c      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
c      COMMON/PYINT1/MINT(400),VINT(400)
c      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
c      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)

CC     SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYDAT4/,/PYSUBS/,/PYPARS/,
CC     &/PYINT1/,/PYINT2/,/PYINT5/

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
         print *,'channel',IDC,KFDP(IDC,1)
        IF(IABS(KFDP(IDC,1)).GE.2) MDME(IDC,1)=MIN(0,MDME(IDC,1))
        print *,' switch',MDME(IDC,1)
c *as*        IF(IABS(KFDP(IDC,1)).GE.6) MDME(IDC,1)=MIN(0,MDME(IDC,1))
  100 CONTINUE
C...Only allow W+,W- decay to quarks (i.e. no leptonic final states).
c *as*      DO 200 IDC=MDCY(24,2),MDCY(24,2)+MDCY(24,3)-1
c *as*        IF(IABS(KFDP(IDC,1)).GE.6) MDME(IDC,1)=MIN(0,MDME(IDC,1))
c *as*  200 CONTINUE   
C...Initialize.
C...Matrixelement       *as* comment 2 lines!
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



c *as* Hadronisation parameter
c      print *,' parj(41) =',parj(41)
      PARJ(41) = DA
c      print *,' parj(41) =',parj(41)

c      print *,' parj(42) =',parj(42)
      PARJ(42) = DB
c      print *,' parj(42) =',parj(42)

c      print *,' parj(21) =',parj(21)
      PARJ(21) = DS
c      print *,' parj(21) =',parj(21)

      ECMS = ECM

c *as*      MSTU(22) = 100000

C>>>ME
C      MSTP(71) = 0
C...ISR off
C      MSTJ(107) = 0

C      MSTJ(101) = -2
C      MSTP(48) = 1


C...DELPHI
C      PARJ(122) = 0.163
C      PARJ(129) = 0.0025
C      PARJ(41)  = 0.903
C      PARJ(42)  = 0.58
C      PARJ(21)  = 0.477
C     PARJ(2)   = 0.277
C      PARJ(1)   = 0.087
C      PARJ(25)  = 0.65
C      PARJ(26)  = 0.23  
C      MSTJ(11) = 3 
C      MSTJ(12) = 3
C...test
C      MSTJ(101)= -1
C      MSTJ(111)= 1

C      MSTJ(46) = 3

C      print *,"Init"
C      MSTP(111) = 0
C      MSTP(71)  = 0
C     ISR
C...Electron carry's the whole energy

C ISR TEST

C NO ISR
C      MSTP(11) = 0
C      MSTP(61)  = 0
C ISR
c *as*      MSTP(11)  = IISR
c *as*      MSTP(61)  = IISR

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

c     do a pythia event      
      call pyevnt

c     do a print out
      if (iev.le.10) call pylist(1)

c     convert pylist to HEPEVT
      call pyhepc(1)	

      return
      end


      SUBROUTINE FHAWFACE(NK,KFJET,PJET)
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








