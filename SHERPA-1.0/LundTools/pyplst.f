C*********************************************************************
      SUBROUTINE PYPLST
C*********************************************************************
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
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
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYDAT4/,/PYSUBS/,/PYPARS/,
     &/PYINT1/,/PYINT2/,/PYINT5/
C
      WRITE(*,*) 'MSEL = ',MSEL
      DO 9999 ITEST=1,200
         WRITE(*,*) 'MSUB(',ITEST,') = ',MSUB(ITEST)
9999  CONTINUE
      DO 10000 ITEST=1,200
         WRITE(*,*) 'MSTP(',ITEST,') = ',MSTP(ITEST),
     &        ' PARP(',ITEST,') = ',PARP(ITEST),
     &        ';    MSTU(',ITEST,') = ',MSTU(ITEST),
     &        ' PARU(',ITEST,') = ',PARU(ITEST)
10000 CONTINUE
      DO 10001 ITEST=1,200
         WRITE(*,*) 'MSTI(',ITEST,') = ',MSTI(ITEST),
     &        ' PARI(',ITEST,') = ',PARI(ITEST),
     &        ';    MSTJ(',ITEST,') = ',MSTJ(ITEST),
     &        ' PARJ(',ITEST,') = ',PARJ(ITEST) 
10001 CONTINUE
      DO 10002 ITEST=1,100
         WRITE(*,*) 'CKIN(',ITEST,') = ',CKIN(ITEST),
     &        ';    CKIN(',ITEST+100,') = ',CKIN(ITEST+100)
10002 CONTINUE
      DO 10004 ITEST=-40,40
         WRITE(*,*) 'KFIN(1,',ITEST,') = ',KFIN(1,ITEST),
     &        ';    KFIN(2,',ITEST,') = ',KFIN(2,ITEST)
10004 CONTINUE
      DO 10005 ITEST=1,400
         WRITE(*,*) 'MINT(',ITEST,') = ',MSTI(ITEST),
     &        ' VINT(',ITEST,') = ',PARI(ITEST)
10005 CONTINUE
      DO 10006 ITEST=1,250
         WRITE(*,*) 'MDCY(',ITEST,',1) = ',MDCY(ITEST,1),
     &        ' MDCY(',ITEST+250,',1) = ',MDCY(ITEST+250,1)
10006 CONTINUE
      RETURN
      END

