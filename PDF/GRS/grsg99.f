*********************************************************************
*                                                                   *
*            LO AND NLO (REAL) PHOTON - PARAMETRIZATIONS            *
*                                                                   *
*                 FOR A DETAILED EXPLANATION SEE :                  *
*                 M. GLUECK, E.REYA, I. Schienbein                  *
*                    Phys. Rev. D60(1999)054019                     *
*                                                                   *
*        PROBLEMS/QUESTIONS TO schien@lpsc.in2p3.fr                 *   
*                                                                   *
*   INPUT:   ISET = number of the parton set :                      *
*              ISET = 1  LEADING ORDER SET                          *
*                        (DATA FILE 'grsg99lo.grid')                *
*              ISET = 2  NEXT-TO-LEADING ORDER MSbar SET            *
*                        (DATA FILE 'grsg99m.grid')                 *
*              ISET = 3  NEXT-TO-LEADING ORDER DISgamma SET         *
*                        (DATA FILE 'grsg99d.grid')                 *
*                                                                   *
*            X  = Bjorken-x       (between  1.E-5 and 1   )         *
*            Q2 = scale in GeV**2 (between  0.4   and 1.E6)         *
*                                                                   *
*   OUTPUT:  UL = X*UP/ALPHA ; DL = x*DOWN/ALPHA                    *   
*            SL = X*STRANGE SEA/ALPHA : GL = X*GLUON/ALPHA          *
*                                                                   *
*            Always x times the distribution is returned            *
*                                                                   *
*   COMMON:  The main program or the calling routine has to have    *
*            a common block  COMMON / INTINI / IINI , and  IINI     *
*            has always to be zero when GRSG99 is called for the    *
*            first time or when 'ISET' has been changed.            *
*                                                                   *
* 16.01.2001                                                        *
*********************************************************************
C... 
      SUBROUTINE GRSG99(ISET,X,Q2,UL,DL,SL,GL)
C...
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*128 mfile
      common/mrinput/mfile
      INTEGER ISET,IINI,IIREAD, ISETSAV 
      PARAMETER (NPART=4, NX=51, NQ=34, NARG=2)
      DIMENSION XUF(NX,NQ), XDF(NX,NQ), XSF(NX,NQ), XGF(NX,NQ)
      DIMENSION PARTON (NPART,NQ,NX-1), QS(NQ), XB(NX)
      DIMENSION XT(NARG), NA(NARG), ARRF(NX+NQ)
      COMMON / INTINI / IINI
      SAVE XUF, XDF, XSF, XGF, NA, ARRF, ISETSAV
*...BJORKEN-X AND Q**2 VALUES OF THE GRID :
       DATA QS / 0.4D0, 0.45D0, 0.5D0, 0.55D0, 0.6D0, 0.75D0,
     1           1.0D0, 1.5D0, 2.0D0, 2.5D0, 4.0D0, 6.4D0, 10.D0,
     2           15.D0, 25.D0, 40.D0, 64.D0, 1.D2, 
     3           1.8D2, 3.2D2, 5.8D2, 1.D3, 1.8D3, 3.2D3, 5.8D3,
     4           1.D4, 1.8D4, 3.2D4, 5.8D4, 1.D5, 1.8D5, 3.2D5, 
     5           5.8D5, 1.D6/
       DATA XB / 1.D-5, 1.5D-5, 2.2D-5, 3.2D-5, 4.8D-5, 7.D-5,
     1           1.D-4, 1.5D-4, 2.2D-4, 3.2D-4, 4.8D-4, 7.D-4,
     2           1.D-3, 1.5D-3, 2.2D-3, 3.2D-3, 4.8D-3, 7.D-3,
     3           1.D-2, 1.5D-2, 2.2D-2, 3.2D-2, 5.0D-2, 7.5D-2,
     4           0.1D0, 0.125D0, 0.15D0, 0.175D0, 0.2D0, 0.225D0,
     5           0.25D0, 0.275D0, 0.3D0, 0.325D0, 0.35D0, 0.375D0, 
     6           0.4D0, 0.45D0, 0.5D0, 0.55D0, 0.6D0, 0.65D0, 0.7D0,  
     7           0.75D0, 0.8D0, 0.85D0, 0.9D0, 0.92,0.95, 0.98, 0.99/

*...CHECK OF X AND Q2 VALUES :
       IF ( (X.LT.1.0D-5) .OR. (X.GT.1.0D0) ) THEN
           WRITE(6,91)
  91       FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE')
           STOP
       ENDIF
       IF ( (Q2.LT.0.4D0) .OR. (Q2.GT.1.D6) ) THEN
           WRITE(6,92)
  92       FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE')
           STOP
       ENDIF
*...INITIALIZATION :
*    SELECTION AND READING OF THE GRID :
*    FILE - NO. = 11 FOR LEADING ORDER ( FIRST NUMBER IN THE
*                                        GRID: 1.884E-03 )
*    FILE - NO. = 22 FOR NEXT-TO-LEADING ORDER MSbar ( FIRST NUMBER IN THE
*                                                      GRID: 1.658E-03 )
*    FILE - NO. = 33 FOR NEXT-TO-LEADING ORDER DISgamma ( FIRST NUMBER IN THE
*                                                         GRID: 1.658E-03 )
      IF (IINI.NE.0) GOTO 16
      IF (ISET.EQ.1) THEN
       ISETSAV=1  
       IIREAD=11
       OPEN(11,FILE=mfile)
      ELSE IF (ISET.EQ.2) THEN
       ISETSAV=2  
       IIREAD=22
       OPEN(22,FILE=mfile)
      ELSE IF (ISET.EQ.3) THEN
       ISETSAV=3  
       IIREAD=33
       OPEN(33,FILE=mfile)
      ELSE
        WRITE(6,93)
  93    FORMAT (2X,'PARTON INTERPOLATION: ISET OUT OF RANGE')
        GOTO 60
      END IF
C     
       DO 15 M = 1, NX-1
       DO 15 N = 1, NQ
       READ(IIREAD,90) PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M),
     1                 PARTON(4,N,M)
  90   FORMAT (4(1PE12.5))
  15   CONTINUE
C
      close(IIREAD)
      IINI = 1
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
        XB0 = XB(IX)
        XB1 = 1.D0-XB(IX)
        XUF(IX,IQ) = PARTON(1,IQ,IX) / (XB1**3 * XB0**0.7)
        XDF(IX,IQ) = PARTON(2,IQ,IX) / (XB1**7 * XB0**0.3)
        XSF(IX,IQ) = PARTON(3,IQ,IX) / (XB1**7 * XB0**0.3)
        XGF(IX,IQ) = PARTON(4,IQ,IX) / (XB1**5 * XB0**0.3)
  20  CONTINUE
        XUF(NX,IQ) = 0.D0
        XDF(NX,IQ) = 0.D0
        XSF(NX,IQ) = 0.D0
        XGF(NX,IQ) = 0.D0
  10  CONTINUE
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = DLOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = DLOG(QS(IQ))
  40  CONTINUE
  16  CONTINUE
      if (ISET .ne. ISETSAV) then 
         print*,'Warning : ISET has been changed'
         print*,'You should reinitialize the GRIDS ! by setting iini=0'
      end if   
*...INTERPOLATION :
      XT(1) = DLOG(X)
      XT(2) = DLOG(Q2)
      UL = FINT(NARG,XT,NA,ARRF,XUF) * (1.D0-X)**3 * X**0.7
      DL = FINT(NARG,XT,NA,ARRF,XDF) * (1.D0-X)**7 * X**0.3
      SL = FINT(NARG,XT,NA,ARRF,XSF) * (1.D0-X)**7 * X**0.3
      GL = FINT(NARG,XT,NA,ARRF,XGF) * (1.D0-X)**5 * X**0.3
 60   RETURN
      END
*
*...CERN LIBRARY ROUTINE E104 (INTERPOLATION) :
*
      FUNCTION FINT(NARG,ARG,NENT,ENT,TABLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ARG(5),NENT(5),ENT(63),TABLE(882)
      DIMENSION D(5),NCOMB(5),IENT(5)
      KD=1
      M=1
      JA=1
         DO 5 I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
         DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
    5 JA=JB+1
      FINT=0.D0
   10 FAC=1.D0
      IADR=KD
      IFADR=1
         DO 15 I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.D0-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
      FINT=FINT+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
         DO 50  K=IL,NARG
   50 NCOMB(K)=1
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END

