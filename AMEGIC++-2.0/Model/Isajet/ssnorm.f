CDECK  ID>, SSNORM.
      SUBROUTINE SSNORM(ID)
C-----------------------------------------------------------------------
C          Normalize branching ratios for ID
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
      INTEGER ID,I
      REAL GAMSUM
C
      GAMSUM=0
      DO 100 I=1,NSSMOD
        IF(ISSMOD(I).EQ.ID) GAMSUM=GAMSUM+GSSMOD(I)
100   CONTINUE
      IF(GAMSUM.EQ.0) RETURN
      DO 200 I=1,NSSMOD
        IF(ISSMOD(I).EQ.ID) BSSMOD(I)=GSSMOD(I)/GAMSUM
200   CONTINUE
      RETURN
      END
