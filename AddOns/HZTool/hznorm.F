      SUBROUTINE HZNORM(IFLAG)
      IMPLICIT NONE
#include "hepevtp.inc"
#include "heracmn.inc"
      INTEGER IFLAG
      DOUBLE PRECISION NEVTNRM
      CHARACTER *6 XXXX
      Data XXXX/'HZNORM'/
      IF (IFLAG.eq.1) then
         CALL HCDIR('//PAWC',' ')
         CALL HMDIR(XXXX,'S')
         CALL HCDIR('//HISTO',' ')
         CALL HMDIR(XXXX,'S')
         CALL HBOOK1(1,'SIGMA',1,0.0,1.0,0.0)
         NEVTNRM=0.0
      ELSE IF(IFLAG.EQ.2) THEN
         CALL HCDIR('//PAWC/'//XXXX,' ')
         CALL HFILL(1,0.5,0.0,WTX)
         NEVTNRM=NEVTNRM+1.0
      ELSE IF(IFLAG.EQ.3) THEN
         CALL HCDIR('//PAWC/'//XXXX,' ')
         CALL HNORMA(1,1.0/NEVTNRM)
      ENDIF
      RETURN
      END
