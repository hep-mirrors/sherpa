C HEPEVT with double precision

      SUBROUTINE INHEPEVT(NEVHEPW,NHEPW,ISTHEPW,IDHEPW,
     &                    JMOHEPW,JDAHEPW,PHEPW,VHEPW)
CC*****************************************************************
CC
CC
CC*****************************************************************

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
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

      NEVHEP=NEVHEPW
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
      RETURN
      END

CC*****************************************************************

      SUBROUTINE OUTHEPEVT()
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
C...HEPEVT commonblock.
      PARAMETER (NMXHEP=4000)
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      DOUBLE PRECISION PHEP,VHEP
      SAVE /HEPEVT/


      PRINT *,'EVENT NUMBER : ',NEVHEP
      PRINT *,'==================================================='
           
      DO IHEP=1,NHEP
        PRINT *,IHEP,IDHEP(IHEP),ISTHEP(IHEP),
     &       JMOHEP(1,IHEP),JMOHEP(2,IHEP),
     &       JDAHEP(1,IHEP),JDAHEP(2,IHEP)
        
        PRINT *,PHEP(1,IHEP),PHEP(2,IHEP),PHEP(3,IHEP),PHEP(4,IHEP),
     &        PHEP(5,IHEP) 
   
        PRINT *,VHEP(1,IHEP),VHEP(2,IHEP),VHEP(3,IHEP),VHEP(4,IHEP)
     &      
    
      ENDDO
      RETURN
      END

C      RETURN
C      END


