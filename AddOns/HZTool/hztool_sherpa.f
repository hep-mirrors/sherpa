      SUBROUTINE HZINIT(XSLO,NCH)
      IMPLICIT NONE
      DOUBLE PRECISION XSLO,NCH
      INTEGER NWPAWC,H,N,ISTAT
      CHARACTER*128 HNAME
      COMMON/HZHNAME/HNAME
include "hepevtp.inc"
include "heracmn.inc"
      PARAMETER (NWPAWC=100000)
      COMMON/PAWC/H(NWPAWC)
      WRITE(6,*) ' Now initializing HZTool'
      WRITE(6,*) ' LO cross section [pb]: ',XSLO
      WRITE(6,*) ' Average charged multi: ',NCH
      WRITE(6,*) ' HBook file ',HNAME
      ENDFILE 6
      CALL HLIMIT(NWPAWC)
      CALL HROPEN(1,'HISTO',HNAME,'N',1024,ISTAT)
      GEN='SHA'
      XSEC=XSLO
      CALL HZXXXX(1)
      RETURN
      END

      SUBROUTINE HZEVNT(WGT)
      IMPLICIT NONE
      DOUBLE PRECISION WGT
include "hepevtp.inc"
include "heracmn.inc"
      CALL HZFILHEP
      WTX=WGT
      CALL HZXXXX(2)
      RETURN
      END

      SUBROUTINE HZFINL()
      IMPLICIT NONE
      INTEGER ICYCLE
include "hepevtp.inc"
include "heracmn.inc"
      WTX=1.0
      CALL HZXXXX(3)
      CALL HCDIR('//PAWC',' ')
      CALL HCDIR('//HISTO',' ')
      CALL HROUT(0,ICYCLE,'T')
      CALL HREND('HISTO')
      RETURN
      END
