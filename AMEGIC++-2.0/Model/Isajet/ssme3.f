CDECK  ID>, SSME3.
      SUBROUTINE SSME3(KTYP,AM,ZI,ZF)
C
C          Give matrix element data for mode most recently saved by
C          SSSAVE. Call this once for each pole in the matrix element,
C          giving the pole type, mass, and couplings. See /DKYSS3/
C          for more comments.
C
C          Assumes SUSY decay product is always FIRST.
C
      IMPLICIT NONE
C
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
C          Data for SUSY 3-body matrix elements. There is a double
C          pointer structure, first to modes, and then to poles that
C          make up the matrix element for that mode:
C          MELEM=-I in /DKYTAB/ points to the mode information:
C            J1SS3(I) = start of pole list for this mode
C            J2SS3(I) = end of pole list for this mode
C            WTSS3(I) = maximum weight for this mode
C          J1SS3<J<J2SS3 points to the corresponding poles:
C            KSS3(J)    = pole type
C            AMSS3(J)   = pole mass
C            ZISS3(2,J) = initial couplings
C            ZFSS3(2,J) = final couplings
C          For gaugino -> gaugino f fbar, the pole types are
C            KSS3=1: spin-1 pole in f-fbar channel
C            KSS3=2: spin-0 pole in gaugino-f channel
C            KSS3=3: spin-0 pole in gaugino-fbar channel
C            KSS3=4: spin-0 pole in f-fbar channel
C          The two couplings are the coefficients of 1,gamma_5 or of
C          gamma_mu,gamma_mu*gamma_5.
C
      INTEGER MXMSS3,MXPSS3
      PARAMETER (MXMSS3=1000)
      PARAMETER (MXPSS3=2000)
      COMMON/DKYSS3/NMSS3,NPSS3,
     $J1SS3(MXMSS3),J2SS3(MXMSS3),WTSS3(MXMSS3),
     $KSS3(MXPSS3),AMSS3(MXPSS3),ZISS3(2,MXPSS3),ZFSS3(2,MXPSS3)
      INTEGER NMSS3,NPSS3,KSS3,J1SS3,J2SS3
      REAL WTSS3,AMSS3
      COMPLEX ZISS3,ZFSS3
C
      INTEGER KTYP,I
      REAL AM
      COMPLEX ZI(2),ZF(2)
C
C          If last SSSAVE failed, then skip the matrix element
C
      IF(.NOT.LSSMOD) RETURN
C
C          If MSSMOD(NSSMOD)=0, have not booked any poles yet for
C          last mode saved. Increment mode counter, and set initial and
C          final poles to next one.
C
      IF(MSSMOD(NSSMOD).EQ.0) THEN
        NMSS3=NMSS3+1
        IF(NMSS3.GT.MXMSS3) THEN
          WRITE(LOUT,*) 'ERROR IN SSME3...TOO MANY MODES=',NMSS3
          STOP99
        ENDIF
        MSSMOD(NSSMOD)=-NMSS3
        J1SS3(NMSS3)=NPSS3+1
        J2SS3(NMSS3)=NPSS3+1
        WTSS3(NMSS3)=0
      ENDIF
C
C          Add pole to list and set second counter to last pole
C
      NPSS3=NPSS3+1
      IF(NPSS3.GT.MXPSS3) THEN
        WRITE(LOUT,*) 'ERROR IN SSME3...TOO MANY POLES=',NPSS3
        STOP99
      ENDIF
      KSS3(NPSS3)=KTYP
      AMSS3(NPSS3)=AM
      DO 100 I=1,2
        ZISS3(I,NPSS3)=ZI(I)
        ZFSS3(I,NPSS3)=ZF(I)
100   CONTINUE
      J2SS3(NMSS3)=NPSS3
C
      RETURN
      END
