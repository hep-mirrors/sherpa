C---HERWIG 6
      SUBROUTINE HWCINT(CPTMIN,PRNSET) 
      DOUBLE PRECISION CPTMIN
      LOGICAL PRNSET
C---COMMON BLOCKS ARE INCLUDED AS FILE herwig6500.inc
      INCLUDE 'herwig6500.inc'
      EXTERNAL HWUDAT
C---INITIALISE OTHER COMMON BLOCKS
      CALL HWIGIN
C---USER CAN RESET PARAMETERS AT
C   THIS POINT, OTHERWISE DEFAULT
C   VALUES IN HWIGIN WILL BE USED.
      PTMIN=CPTMIN
C---DON'T INCLUDE SPIN CORRELATIONS
      SYSPIN=.FALSE.
      PRNDEF=.FALSE.
      IF (PRNSET) PRNDEF=.TRUE.
C---COMPUTE PARAMETER-DEPENDENT CONSTANTS
      CALL HWUINC
C---CALL HWUSTA TO MAKE ANY PARTICLE STABLE
C      CALL HWUSTA('PI0     ')
C---USER'S INITIAL CALCULATIONS
      CALL HWABEG
C---INITIALISE ELEMENTARY PROCESS
      CALL HWEINI
      RETURN
      END
C---SINGLE EVENT GENERATION
      SUBROUTINE HWCGSE
C---INITIALISE EVENT
      CALL HWUINE
C---GENERATE HARD SUBPROCESS
      CALL HWEPRO
C---GENERATE PARTON CASCADES
      CALL HWBGEN
C---DO HEAVY OBJECT DECAYS
      CALL HWDHOB
C---DO CLUSTER FORMATION
      CALL HWCFOR
C---DO CLUSTER DECAYS
      CALL HWCDEC
C---DO UNSTABLE PARTICLE DECAYS
      CALL HWDHAD
C---DO HEAVY FLAVOUR HADRON DECAYS
      CALL HWDHVY
C---ADD SOFT UNDERLYING EVENT IF NEEDED
      CALL HWMEVT
C---FINISH EVENT
      CALL HWUFNE
      RETURN
      END
C---TERMINATE ELEMENTARY PROCESS
      SUBROUTINE HWCTRM
      CALL HWEFIN
C---USER'S TERMINAL CALCULATIONS
      CALL HWAEND
      STOP
      END
C---DUMMY ROUTINES
      SUBROUTINE HWAINI
      RETURN
      END
      SUBROUTINE HWABEG
      RETURN
      END
      SUBROUTINE HWAEND
      RETURN
      END
      SUBROUTINE HWHVVJ
      RETURN
      END
C*********************************************************************
 
C...UPEVNT
C...Dummy routine, to be replaced by a user implementing external
C...processes. Depending on cross section model chosen, it either has
C...to generate a process of the type IDPRUP requested, or pick a type
C...itself and generate this event. The event is to be stored in the
C...HEPEUP commonblock, including (often) an event weight.
 
      SUBROUTINE UPEVNT
 
      RETURN
      END
 
C*********************************************************************
 
C...PDFSET
C...Dummy routine, to be removed when PDFLIB is to be linked.
 
      SUBROUTINE PDFSET(PARM,VALUE)
 
      RETURN
      END
 
C*********************************************************************
 
C...STRUCTM
C...Dummy routine, to be removed when PDFLIB is to be linked.
 
      SUBROUTINE STRUCTM(XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
 
      RETURN
      END
 
C*********************************************************************
 
C...UPINIT
C...Dummy routine, to be replaced by a user implementing external
C...processes. Is supposed to fill the HEPRUP commonblock with info
C...on incoming beams and allowed processes.
 
      SUBROUTINE UPINIT
 
      RETURN
      END
 
C*********************************************************************
