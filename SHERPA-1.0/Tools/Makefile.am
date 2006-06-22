pkglib_LTLIBRARIES       =	libSherpaTools.la

if CLHEP_SUPPORT
HEPMCSOURCES=HepMC_Interface.C
HEPMCHEADERS=HepMC_Interface.H
endif

libSherpaTools_la_SOURCES     =	\
				Kinematics.C \
				Input_Output_Handler.C \
			    	$(HEPMCSOURCES) \
                                Event_Reader.C

pkginclude_HEADERS        =	\
			Kinematics.H \
			Input_Output_Handler.H \
		 	$(HEPMCHEADERS) \
			Event_Reader.H

AM_CXXFLAGS   =	@SHERPAFLAGS@ @CONDITIONAL_ROOTFLAGS@

INCLUDES      =	@SHERPAINCS@ @AMISICINCS@ @APACICINCS@ @MODELINCS@ \
		@PDFINCS@ @ATOOLSINCS@ @CONDITIONAL_CLHEPINCS@ @CONDITIONAL_ROOTINCS@

EXTRA_DIST    =	@SVNINCLUDE@
