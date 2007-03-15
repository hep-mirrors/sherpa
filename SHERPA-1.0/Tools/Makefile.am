pkglib_LTLIBRARIES       =	libSherpaTools.la

if CLHEP_SUPPORT
HEPMCSOURCES=HepMC_Interface.C
HEPMCHEADERS=HepMC_Interface.H
else
HEPMCEXTRADIST=HepMC_Interface.H HepMC_Interface.C
endif

libSherpaTools_la_SOURCES     =	\
				Kinematics.C \
				Input_Output_Handler.C \
			    	$(HEPMCSOURCES) \
			    	HepEvt_Interface.C \
				Pythia_HepEvt_Translator.C \
				Herwig_HepEvt_Translator.C \
                                Event_Reader.C

pkginclude_HEADERS    =	\
			Kinematics.H \
			Input_Output_Handler.H \
		 	$(HEPMCHEADERS) \
			HepEvt_Interface.H \
			Pythia_HepEvt_Translator.H \
			Herwig_HepEvt_Translator.H \
			Event_Reader.H

AM_CXXFLAGS   =	@SHERPAFLAGS@ @CONDITIONAL_ROOTFLAGS@

INCLUDES      =	@SHERPAINCS@ @AMISICINCS@ @APACICINCS@ @MODELINCS@ \
		@PDFINCS@ @ATOOLSINCS@ @CONDITIONAL_CLHEPINCS@ @CONDITIONAL_ROOTINCS@

EXTRA_DIST    =	@SVNINCLUDE@ $(HEPMCEXTRADIST)
