pkglib_LTLIBRARIES       =	libSherpaTools.la

if HEPMC2_SUPPORT
HEPMC2SOURCES=HepMC2_Interface.C Output_HepMC2_Genevent.C
HEPMC2HEADERS=HepMC2_Interface.H Output_HepMC2_Genevent.H
else
HEPMC2EXTRADIST=HepMC2_Interface.H HepMC2_Interface.C \
Output_HepMC2_Genevent.H Output_HepMC2_Genevent.C
endif

libSherpaTools_la_SOURCES     =	\
				Input_Output_Handler.C \
				Output_Base.C \
				Output_Sherpa.C \
				Output_HepEvt.C \
				Output_D0_HepEvt.C \
			    	$(HEPMC2SOURCES) \
			    	HepEvt_Interface.C \
				Pythia_HepEvt_Translator.C \
                                Event_Reader.C \
				Analysis_Interface.C

localinc_HEADERS      =	\
			Input_Output_Handler.H \
			Output_Base.H \
			Output_Sherpa.H \
			Output_HepEvt.H \
			Output_D0_HepEvt.H \
		 	$(HEPMC2HEADERS) \
			HepEvt_Interface.H \
			Pythia_HepEvt_Translator.H \
			Event_Reader.H \
			Analysis_Interface.H

libSherpaTools_la_CXXFLAGS = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@

libSherpaTools_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@ @CONDITIONAL_HEPMC2INCS@

EXTRA_DIST    = @SVNINCLUDE@ $(HEPMC2EXTRADIST)
