pkglib_LTLIBRARIES       =	libSherpaTools.la

if HEPMC2_SUPPORT
HEPMC2SOURCES=HepMC2_Interface.C
HEPMC2HEADERS=HepMC2_Interface.H
else
HEPMC2EXTRADIST=HepMC2_Interface.H HepMC2_Interface.C
endif

SVNTAG = SHERPA/Tools
include ../../svn.make

libSherpaTools_la_SOURCES     =	\
				SVN_Info.C \
				Output_Base.C \
				Event_Reader_Base.C \
				HepEvt_Interface.C \
			    	$(HEPMC2SOURCES) \
				Pythia_HepEvt_Translator.C \
				Output_HepEvt.C \
				Output_LHEF.C \
				Analysis_Interface.C

localinc_HEADERS      =	\
			Definitions.H \
			Output_Base.H \
			Event_Reader_Base.H \
			HepEvt_Interface.H \
		 	$(HEPMC2HEADERS) \
			Pythia_HepEvt_Translator.H \
			Output_HepEvt.H \
			Output_LHEF.H \
			Analysis_Interface.H

libSherpaTools_la_LIBADD = @CONDITIONAL_HEPMC2LIBS@

libSherpaTools_la_CXXFLAGS = $(AM_CXXFLAGS)

libSherpaTools_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_HEPMC2INCS@

EXTRA_DIST    = @SVNINCLUDE@ $(HEPMC2EXTRADIST)
