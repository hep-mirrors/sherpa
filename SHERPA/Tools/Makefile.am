pkglib_LTLIBRARIES       =	libSherpaTools.la

if HEPMC2_SUPPORT
HEPMC2SOURCES=HepMC2_Interface.C Output_HepMC2_Genevent.C Output_HepMC2_Short.C
HEPMC2HEADERS=HepMC2_Interface.H Output_HepMC2_Genevent.H Output_HepMC2_Short.H
else
HEPMC2EXTRADIST=HepMC2_Interface.H HepMC2_Interface.C \
Output_HepMC2_Genevent.H Output_HepMC2_Genevent.C \
Output_HepMC2_Short.H Output_HepMC2_Short.C
endif

SVNTAG = SHERPA/Tools
include ../../svn.make

libSherpaTools_la_SOURCES     =	\
				SVN_Info.C \
				Input_Output_Handler.C \
				Output_Base.C \
				Output_Sherpa.C \
				Output_HepEvt.C \
				Output_RootNtuple.C \
				Output_D0_HepEvt.C \
				Output_LHEF.C \
			    	$(HEPMC2SOURCES) \
			    	HepEvt_Interface.C \
				Pythia_HepEvt_Translator.C \
                                Event_Reader.C \
				Analysis_Interface.C \
				RootNtuple_Reader.C

localinc_HEADERS      =	\
			Input_Output_Handler.H \
			Output_Base.H \
			Output_Sherpa.H \
			Output_HepEvt.H \
			Output_D0_HepEvt.H \
			Output_RootNtuple.H \
			Output_LHEF.H \
		 	$(HEPMC2HEADERS) \
			HepEvt_Interface.H \
			Pythia_HepEvt_Translator.H \
			Event_Reader_Base.H \
			Event_Reader.H \
			Analysis_Interface.H \
			RootNtuple_Reader.H

libSherpaTools_la_LIBADD = @CONDITIONAL_ROOTLIBS@

libSherpaTools_la_CXXFLAGS = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@

libSherpaTools_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@ @CONDITIONAL_HEPMC2INCS@

EXTRA_DIST    = @SVNINCLUDE@ $(HEPMC2EXTRADIST)
