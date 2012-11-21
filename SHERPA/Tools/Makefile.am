pkglib_LTLIBRARIES       =	libSherpaTools.la

SVNTAG = SHERPA/Tools
include ../../svn.make

libSherpaTools_la_SOURCES     =	\
				SVN_Info.C \
				Output_Base.C \
				Event_Reader_Base.C \
				HepEvt_Interface.C \
				Pythia_HepEvt_Translator.C \
				Output_HepEvt.C \
				Output_LHEF.C \
				Analysis_Interface.C

localinc_HEADERS      =	\
			Definitions.H \
			Output_Base.H \
			Event_Reader_Base.H \
			HepEvt_Interface.H \
			Pythia_HepEvt_Translator.H \
			Output_HepEvt.H \
			Output_LHEF.H \
			Analysis_Interface.H

libSherpaTools_la_CXXFLAGS = $(AM_CXXFLAGS)

libSherpaTools_la_CPPFLAGS = $(AM_CPPFLAGS)

EXTRA_DIST    = @SVNINCLUDE@
