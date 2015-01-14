pkglib_LTLIBRARIES       =	libSherpaTools.la

SVNTAG = SHERPA/Tools
include ../../svn.make

libSherpaTools_la_SOURCES     =	\
				SVN_Info.C \
				Output_Base.C \
				Event_Reader_Base.C \
				HepEvt_Interface.C \
			    	HepMC2_Interface.C \
				Pythia_HepEvt_Translator.C \
				Output_HepEvt.C \
				Output_LHEF.C \
                                Analysis_Interface.C \
                                Scale_Variations.C

localinc_HEADERS      =	\
			Definitions.H \
			Output_Base.H \
			Event_Reader_Base.H \
			HepEvt_Interface.H \
		 	HepMC2_Interface.H \
			Pythia_HepEvt_Translator.H \
			Output_HepEvt.H \
			Output_LHEF.H \
                        Analysis_Interface.H \
                        Scale_Variations.H

libSherpaTools_la_LIBADD = @CONDITIONAL_HEPMC2LIBS@ @CONDITIONAL_LHAPDFLIBS@

libSherpaTools_la_CXXFLAGS = $(AM_CXXFLAGS)

libSherpaTools_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_GZIPINCS@ @CONDITIONAL_HEPMC2INCS@ @CONDITIONAL_LHAPDFINCS@

