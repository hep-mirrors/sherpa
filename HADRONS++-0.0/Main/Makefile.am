if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES = libHadronsMain.la

libHadronsMain_la_SOURCES = \
		Hadron_Decays.C \
		Decay_Table_Reader.C \
		Hadron_Decay_Channel.C 
noinst_HEADERS = \
		Hadron_Decays.H \
		Decay_Table_Reader.H \
		Hadron_Decay_Channel.H 

INCLUDES = \
	-I../../HADRONS++-0.0/ME_Library \
	-I../../HADRONS++-0.0/PS_Library \
	-I../../PHASIC++-1.0/Main \
	-I../../ATOOLS-2.0/Phys \
	-I../../ATOOLS-2.0/Math \
	-I../../ATOOLS-2.0/Org 
EXTRA_DIST    = $(CVSINCLUDE)
