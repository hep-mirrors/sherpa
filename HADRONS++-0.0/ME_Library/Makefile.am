if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES = libHadronsMEs.la

libHadronsMEs_la_SOURCES = \
		HD_ME_Base.C \
		HD_ME_Selector.C \
		Three_Body_MEs.C \
		Two_Body_MEs.C 

noinst_HEADERS = \
		HD_ME_Base.H \
		HD_ME_Selector.H \
		Three_Body_MEs.H \
		Two_Body_MEs.H 

INCLUDES = \
	-I../../HADRONS++-0.0/PS_Library \
	-I../../PHASIC++-1.0/Main \
	-I../../ATOOLS-2.0/Phys \
	-I../../ATOOLS-2.0/Math \
	-I../../ATOOLS-2.0/Org 
EXTRA_DIST    = $(CVSINCLUDE)
