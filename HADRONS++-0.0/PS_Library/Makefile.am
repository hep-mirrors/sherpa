if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES = libHadronsPSs.la

libHadronsPSs_la_SOURCES = \
		HD_PS_Base.C \
		Three_Body_PSs.C \
		Two_Body_PSs.C  \
		Four_Body_PSs.C \
		ResonanceFlavour.C

include_HEADERS = \
		HD_PS_Base.H \
		Three_Body_PSs.H \
		Two_Body_PSs.H  \
		Four_Body_PSs.H \
		ResonanceFlavour.H

INCLUDES = \
	-I../../HADRONS++-0.0/Main \
	-I../../HADRONS++-0.0/ME_Library \
	-I../../PHASIC++-1.0/Main \
	-I../../ATOOLS-2.0/Phys \
	-I../../ATOOLS-2.0/Math \
	-I../../ATOOLS-2.0/Org \
	-I../../MODEL-1.0/Main
EXTRA_DIST    = $(CVSINCLUDE)
