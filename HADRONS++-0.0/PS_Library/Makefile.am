if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES = libHadronsPSs.la

libHadronsPSs_la_SOURCES = \
		HD_PS_Base.C \
		Three_Body_PSs.C \
		Two_Body_PSs.C  \
		Four_Body_PSs.C \
		Channel_Generator_Base.C \
		Channel_Generator.C \
		Phase_Space_Generator.C \
		Channel_Generator_Decays.C \
		ResonanceFlavour.C

noinst_HEADERS = \
		HD_PS_Base.H \
		Three_Body_PSs.H \
		Two_Body_PSs.H  \
		Four_Body_PSs.H \
		Channel_Generator_Base.H \
		Channel_Generator.H \
		Phase_Space_Generator.H \
		Channel_Generator_Decays.H \
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
