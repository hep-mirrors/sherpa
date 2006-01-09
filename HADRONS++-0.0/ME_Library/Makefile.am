if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

DEFS          =

lib_LTLIBRARIES = libHadronsMEs.la

libHadronsMEs_la_SOURCES = \
		HD_ME_Base.C \
		HD_ME_Selector.C \
		Three_Body_MEs.C \
		Two_Body_MEs.C \
		B_Meson_Decay_MEs.C \
		K_Meson_Decay_MEs.C \
		Tau_Decay_MEs.C  \
		Traces.C \
		XYZFuncs.C \
		Tools.C

pkginclude_HEADERS = \
		HD_ME_Base.H \
		HD_ME_Selector.H \
		Three_Body_MEs.H \
		Two_Body_MEs.H  \
		B_Meson_Decay_MEs.H \
		K_Meson_Decay_MEs.H \
		Tau_Decay_MEs.H \
		Traces.H \
		XYZFuncs.H \
		Tools.H

INCLUDES = \
	-I../../HADRONS++-0.0/PS_Library \
	-I../../PHASIC++-1.0/Main \
	-I../../ATOOLS-2.0/Phys \
	-I../../ATOOLS-2.0/Math \
	-I../../ATOOLS-2.0/Org \
    -I../../MODEL-1.0/Main
EXTRA_DIST    = $(CVSINCLUDE)
