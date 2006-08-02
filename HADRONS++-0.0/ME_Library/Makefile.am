pkglib_LTLIBRARIES = libHadronsMEs.la

libHadronsMEs_la_SOURCES = \
		HD_ME_Base.C \
		Three_Body_MEs.C \
		Two_Body_MEs.C \
		B_Meson_Decay_MEs.C \
		K_Meson_Decay_MEs.C \
		Tau_Decay_MEs.C  \
		Traces.C \
		XYZFuncs.C

pkginclude_HEADERS = \
		HD_ME_Base.H \
		Three_Body_MEs.H \
		Two_Body_MEs.H  \
		B_Meson_Decay_MEs.H \
		K_Meson_Decay_MEs.H \
		Tau_Decay_MEs.H \
		Traces.H \
		XYZFuncs.H



#AM_CXXFLAGS = @SHERPAFLAGS@

INCLUDES = @HADRONSINCS@ @ATOOLSINCS@ @MODELINCS@
EXTRA_DIST    = @SVNINCLUDE@
