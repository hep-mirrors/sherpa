pkglib_LTLIBRARIES       =	libExtraXSNLO.la

libExtraXSNLO_la_SOURCES      =	Virtual_ME2_Base.C \
				WhiteHat_Virtual.C \
				DY_QCD_Virtual.C \
				DIS1jet_QCD_Virtual.C \
				ee3jet_QCD_Virtual.C

localinc_HEADERS          =	Virtual_ME2_Base.H

libExtraXSNLO_la_LIBADD = @CONDITIONAL_WHITEHATLIBS@

libExtraXSNLO_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_WHITEHATINCS@
