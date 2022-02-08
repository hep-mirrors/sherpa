
pkglib_LTLIBRARIES    =	libSherpaPythia.la

if PYTHIA8_SUPPORT
libSherpaPythia_la_SOURCES = \
	Git_Info.C \
	Pythia_Jet_Criterion.C \
	Pythia8_Hadronisation.C

libSherpaPythia_la_LIBADD   =	@CONDITIONAL_PYTHIA8LDADD@
libSherpaPythia_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_PYTHIA8CPPFLAGS@

else
libSherpaPythia_la_SOURCES = \
	Git_Info.C \
	Pythia_Jet_Criterion.C
libSherpaPythia_la_CPPFLAGS = $(AM_CPPFLAGS)

EXTRA_DIST = Pythia8_Hadronisation.C
endif

GITTAG = AddOns/Pythia
include ../../git.make

