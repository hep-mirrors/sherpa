
GITTAG = AddOns/MCFM
include ../../git.make

MCFM_SOURCES  =	\
		MCFM_Interface.C

if MCFM_SUPPORT
pkglib_LTLIBRARIES    =	libSherpaMCFM.la

libSherpaMCFM_la_SOURCES      =	$(MCFM_SOURCES)

libSherpaMCFM_la_LIBADD       =	@CONDITIONAL_MCFMLIBS@

libSherpaMCFM_la_CPPFLAGS     =	$(AM_CPPFLAGS) -fopenmp @CONDITIONAL_MCFMINCS@

MCFM_EXTRA_DIST =
else
MCFM_EXTRA_DIST = $(MCFM_SOURCES)
endif

EXTRA_DIST = $(MCFM_EXTRA_DIST)
