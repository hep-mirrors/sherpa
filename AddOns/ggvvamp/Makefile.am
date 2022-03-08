GITTAG = AddOns/ggvvamp
include ../../git.make

GGVVAMP_SOURCES=Virtual_Amp.C \
				Precise_Amp.C

GGVVAMP_ADDS=Virtual_Amp.H \
			 Precise_Amp.H

if GGVVAMP_SUPPORT
pkglib_LTLIBRARIES = libSherpaggvvamp.la
libSherpaggvvamp_la_SOURCES = $(GGVVAMP_SOURCES)
libSherpaggvvamp_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_GGVVAMPINCS@
libSherpaggvvamp_la_LIBADD   =	@CONDITIONAL_GGVVAMPLIBS@
GGVVAMP_EXTRA_DIST = $(GGVVAMP_ADDS)
else
GGVVAMP_EXTRA_DIST = $(GGVVAMP_SOURCES) $(GGVVAMP_ADDS)
endif

EXTRA_DIST = $(GGVVAMP_EXTRA_DIST)
