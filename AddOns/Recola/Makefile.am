RECOLA_SOURCES      =	\
			Recola_Virtual.C \
			Recola_Born.C \
			Recola_Interface.C

RECOLA_ADDS         =	\
			Recola_Virtual.H \
			Recola_Born.H \
			Recola_Interface.H  

if RECOLA_SUPPORT
pkglib_LTLIBRARIES    =	libSherpaRecola.la

libSherpaRecola_la_SOURCES  = $(RECOLA_SOURCES)

libSherpaRecola_la_CPPFLAGS =	$(AM_CPPFLAGS) @CONDITIONAL_RECOLAINCS@

RECOLA_EXTRA_DIST = $(RECOLA_ADDS)
else
RECOLA_EXTRA_DIST = $(RECOLA_SOURCES) $(RECOLA_ADDS)
endif

EXTRA_DIST = $(RECOLA_EXTRA_DIST)
