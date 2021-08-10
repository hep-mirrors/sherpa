RECOLA_SOURCES= \
				Recola_Interface.C \
				Recola_Born.C \
				Recola_Virtual.C 

RECOLA_ADDS= 	\
				Recola_Interface.H \ 
				Recola_Born.H \
				Recola_Virtual.H 

if RECOLA_SUPPORT
pkglib_LTLIBRARIES = libSherpaRecola.la

libSherpaRecola_la_SOURCES  = $(RECOLA_SOURCES)

libSherpaRecola_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_RECOLAINCS@

RECOLA_EXTRA_DIST = $(RECOLA_ADDS)
else
	RECOLA_EXTRA_DIST = $(RECOLA_SOURCES) $(RECOLA_ADDS)
endif

EXTRA_DIST = $(RECOLA_EXTRA_DIST)

