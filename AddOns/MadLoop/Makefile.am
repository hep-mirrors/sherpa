MADLOOP_SOURCES      =	\
			ML5_Interface.C

if MADLOOP_SUPPORT
pkglib_LTLIBRARIES    =	libSherpaML5.la

libSherpaML5_la_SOURCES  = $(MADLOOP_SOURCES)

libSherpaML5_la_CPPFLAGS = $(AM_CPPFLAGS)
else
MADLOOP_EXTRA_DIST = $(MADLOOP_SOURCES)
endif

EXTRA_DIST = $(MADLOOP_EXTRA_DIST)

nodist_bin_SCRIPTS = sconsloops
