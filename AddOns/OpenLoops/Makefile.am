
OPENLOOPS_SOURCES      =	\
			OpenLoops_Virtual.C \
			OpenLoops_Born.C \
			OpenLoops_Interface.C

OPENLOOPS_ADDS         =	\
			OpenLoops_Virtual.H \
			OpenLoops_Born.H \
			OpenLoops_Interface.H

if OPENLOOPS_SUPPORT
pkglib_LTLIBRARIES    =	libSherpaOpenLoops.la

libSherpaOpenLoops_la_SOURCES  = $(OPENLOOPS_SOURCES)

OPENLOOPS_EXTRA_DIST = $(OPENLOOPS_ADDS)
else
OPENLOOPS_EXTRA_DIST = $(OPENLOOPS_SOURCES) $(OPENLOOPS_ADDS)
endif

EXTRA_DIST = @SVNINCLUDE@ $(OPENLOOPS_EXTRA_DIST)
