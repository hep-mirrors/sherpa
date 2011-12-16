
LHOLE_SOURCES  =	LH_OLE_Interface.C \
                        LH_OLE_Communicator.C

LHOLE_ADDS     =	LH_OLE_Communicator.H

if USING__LHOLE
pkglib_LTLIBRARIES    =	libSherpaLHOLE.la

libSherpaLHOLE_la_SOURCES      =	$(LHOLE_SOURCES)

libSherpaLHOLE_la_CPPFLAGS     =	$(AM_CPPFLAGS)

LHOLE_EXTRA_DIST = $(LHOLE_ADDS)
else
LHOLE_EXTRA_DIST = $(LHOLE_SOURCES) $(LHOLE_ADDS)
endif

EXTRA_DIST = @SVNINCLUDE@ $(LHOLE_EXTRA_DIST)
