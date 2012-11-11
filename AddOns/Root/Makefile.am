
ROOTOUT_SOURCES = Output_RootNtuple.C

ROOTOUT_ADDS = Output_RootNtuple.H

ROOTIN_SOURCES = RootNtuple_Reader.C

ROOTIN_ADDS = RootNtuple_Reader.H

if ROOT_SUPPORT
pkglib_LTLIBRARIES    =	libSherpaRootOutput.la libSherpaRootInput.la

libSherpaRootOutput_la_SOURCES        =	$(ROOTOUT_SOURCES)

libSherpaRootOutput_la_LIBADD         =	@CONDITIONAL_ROOTLIBS@

libSherpaRootOutput_la_CPPFLAGS       =	$(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@

libSherpaRootInput_la_SOURCES         =	$(ROOTIN_SOURCES)

libSherpaRootInput_la_LIBADD          =	@CONDITIONAL_ROOTLIBS@

libSherpaRootInput_la_CPPFLAGS        =	$(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@

ROOT_EXTRA_DIST = @SVNINCLUDE@ $(ROOTOUT_ADDS) $(ROOTIN_ADDS)
else
ROOT_EXTRA_DIST = $(ROOTOUT_SOURCES) $(ROOTOUT_ADDS) $(ROOTIN_SOURCES) $(ROOTIN_ADDS)
endif

EXTRA_DIST = @SVNINCLUDE@ $(ROOT_EXTRA_DIST)
