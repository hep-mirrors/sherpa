
ROOTOUT_SOURCES = Output_RootNtuple.C

ROOTOUT_ADDS = Output_RootNtuple.H

ROOTIN_SOURCES = RootNtuple_Reader.C

ROOTIN_ADDS = RootNtuple_Reader.H

if ROOT_SUPPORT
pkglib_LTLIBRARIES    =	libSherpaRootOutput.la libSherpaERootOutput.la libSherpaEDRootOutput.la \
		libSherpaRootInput.la libSherpaERootInput.la libSherpaEDRootInput.la

libSherpaRootOutput_la_SOURCES        =	$(ROOTOUT_SOURCES)
libSherpaERootOutput_la_SOURCES       =	$(ROOTOUT_SOURCES)
libSherpaEDRootOutput_la_SOURCES      =	$(ROOTOUT_SOURCES)

libSherpaRootOutput_la_LIBADD         =	@CONDITIONAL_ROOTLIBS@
libSherpaERootOutput_la_LIBADD        =	@CONDITIONAL_ROOTLIBS@
libSherpaEDRootOutput_la_LIBADD       =	@CONDITIONAL_ROOTLIBS@

libSherpaRootOutput_la_CPPFLAGS       =	$(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@
libSherpaERootOutput_la_CPPFLAGS      =	$(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@
libSherpaEDRootOutput_la_CPPFLAGS     =	$(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@

libSherpaRootOutput_la_CXXFLAGS       = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@
libSherpaERootOutput_la_CXXFLAGS      = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@
libSherpaEDRootOutput_la_CXXFLAGS     = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@

libSherpaRootInput_la_SOURCES         =	$(ROOTIN_SOURCES)
libSherpaERootInput_la_SOURCES        =	$(ROOTIN_SOURCES)
libSherpaEDRootInput_la_SOURCES       =	$(ROOTIN_SOURCES)

libSherpaRootInput_la_LIBADD          =	@CONDITIONAL_ROOTLIBS@
libSherpaERootInput_la_LIBADD         =	@CONDITIONAL_ROOTLIBS@
libSherpaEDRootInput_la_LIBADD        =	@CONDITIONAL_ROOTLIBS@

libSherpaRootInput_la_CPPFLAGS        =	$(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@
libSherpaERootInput_la_CPPFLAGS       =	$(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@
libSherpaEDRootInput_la_CPPFLAGS      =	$(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@

libSherpaRootInput_la_CXXFLAGS       = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@
libSherpaERootInput_la_CXXFLAGS      = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@
libSherpaEDRootInput_la_CXXFLAGS     = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@

ROOT_EXTRA_DIST = $(ROOTOUT_ADDS) $(ROOTIN_ADDS)
else
ROOT_EXTRA_DIST = $(ROOTOUT_SOURCES) $(ROOTOUT_ADDS) $(ROOTIN_SOURCES) $(ROOTIN_ADDS)
endif

EXTRA_DIST = $(ROOT_EXTRA_DIST)
