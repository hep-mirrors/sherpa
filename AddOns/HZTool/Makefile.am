SVNTAG = AddOns__HZTool
include ../../svn.make

HZTOOL_SOURCES        =	\
			SVN_Info.C \
			hznorm.F \
			hztool_sherpa.F \
			HZTool_Interface.C
HZTOOL_ADDS           = \
			HZTool_Wrapper.H

AM_FFLAGS     =	-fno-automatic @CONDITIONAL_HZTOOLINCS@

if HZTOOL_SUPPORT
if CERNLIB_SUPPORT
pkglib_LTLIBRARIES    =	libHZToolSherpa.la

libHZToolSherpa_la_SOURCES    =	$(HZTOOL_SOURCES)

libHZToolSherpa_la_LIBADD     =	@CONDITIONAL_HZTOOLLIBS@

HZTOOL_EXTRA_DIST    = $(HZTOOL_ADDS)
else
HZTOOL_EXTRA_DIST    = $(HZTOOL_ADDS) $(HZTOOL_SOURCES)
endif
else
HZTOOL_EXTRA_DIST    = $(HZTOOL_ADDS) $(HZTOOL_SOURCES)
endif

EXTRA_DIST = @SVNINCLUDE@ $(HZTOOL_EXTRA_DIST)
