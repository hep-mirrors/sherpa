
HZTOOL_SOURCES        =	\
			hznorm.F \
			hztool_sherpa.F \
			HZTool_Interface.C
HZTOOL_ADDS           = \
			HZTool_Wrapper.H

AM_FFLAGS     =	-fno-automatic @CONDITIONAL_HZTOOLINCS@

if HZTOOL_SUPPORT
if CERNLIB_SUPPORT
pkglib_LTLIBRARIES    =	libSherpaHZToolAnalysis.la

libSherpaHZToolAnalysis_la_SOURCES    =	$(HZTOOL_SOURCES)

libSherpaHZToolAnalysis_la_LIBADD     =	@CONDITIONAL_HZTOOLLIBS@ @CONDITIONAL_CERNLIBLIBS@

HZTOOL_EXTRA_DIST    = $(HZTOOL_ADDS)
else
HZTOOL_EXTRA_DIST    = $(HZTOOL_ADDS) $(HZTOOL_SOURCES)
endif
else
HZTOOL_EXTRA_DIST    = $(HZTOOL_ADDS) $(HZTOOL_SOURCES)
endif

EXTRA_DIST = $(HZTOOL_EXTRA_DIST)
