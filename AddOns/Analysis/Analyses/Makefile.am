
ANAANA_SOURCES        =	\
			Analysis_Base.C \
			BFKL_Analysis.C \
			WPolarization_Analysis.C

ANAANA_ADDS           =	\
			Analysis_Base.H

if USING__Analysis
pkglib_LTLIBRARIES    =	libSherpaAnalyses.la

libSherpaAnalyses_la_SOURCES  =	$(ANAANA_SOURCES)

ANAANA_EXTRA_DIST = $(ANAANA_ADDS)
else
ANAANA_EXTRA_DIST = $(ANAANA_SOURCES) $(ANAANA_ADDS)
endif

EXTRA_DIST = @SVNINCLUDE@ $(ANAANA_EXTRA_DIST)
