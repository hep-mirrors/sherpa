GITTAG = AddOns/Analysis/Analyses
include ../../../git.make

ANAANA_SOURCES        =	\
			Git_Info.C \
			Analysis_Base.C \
			WPolarization_Analysis.C \
			ATLAS_CONF_2011_038.C \
			xQ2_Analysis.C

ANAANA_ADDS           =	\
			Analysis_Base.H

if USING__Analysis
pkglib_LTLIBRARIES    =	libSherpaAnalyses.la

libSherpaAnalyses_la_SOURCES  =	$(ANAANA_SOURCES)

localinc_HEADERS      = $(ANAANA_ADDS)

ANAANA_EXTRA_DIST = $(ANAANA_ADDS)
else
ANAANA_EXTRA_DIST = $(ANAANA_SOURCES) $(ANAANA_ADDS)
endif

EXTRA_DIST = $(ANAANA_EXTRA_DIST)
