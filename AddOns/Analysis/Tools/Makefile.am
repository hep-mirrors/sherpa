GITTAG = AddOns/Analysis/Tools
include ../../../git.make

ANATOOLS_SOURCES      =	\
			Git_Info.C \
			Lorentz.C \
			Breit_Frame.C \
			MET.C \
			Blob_Selector.C \
			List_Creator.C \
			List_Merger.C \
			Particle_Ordering.C \
			Particle_Qualifier.C \
			Particle_Selector.C 

ANATOOLS_ADDS         =	\
			Particle_Qualifier.H 

if USING__Analysis
pkglib_LTLIBRARIES    =	libSherpaAnalysisTools.la

libSherpaAnalysisTools_la_SOURCES     = $(ANATOOLS_SOURCES)

localinc_HEADERS      = $(ANATOOLS_ADDS)

libSherpaAnalysisTools_la_CXXFLAGS = $(AM_CXXFLAGS)

libSherpaAnalysisTools_la_CPPFLAGS = $(AM_CPPFLAGS)

ANATOOLS_EXTRA_DIST   =	$(ANATOOLS_ADDS)
else
ANATOOLS_EXTRA_DIST   =	$(ANATOOLS_SOURCES) $(ANATOOLS_ADDS)
endif

EXTRA_DIST = $(ANATOOLS_EXTRA_DIST)
