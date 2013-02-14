ANATOOLS_SOURCES      =	\
			Lorentz.C \
			Breit_Frame.C \
			MET.C \
			Blob_Selector.C \
			List_Creator.C \
			List_Merger.C \
			List_Modifier.C \
			Particle_Ordering.C \
			Particle_Qualifier.C \
			Universal_Selector.C \
			Particle_Selector.C 

ANATOOLS_ADDS         =	\
			Particle_Qualifier.H 

if USING__Analysis
pkglib_LTLIBRARIES       =	libSherpaAnalysisTools.la

libSherpaAnalysisTools_la_SOURCES     = $(ANATOOLS_SOURCES)

libSherpaAnalysisTools_la_CXXFLAGS = $(AM_CXXFLAGS)

libSherpaAnalysisTools_la_CPPFLAGS = $(AM_CPPFLAGS)

ANATOOLS_EXTRA_DIST   =	$(ANATOOLS_ADDS)
else
ANATOOLS_EXTRA_DIST   =	$(ANATOOLS_SOURCES) $(ANATOOLS_ADDS)
endif

EXTRA_DIST = @SVNINCLUDE@ $(ANATOOLS_EXTRA_DIST)
