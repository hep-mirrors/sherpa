SHERPASYS     =	../..

AMISICFLAGS   =	# -pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
PDFINCS       =	-I$(SHERPASYS)/PDF-1.0/Main
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
EXTRAXSDIR    =	$(SHERPASYS)/EXTRA_XS-1.0
EXTRAXSINCS   =	-I$(EXTRAXSDIR)/Main -I$(EXTRAXSDIR)/Two2Two

if SHERPA_SUPPORT
SHERPADIR     = $(SHERPASYS)/SHERPA-1.0
SHERPAINCS    =	-I$(SHERPADIR)/PerturbativePhysics -I$(SHERPADIR)/SoftPhysics
endif

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       = libAmisicTools.la
libAmisicTools_la_SOURCES     = \
				Amisic_Histogram.C \
				Fixed_Variable_Channel.C \
				Environment.C \
				Semihard_QCD.C 
include_HEADERS        = \
			Amisic_Histogram.H \
			Fixed_Variable_Channel.H \
			Environment.H \
			Semihard_QCD.H 
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(AMISICFLAGS)
INCLUDES      =	$(SHERPAINCS) $(EXTRAXSINCS) $(PHASICINCS) $(PDFINCS) \
		$(BEAMINCS) $(MODELINCS) $(ATOOLSINCS)
