SHERPASYS     =	../..

AMISICFLAGS   =	-pedantic -Wall -Winline

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
SHERPAFLAGS   =	-DUSING__Sherpa 
else
endif

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       = libAmisic.la
libAmisic_la_SOURCES  = \
			MI_Base.C \
			Environment.C \
			Amisic_Histogram.C \
		       	Grid_Creator.C \
			Profile_Function_Base.C \
			Profile_Function.C \
			Fixed_Variable_Channel.C \
			Semihard_QCD.C \
			Simple_Chain.C \
			Amisic.C 
noinst_HEADERS        = \
			Amisic_Histogram.H \
			Grid_Handler.H \
			Grid_Handler.C \
			Grid_Creator.H \
			MI_Base.H \
			Environment.H \
			Profile_Function_Base.H \
			Profile_Function.H \
			Fixed_Variable_Channel.H \
			Semihard_QCD.H \
			Simple_Chain.H \
			Amisic.H 
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(AMISICFLAGS) $(SHERPAFLAGS)
INCLUDES      =	$(SHERPAINCS) $(EXTRAXSINCS) $(PHASICINCS) $(PDFINCS) \
		$(BEAMINCS) $(MODELINCS) $(ATOOLSINCS)
