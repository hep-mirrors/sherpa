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
SHERPAINCS    =	-I$(SHERPASYS)/SHERPA-1.0/PerturbativePhysics
USING__Sherpa = TRUE
endif
if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       = libAmisic.la
libAmisic_la_SOURCES  = \
			MI_Base.C \
			Environment.C \
			Simple_Chain.C \
			Amisic.C 
noinst_HEADERS        = \
			Grid_Handler.H \
			Grid_Handler.IC \
			Grid_Handler.C \
			Grid_Creator_Base.H \
			Grid_Creator_Base.IC \
		       	Grid_Creator_Base.C \
			Grid_Creator.H \
			Grid_Creator.IC \
		       	Grid_Creator.C \
			Fixed_Variable_Channel.H \
			Fixed_Variable_Channel.IC \
			Fixed_Variable_Channel.C \
			MI_Base.H \
			MI_Base.IC \
			Environment.H \
			Simple_Chain.H \
			Simple_Chain.IC \
			Amisic.H \
			Amisic.IC 
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(AMISICFLAGS)
INCLUDES      =	$(SHERPAINCS) $(EXTRAXSINCS) $(PHASICINCS) $(PDFINCS) \
		$(BEAMINCS) $(MODELINCS) $(ATOOLSINCS)
