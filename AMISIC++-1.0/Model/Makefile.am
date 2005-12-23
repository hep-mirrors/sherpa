SHERPASYS     =	../..

AMISICFLAGS   =	# -pedantic -Wall -Winline
DEFS          =

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
PDFDIR        =	$(SHERPASYS)/PDF-1.0
PDFINCS       =	-I$(PDFDIR)/Main -I$(PDFDIR)/Remnant
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
EXTRAXSDIR    =	$(SHERPASYS)/EXTRA_XS-1.0
EXTRAXSINCS   =	-I$(EXTRAXSDIR)/Main -I$(EXTRAXSDIR)/Two2Two
AMISICINCS    = -I../Tools -I../Main

if SHERPA_SUPPORT
SHERPADIR     = $(SHERPASYS)/SHERPA-1.0
SHERPAINCS    =	-I$(SHERPADIR)/PerturbativePhysics
endif

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       = libAmisicModel.la
libAmisicModel_la_SOURCES     = \
			       	Grid_Creator.C \
				Profile_Function_Base.C \
				Profile_Function.C \
				Simple_Chain.C \
				Reggeon_Trajectory.C \
				Simple_String.C 
include_HEADERS        = \
			Grid_Creator.H \
			Profile_Function_Base.H \
			Profile_Function.H \
			Simple_Chain.H \
			Reggeon_Trajectory.H \
			Simple_String.H 
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(AMISICFLAGS)
INCLUDES      =	$(SHERPAINCS) $(AMISICINCS) $(EXTRAXSINCS) $(PHASICINCS) \
		$(PDFINCS) $(BEAMINCS) $(MODELINCS) $(ATOOLSINCS)
