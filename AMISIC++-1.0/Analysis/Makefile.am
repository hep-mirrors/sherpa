SHERPASYS     =	../..

AMISICFLAGS   =	-pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
MODELINC      =	-I$(SHERPASYS)/MODEL-1.0/Main

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES               =	libAmisicAnalysis.la
libAmisicAnalysis_la_SOURCES  = \
				Primordial_KT.C \
				Jet_PT.C \
				Jet_Y.C \
				Jet_Multiplicity.C \
				Jet_Observables.C
noinst_HEADERS        =	\
			Primordial_KT.H \
			Jet_PT.H \
			Jet_Y.H \
			Jet_Multiplicity.H \
			Jet_Observables.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(AMISICFLAGS)
INCLUDES      = $(MODELINC) $(ATOOLSINCS)
