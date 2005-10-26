SHERPASYS     = ../..

SHERPAFLAGS   =	-pedantic -Wall 

ATOOLSDIR     = $(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    = -I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org 
PDFINCS       =	-I$(SHERPASYS)/MODEL-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/PDF-1.0/Main
APACICINCS    =	-I../Showers

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES = libApacicMain.la
libApacicMain_la_SOURCES      = \
				Knot.C \
				Knot_List.C \
				Tree.C \
				Apacic.C
include_HEADERS       =	\
			Knot.H \
			Knot_List.H \
			Tree.H \
			Apacic.H
AM_CXXFLAGS   = $(SHERPAFLAGS)
INCLUDES = $(APACICINCS) $(PDFINCS) $(MODELINCS) $(ATOOLSINCS)
EXTRA_DIST    = $(CVSINCLUDE)
