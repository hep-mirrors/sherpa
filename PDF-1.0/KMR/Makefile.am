SHERPASYS     =	../..

AMISICFLAGS   =	-pedantic -Wall -Winline

PDFINCS       =	-I../Main -I../Sudakov
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    = -I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
APACICINCS    =	-I$(SHERPASYS)/APACIC++-2.0/Showers

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       = libDUPDF.la
libDUPDF_la_SOURCES   =	\
			LL_Branching.C \
			LL_Single_Sudakov.C \
			LL_Sudakov.C \
			Doubly_Unintegrated_PDF.C
noinst_HEADERS        = \
			LL_Branching.H \
			LL_Single_Sudakov.H \
			LL_Sudakov.H \
			Doubly_Unintegrated_PDF.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(AMISICFLAGS)
INCLUDES      =	$(APACICINCS) $(PDFINCS) $(MODELINCS) $(ATOOLSINCS)
