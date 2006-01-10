SHERPASYS     =	../..

PDFFLAGS      =	# -pedantic -Wall -Winline
DEFS          =

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
APACICINCS    =	-I$(SHERPASYS)/APACIC++-2.0/Showers
PDFINCS       =	-I../LHAPDF -I../Main

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

pkglib_LTLIBRARIES       =	libSudakov.la
libSudakov_la_SOURCES =	\
			NLL_Branching_Probability_Base.C \
			NLL_Branching_Probabilities.C \
			NLL_Sudakov_Base.C \
			NLL_Single_Sudakov.C \
			NLL_Combined_Sudakov.C \
			NLL_Sudakov.C \
			NLL_JetRate.C
pkginclude_HEADERS        =	\
			NLL_Branching_Probability_Base.H \
			NLL_Branching_Probabilities.H \
			NLL_Sudakov_Base.H \
			NLL_Single_Sudakov.H \
			NLL_Combined_Sudakov.H \
			NLL_Sudakov.H \
			NLL_JetRate.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(PDFFLAGS)
INCLUDES      =	$(APACICINCS) $(PDFINCS) $(MODELINCS) $(ATOOLSINCS)

