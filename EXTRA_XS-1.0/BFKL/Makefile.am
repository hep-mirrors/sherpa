SHERPASYS     =	../..

EXTRAXSFLAGS  =	-pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
PDFINCS       =	-I$(SHERPASYS)/PDF-1.0/Main
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
EXTRAXSINCS   =	-I../Main

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       =	libExtraXSBFKL.la
libExtraXSBFKL_la_SOURCES     =	\
				LEV_Base.C \
				LEV_Group.C \
				LO_LEV.C \
				BFKL_Sudakov.C \
				Ladder.C
include_HEADERS       =	\
			LEV_Base.H \
			LEV_Group.H \
			LO_LEV.H \
			BFKL_Sudakov.H \
			Ladder.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(EXTRAXSFLAGS)
INCLUDES      =	$(EXTRAXSINCS) $(PHASICINCS) $(PDFINCS) \
		$(BEAMINCS) $(MODELINCS) $(ATOOLSINCS)
