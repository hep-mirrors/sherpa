SHERPASYS     =	../..

EXTRAXSFLAGS  =	-pedantic -Wall -Winline

SHERPAINCS    =	-I$(SHERPASYS)/SHERPA-1.0/SoftPhysics
EXTRAXSINCS   =	-I../Two2Two
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
PDFINCS       =	-I$(SHERPASYS)/PDF-1.0/Main
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES = libExtraXS.la
libExtraXS_la_SOURCES =	\
			XS_Base.C \
			Single_XS.C \
			XS_Group.C \
			Simple_XS.C \
			XS_Selector.C
noinst_HEADERS        =	\
			XS_Base.H \
			Single_XS.H \
			XS_Group.H \
			Simple_XS.H \
			XS_Selector.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(EXTRAXSFLAGS)
INCLUDES      =	$(SHERPAINCS) $(EXTRAXSINCS) $(PHASICINCS) $(PDFINCS) \
		$(BEAMINCS) $(MODELINCS) $(ATOOLSINCS)
