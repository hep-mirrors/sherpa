SHERPASYS     =	../..

EXTRAXSFLAGS  =	# -pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
PDFINCS       =	-I$(SHERPASYS)/PDF-1.0/Main
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
EXTRAXSINCS   =	-I../Two2Two

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES = libExtraXS.la
libExtraXS_la_SOURCES =	\
			XS_Base.C \
			Single_XS.C \
			XS_Group.C \
			SimpleXSecs.C \
			XS_Selector.C
noinst_HEADERS        =	\
			XS_Base.H \
			Single_XS.H \
			XS_Group.H \
			SimpleXSecs.H \
			XS_Selector.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(EXTRAXSFLAGS)
INCLUDES      =	$(EXTRAXSINCS) $(PHASICINCS) $(PDFINCS) \
		$(BEAMINCS) $(MODELINCS) $(ATOOLSINCS)
