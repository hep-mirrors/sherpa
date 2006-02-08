SHERPASYS     =	../..

EXTRAXSFLAGS  =	# -pedantic -Wall -Winline
DEFS          =

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
PDFINCS       =	-I$(SHERPASYS)/PDF-1.0/Main
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
EXTRAXSINCS   =	-I../Main

if CVS_INCLUDE
CVSINCLUDE    =	.svn
endif

pkglib_LTLIBRARIES       =	libExtraXS2_2.la
libExtraXS2_2_la_SOURCES      =	\
				XS_QCD.C \
				XS_Drell_Yan.C \
				XS_4F.C \
				QED_Processes.C \
				QCD_Processes.C
pkginclude_HEADERS        =	\
			XS_QCD.H \
			XS_Drell_Yan.H \
			XS_4F.H \
			QED_Processes.H \
			QCD_Processes.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(EXTRAXSFLAGS)
INCLUDES      =	$(EXTRAXSINCS) $(PHASICINCS) $(PDFINCS) \
		$(BEAMINCS) $(MODELINCS) $(ATOOLSINCS)
