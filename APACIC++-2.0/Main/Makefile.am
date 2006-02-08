SHERPASYS     = ../..

SHERPAFLAGS   =	# -pedantic -Wall 
DEFS          =

ATOOLSDIR     = $(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    = -I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org 
PDFINCS       =	-I$(SHERPASYS)/MODEL-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/PDF-1.0/Main
APACICINCS    =	-I../Showers

if CVS_INCLUDE
CVSINCLUDE    =	.svn
endif

pkglib_LTLIBRARIES = libApacicMain.la
libApacicMain_la_SOURCES      = \
				Knot.C \
				Tree.C \
				Apacic.C
pkginclude_HEADERS       =	\
			Knot.H \
			Tree.H \
			Veto_Info.H \
			Apacic.H
AM_CXXFLAGS   = $(SHERPAFLAGS)
INCLUDES = $(APACICINCS) $(PDFINCS) $(MODELINCS) $(ATOOLSINCS)
EXTRA_DIST    = $(CVSINCLUDE)
