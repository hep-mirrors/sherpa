if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

SHERPAFLAGS   =	-pedantic -Wall 
#-Winline 

lib_LTLIBRARIES = libApacicMain.la
libApacicMain_la_SOURCES = Apacic.C Knot.C Knot_List.C Tree.C
AM_CXXFLAGS   = $(SHERPAFLAGS)
INCLUDES = -I../../APACIC++-2.0/Showers \
	   -I../../PDF-1.0/Main \
	   -I../../MODEL-1.0/Main \
	   -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org 
noinst_HEADERS = Apacic.H Knot.H Knot_List.H Tree.H
EXTRA_DIST    = $(CVSINCLUDE)
