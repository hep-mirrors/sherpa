if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif
lib_LTLIBRARIES = libApacicMain.la
libApacicMain_la_SOURCES = Hard_Interface.C Knot.C Knot_List.C Tree.C
INCLUDES = -I../../APACIC++-2.0/Showers \
	   -I../../PDF-1.0/Main \
	   -I../../MODEL-1.0/Main \
	   -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org 
noinst_HEADERS = Hard_Interface.H Knot.H Knot_List.H Tree.H
EXTRA_DIST    = $(CVSINCLUDE)
