if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES = libAhadicMain.la
libAhadicMain_la_SOURCES = \
			Ahadic.C 
INCLUDES = -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org \
	   -I../Tools \
	   -I../Formation \
	   -I../Decays
include_HEADERS = \
		Ahadic.H 
EXTRA_DIST    = $(CVSINCLUDE)
