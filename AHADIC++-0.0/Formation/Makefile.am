if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES = libAhadicFormation.la
libAhadicFormation_la_SOURCES = \
			Cluster_Formation_Handler.C \
			Gluon_Decayer.C \
			Cluster_Former.C \
			Colour_Reconnections.C
INCLUDES = -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org \
	   -I../Tools \
	   -I../Main
noinst_HEADERS = \
		Cluster_Formation_Handler.H \
		Gluon_Decayer.H \
		Cluster_Former.H \
		Colour_Reconnections.H
EXTRA_DIST    = $(CVSINCLUDE)
