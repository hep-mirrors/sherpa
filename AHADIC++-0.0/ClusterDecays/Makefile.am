if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif
lib_LTLIBRARIES = libClusterDecays.la
libClusterDecays_la_SOURCES = \
		Cluster_Decay_Handler.C \
		Cluster_Decay_Mode.C
INCLUDES = -I../Main \
	   -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org 
noinst_HEADERS = \
		Cluster_Decay_Handler.H \
		Cluster_Decay_Mode.H
EXTRA_DIST    = $(CVSINCLUDE)
