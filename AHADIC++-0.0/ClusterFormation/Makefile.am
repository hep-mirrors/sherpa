if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif
lib_LTLIBRARIES = libClusterFormation.la
libClusterFormation_la_SOURCES = \
				Cluster_Formation_Handler.C \
				Cluster_Former.C \
				Gluon_Decayer.C 
INCLUDES = -I../Main \
	   -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org 
noinst_HEADERS = \
		Cluster_Formation_Handler.H \
		Cluster_Former.H \
		Gluon_Decayer.H 
EXTRA_DIST    = $(CVSINCLUDE)
