if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES = libAhadicMain.la
libAhadicMain_la_SOURCES = \
			Ahadic.C \
			Cluster.C \
			Hadron_List.C \
			Clusters_2_Hadrons.C
INCLUDES = -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org \
	   -I../ClusterFormation \
	   -I../ClusterDecays
noinst_HEADERS = \
		Ahadic.H \
		Cluster.H \
		Hadron_List.H \
		Clusters_2_Hadrons.H
EXTRA_DIST    = $(CVSINCLUDE)
