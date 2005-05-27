if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES = libAhadicTools.la
libAhadicTools_la_SOURCES = \
			Cluster.C \
			Constituents.C \
			Hadron_Wave_Function.C \
			Hadron_Multiplet.C \
			Transitions.C \
			Cluster_Transformer.C \
			Hadronisation_Parameters.C
INCLUDES = -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org 
noinst_HEADERS = \
		Cluster.H \
		Constituents.H \
		Hadron_Wave_Function.H \
		Hadron_Multiplet.H \
		Transitions.H \
		Cluster_Transformer.H \
		Hadronisation_Parameters.H
EXTRA_DIST    = $(CVSINCLUDE)
