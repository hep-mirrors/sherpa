if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES = libAhadicDecays.la
libAhadicDecays_la_SOURCES = \
			Cluster_Decay_Handler.C \
			Cluster_Decay_Analysis.C \
			Cluster_Decayer_Base.C \
			Cluster_Part.C \
			Hadron_Part.C \
			Pair_Popper.C

INCLUDES = -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org \
	   -I../../MODEL-1.0/Main \
	   -I../Tools \
	   -I../Main
noinst_HEADERS = \
		Cluster_Decay_Handler.H \
		Cluster_Decay_Analysis.H \
		Cluster_Decayer_Base.H \
		Cluster_Part.H \
		Hadron_Part.H  \
		Pair_Popper.H

EXTRA_DIST    = $(CVSINCLUDE)
