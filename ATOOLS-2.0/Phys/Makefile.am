lib_LTLIBRARIES = libToolsPhys.la

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

libToolsPhys_la_SOURCES = Flavour.C Particle.C Decay_Table.C \
                          Flavour_List.C Particle_List.C \
			  Cut_Data.C Selector.C  Combined_Selector.C \
			  Jet_Finder.C Cone_Finder.C Standard_Selector.C \
			  Kt_Algorithm.C Blob.C Blob_List.C Flow.C  
noinst_HEADERS = Flavour.H Particle.H Decay_Table.H \
		 Flavour_List.H Particle_List.H \
		 Cut_Data.H Selector.H Jet_Finder.H Cone_Finder.H \
		 Combined_Selector.H Standard_Selector.H \
		 Kt_Algorithm.H Blob.H Blob_List.H Flow.H Event.H
INCLUDES = -I../Org \
	   -I../Math \
	   -I../../MODEL-1.0/Main

EXTRA_DIST    = $(CVSINCLUDE)