lib_LTLIBRARIES = libToolsPhys.la

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

libToolsPhys_la_SOURCES = Flavour.C Particle.C Decay_Table.C \
                          Flavour_List.C Particle_List.C \
			  Cut_Data.C Selector.C  Combined_Selector.C \
			  Jet_Finder.C Cone_Finder.C Standard_Selector.C \
			  Primitive_Analysis.C Primitive_Observable.C Shower_Observables.C \
			  Blob.C Blob_List.C Flow.C  
noinst_HEADERS = Flavour.H Particle.H Flavour_List.H Particle_List.H Decay_Table.H \
		 Cut_Data.H Selector.H Jet_Finder.H Cone_Finder.H \
		 Combined_Selector.H Standard_Selector.H \
		 Primitive_Analysis.H Primitive_Observable.H Primitive_Observable_Base.H \
                 Shower_Observables.H \
		 Flow.H Blob.H Blob_List.H Event.H  
INCLUDES = -I../Org \
	   -I../Math \
	   -I../../MODEL-1.0/Main

EXTRA_DIST    = $(CVSINCLUDE)