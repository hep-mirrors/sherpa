SHERPASYS = ../..

ATOOLSFLAGS   =	# -pedantic -Wall -Winline

ATOOLSINCS    =	-I../Math -I../Org
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main 

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       =	libToolsPhys.la
libToolsPhys_la_SOURCES       =	\
				Flavour.C \
				Particle.C \
				Decay_Table.C \
				Flavour_List.C \
				Particle_List.C \
				Cut_Data.C \
				Selector.C \
				Combined_Selector.C \
				Jet_Finder.C \
				Cone_Finder.C \
				Durham_Algorithm.C \
				Kt_Algorithm.C \
				Standard_Selector.C \
				Blob.C \
				Blob_List.C \
				Momentum_Shifter.C \
				HepEvt_Interface.C \
				Flow.C  
noinst_HEADERS        =	\
			Flavour.H \
			Particle.H \
			Flavour_List.H \
			Particle_List.H \
			Decay_Table.H \
			Remnant_Info.H \
			Cut_Data.H \
			Selector.H \
			Jet_Finder.H \
			Cone_Finder.H \
			Durham_Algorithm.H \
			Kt_Algorithm.H \
			Combined_Selector.H \
			Standard_Selector.H \
			Flow.H \
			Blob.H \
			Blob_List.H \
			Momentum_Shifter.H \
			HepEvt_Interface.H \
			Event.H  
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(ATOOLSFLAGS)
INCLUDES      = $(MODELINCS) $(ATOOLSINCS)
