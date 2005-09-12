SHERPASYS = ../..

ATOOLSFLAGS   =	-pedantic -Wall -Winline

ATOOLSINCS    =	-I../Math -I../Org
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main 

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       =	libToolsPhys.la
libToolsPhys_la_SOURCES       =	\
				Flavour.C \
				Particle.C \
				Particle_Qualifier.C \
			        Simple_Polarisation_Info.C \
				Decay_Table.C \
				Flavour_List.C \
				Particle_List.C \
				Cut_Data.C \
				Selector.C \
				Combined_Selector.C \
				Jet_Finder.C \
				Cone_Finder.C \
				Durham_Algorithm.C \
				Cambridge_Algorithm.C \
				Kt_Algorithm.C \
				Standard_Selector.C \
				Blob.C \
				Blob_List.C \
				Parton_Finder.C \
				Color_Tester.C \
				Momentum_Shifter.C \
				HepEvt_Interface.C \
				Flow.C \
				Color.C
noinst_HEADERS        =	\
			Flavour.H \
			Particle.H \
			Particle_Qualifier.H \
			Simple_Polarisation_Info.H \
			Flavour_List.H \
			Particle_List.H \
			Decay_Table.H \
			Cut_Data.H \
			Selector.H \
			Jet_Finder.H \
			Cone_Finder.H \
			Durham_Algorithm.H \
			Cambridge_Algorithm.H \
			Kt_Algorithm.H \
			Combined_Selector.H \
			Standard_Selector.H \
			Flow.H \
			Blob.H \
			Blob_List.H \
			Parton_Finder.H \
			Color_Tester.H \
			Momentum_Shifter.H \
			HepEvt_Interface.H \
			Event.H \
			Color.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(ATOOLSFLAGS)
INCLUDES      = $(MODELINCS) $(ATOOLSINCS)
