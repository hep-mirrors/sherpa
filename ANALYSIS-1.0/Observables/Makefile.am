SHERPASYS     =	../..

ANALYSISFLAGS =	-pedantic -Wall #-Winline

ANALYSISINCS  =	-I../Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFINCS       =	-I$(SHERPASYS)/PDF-1.0/Main
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if ROOT_SUPPORT
ROOTINCS      =	-I$(ROOTSYS)/include
ROOTFLAGS     = -DROOT_SUPPORT -Wno-long-long
endif

lib_LTLIBRARIES       =	libObservables.la
libObservables_la_SOURCES     =	\
				HT.C \
				Multiplicity.C \
			    	Statistics_Observable.C \
				One_Particle_Selector.C \
				One_Particle_Extractor.C \
				All_Particle_Selector.C \
			    	MI_Observables.C \
			    	ISR_Statistics.C \
			    	Shower_Statistics.C \
			    	One_Particle_Observables.C \
			    	Two_Particle_Observables.C \
			    	Four_Particle_Observables.C \
			    	Calorimeter_Cone.C \
			    	Jet_Cone_Distribution.C \
			    	Jet_Observables.C \
				Momentum_Tensor.C \
				CD_Parameters.C \
				Sphericity_Observables.C \
			    	Event_Shapes_EE.C \
			    	Shape_Observables_EE.C \
				Jet_Mass_and_Broadening.C \
			    	Scaled_Observables.C \
				Four_Jet_Angles.C \
				Field_Stuart_Cone_Finder.C
noinst_HEADERS        =	\
			HT.H \
			Multiplicity.H \
			Statistics_Observable.H \
			One_Particle_Selector.H \
			One_Particle_Extractor.H \
			All_Particle_Selector.H \
			MI_Observables.H \
			ISR_Statistics.H \
			Shower_Statistics.H \
			One_Particle_Observables.H \
			Two_Particle_Observables.H \
			Four_Particle_Observables.H \
			Calorimeter_Cone.H \
			Jet_Cone_Distribution.H \
			Jet_Observables.H \
			Momentum_Tensor.H \
			CD_Parameters.H \
			Sphericity_Observables.H \
			Event_Shapes_EE.H \
			Shape_Observables_EE.H \
			Jet_Mass_and_Broadening.H \
			Scaled_Observables.H \
			Four_Jet_Angles.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(ANALYSISFLAGS) $(ROOTFLAGS)
INCLUDES      =	$(ANALYSISINCS) $(MODELINCS) $(PDFINCS) \
		$(ATOOLSINCS) $(ROOTINCS)
