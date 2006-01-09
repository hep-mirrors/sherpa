SHERPASYS     =	../..

ANALYSISFLAGS =	# -pedantic -Wall #-Winline
DEFS          =

ANALYSISINCS  =	-I../Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFINCS       =	-I$(SHERPASYS)/PDF-1.0/Main
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
APACICINCS    = -I$(SHERPASYS)/APACIC++-2.0/Main

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if ROOT_SUPPORT
ROOTINCS      =	-I@ROOTINCDIR@
ROOTFLAGS     = -Wno-long-long
EXTRAOBS      =	Dalitz_Observables.C
endif

lib_LTLIBRARIES       =	libObservables.la
libObservables_la_SOURCES     =	\
				$(EXTRAOBS) \
				HT.C \
				PSM_Observables.C \
				Multiplicity.C \
			    	Statistics_Observable.C \
				One_Particle_Selector.C \
				One_Particle_Extractor.C \
				All_Particle_Selector.C \
				Universal_Selector.C \
			    	MI_Observables.C \
			    	ISR_Statistics.C \
			    	Shower_Statistics.C \
			    	SV_Statistics.C \
			    	One_Particle_Observables.C \
			    	Two_Particle_Observables.C \
                                Three_Particle_Observables.C \
			    	Four_Particle_Observables.C \
				Midpoint_Cone.C \
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
				Field_Stuart_Cone_Finder.C \
				Blob_Selector.C 
pkginclude_HEADERS        =	\
			HT.H \
			PSM_Observables.H \
			Multiplicity.H \
			Statistics_Observable.H \
			One_Particle_Selector.H \
			One_Particle_Extractor.H \
			All_Particle_Selector.H \
			Universal_Selector.H \
			MI_Observables.H \
			ISR_Statistics.H \
			Shower_Statistics.H \
			One_Particle_Observables.H \
			Two_Particle_Observables.H \
			Three_Particle_Observables.H \
			Four_Particle_Observables.H \
			Midpoint_Cone.H \
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
		$(APACICINCS) $(ATOOLSINCS) $(ROOTINCS)
