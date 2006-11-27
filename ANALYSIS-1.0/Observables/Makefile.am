if ROOT_SUPPORT
EXTRAOBS      =	Dalitz_Observables.C
endif

pkglib_LTLIBRARIES       =	libObservables.la

libObservables_la_SOURCES     =	\
				$(EXTRAOBS) \
				List_Modifier.C \
				Global_Calculator.C \
				HT.C \
				PSM_Observables.C \
				Multiplicity.C \
			    	Statistics_Observable.C \
				One_Particle_Selector.C \
				Two_Particle_Selector.C \
				Four_Particle_Calculator.C \
				One_Particle_Extractor.C \
				All_Particle_Selector.C \
				Universal_Selector.C \
			    	MI_Observables.C \
			    	ISR_Statistics.C \
			    	Shower_Statistics.C \
			    	SV_Statistics.C \
			    	One_Particle_Observables.C \
			    	Two_Particle_Observables.C \
				Two_Particle_Observables_2.C \
                                Three_Particle_Observables.C \
			    	Four_Particle_Observables.C \
			    	Four_Particle_Observables_2.C \
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
				Blob_Selector.C \
				Particle_Ordering.C
pkginclude_HEADERS        =	\
			HT.H \
			PSM_Observables.H \
			Multiplicity.H \
			Statistics_Observable.H \
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

AM_CXXFLAGS   =	@SHERPAFLAGS@ @CONDITIONAL_ROOTFLAGS@

INCLUDES      =	@ANALYSISINCS@ @MODELINCS@ @PDFINCS@ \
		@APACICINCS@ @ATOOLSINCS@ @CONDITIONAL_ROOTINCS@

EXTRA_DIST    = @SVNINCLUDE@
