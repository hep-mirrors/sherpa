if ROOT_SUPPORT
EXTRAOBS      =	Dalitz_Observables.C
else
OBSEXTRADIST = Dalitz_Observables.C
endif

pkglib_LTLIBRARIES       =	libObservables.la

libObservables_la_SOURCES     =	\
				$(EXTRAOBS) \
				Global_Calculator.C \
			 	Primitive_Observable_Base.C \
			 	Normalized_Observable.C \
			 	One_Variable_Observable.C \
				HT.C \
				PSM_Observables.C \
				Multiplicity.C \
			    	Statistics_Observable.C \
			    	MI_Observables.C \
			    	One_Particle_Observables.C \
			    	One_Particle_Observables_2.C \
			    	Two_Particle_Observables.C \
				Two_Particle_Observables_2.C \
                                Three_Particle_Observables.C \
			    	Four_Particle_Observables.C \
			    	Four_Particle_Observables_2.C \
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
				General_Observable.C \
				BFKL_Observables.C \
				Atlas_TDR_SL_Top_Fitter.C \
				ChargedDEta.C
pkginclude_HEADERS    =	\
			Primitive_Observable_Base.H \
			Normalized_Observable.H \
			HT.H \
			PSM_Observables.H \
			Multiplicity.H \
			MI_Observables.H \
			One_Particle_Observables.H \
			Two_Particle_Observables.H \
			Three_Particle_Observables.H \
			Four_Particle_Observables.H \
			Jet_Cone_Distribution.H \
			Jet_Observables.H \
			Momentum_Tensor.H \
			CD_Parameters.H \
			Sphericity_Observables.H \
			Event_Shapes_EE.H \
			Shape_Observables_EE.H \
			Jet_Mass_and_Broadening.H \
			Scaled_Observables.H \
			Four_Jet_Angles.H \
			ChargedDEta.H

AM_CXXFLAGS   =	@SHERPAFLAGS@ @CONDITIONAL_ROOTFLAGS@

INCLUDES      =	@ANALYSISINCS@ @MODELINCS@ @ATOOLSINCS@ @CONDITIONAL_ROOTINCS@

EXTRA_DIST    = @SVNINCLUDE@ $(OBSEXTRADIST)
