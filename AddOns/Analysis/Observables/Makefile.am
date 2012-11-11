if ROOT_SUPPORT
EXTRAOBS      =	Dalitz_Observables.C
else
OBSEXTRADIST = Dalitz_Observables.C
endif

ANAOBS_SOURCES        =	\
			Blobdata.C \
			$(EXTRAOBS) \
			Global_Calculator.C \
			Primitive_Observable_Base.C \
			Normalized_Observable.C \
			One_Variable_Observable.C \
			HT.C \
			TotalMass.C \
			LMass.C \
			EV_C_Observables.C \
			EV_C2_Observables.C \
			EV_C4_Observables.C \
			EV_SC_Observables.C \
			PSM_Observables.C \
			Multiplicity.C \
			MI_Observables.C \
			One_Particle_Observables.C \
			One_Particle_Observables_2.C \
			Conditional_One_Particle_Observables.C \
			Two_Particle_Observables.C \
			Two_Particle_Observables_2.C \
                        Three_Particle_Observables.C \
			Four_Particle_Observables.C \
			Four_Particle_Observables_2.C \
			Six_Particle_Observables.C \
			Soft_Photon_Observables.C \
			Jet_Cone_Distribution.C \
			Jet_Observables.C \
			Jet_NJet_Observables.C \
			Jet_X_Observables.C \
			Universal_Observable.C \
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
			ChargedDEta.C

ANAOBS_ADDS           =	\
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
			Six_Particle_Observables.H \
			Soft_Photon_Observables.H \
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

if USING__Analysis
pkglib_LTLIBRARIES       =	libSherpaObservables.la

libSherpaObservables_la_SOURCES       =	$(ANAOBS_SOURCES)

libSherpaObservables_la_CXXFLAGS = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@

libSherpaObservables_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@

ANAOBS_EXTRA_DIST = $(ANAOBS_ADDS)
else
ANAOBS_EXTRA_DIST = $(ANAOBS_SOURCES) $(ANAOBS_ADDS)
endif

EXTRA_DIST = @SVNINCLUDE@ $(OBSEXTRADIST) $(ANAOBS_EXTRA_DIST)
