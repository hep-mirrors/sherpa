set(Observables_SOURCES Event_Output.C
                        Blobdata.C
                        Global_Calculator.C
                        Primitive_Observable_Base.C
                        Normalized_Observable.C
                        One_Variable_Observable.C
                        HT.C
                        XS.C
                        TotalMass.C
                        LMass.C
                        EV_C_Observables.C
                        EV_C2_Observables.C
                        EV_C4_Observables.C
                        EV_SC_Observables.C
                        PSM_Observables.C
                        Multiplicity.C
                        One_Particle_Observables.C
                        One_Particle_Observables_2.C
                        Conditional_One_Particle_Observables.C
                        Two_Particle_Observables.C
                        Two_Particle_Observables_2.C
                        Three_Particle_Observables.C
                        Four_Particle_Observables.C
                        Four_Particle_Observables_2.C
                        Six_Particle_Observables.C
                        Soft_Photon_Observables.C
                        Jet_Cone_Distribution.C
                        Jet_Observables.C
                        Jet_NJet_Observables.C
                        Jet_X_Observables.C
                        Momentum_Tensor.C
                        CD_Parameters.C
                        Sphericity_Observables.C
                        Event_Shapes_EE.C
                        Shape_Observables_EE.C
                        Jet_Mass_and_Broadening.C
                        Scaled_Observables.C
                        Four_Jet_Angles.C
                        ChargedDEta.C)

set(Observables_HEADERS Primitive_Observable_Base.H
                        Normalized_Observable.H
                        HT.H
                        PSM_Observables.H
                        Multiplicity.H
                        MI_Observables.H
                        One_Particle_Observables.H
                        Two_Particle_Observables.H
                        Three_Particle_Observables.H
                        Four_Particle_Observables.H
                        Six_Particle_Observables.H
                        Soft_Photon_Observables.H
                        Jet_Cone_Distribution.H
                        Jet_Observables.H
                        Momentum_Tensor.H
                        CD_Parameters.H
                        Sphericity_Observables.H
                        Event_Shapes_EE.H
                        Shape_Observables_EE.H
                        Jet_Mass_and_Broadening.H
                        Scaled_Observables.H
                        Four_Jet_Angles.H
                        ChargedDEta.H)

sherpa_create_git_info(Observables ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaObservables SHARED ${Observables_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaObservables SherpaObservables)
