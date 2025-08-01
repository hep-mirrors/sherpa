set(AmisicTools_SOURCES MI_Parameters.C
                  MPI_Scale_Setter.C
                  MPI_KFactor_Setter.C
                  Hadronic_XSec_Calculator.C
                  NonPerturbative_XSecs.C
                  Beam_Integrator.C
                  Matter_Overlap.C
                  Interaction_Probability.C
                  Over_Estimator.C
                  Lookup_Tables.C)

set(AmisicTools_HEADERS MI_Parameters.H
                  MPI_Scale_Setter.H
                  MPI_KFactor_Setter.H
                  Hadronic_XSec_Calculator.H
                  NonPerturbative_XSecs.H
                  Beam_Integrator.H
                  Matter_Overlap.H
                  Interaction_Probability.H
                  Over_Estimator.H
                  Lookup_Tables.H)

sherpa_create_git_info(AmisicTools ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(AmisicTools SHARED ${AmisicTools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )

sherpa_handle_shared_library(AmisicTools AmisicTools)
add_dependencies(AmisicTools AmisicTools_Git_Info)
