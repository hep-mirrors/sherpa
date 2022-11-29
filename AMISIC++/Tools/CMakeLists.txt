##### Automatically or semiautomaticaly generated, do not edit.#########
#Tools/Makefile.am
 
set(GITTAG  AMISIC++/Tools)

set(Tools_SOURCES MI_Parameters.C
                  MPI_Scale_Setter.C
                  MPI_KFactor_Setter.C
                  Hadronic_XSec_Calculator.C
                  Impact_Parameter.C
                  Matter_Overlap.C
                  Interaction_Probability.C
                  Over_Estimator.C)

set(Tools_HEADERS MI_Parameters.H
                  MPI_Scale_Setter.H
                  MPI_KFactor_Setter.H
                  Hadronic_XSec_Calculator.H
                  Impact_Parameter.H
                  Matter_Overlap.H
                  Interaction_Probability.H
                  Over_Estimator.H)

sherpa_create_git_info(Tools ../Tools/)
list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
list(APPEND Tools_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(AmisicTools SHARED ${Tools_SOURCES} )
sherpa_handle_shared_library(AmisicTools AmisicTools)
