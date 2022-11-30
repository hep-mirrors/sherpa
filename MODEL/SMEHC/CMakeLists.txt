set(SMEHC_SOURCES Effective_Higgs_Coupling.C
                  Model.C
                  SVV_LC.C
                  SVVV_LC.C
                  SVVP_LC.C
                  SVVVV_LC.C)

set(SMEHC_HEADERS Effective_Higgs_Coupling.H)

sherpa_create_git_info(SMEHC)
 
add_library(SherpaSMEHC SHARED ${SMEHC_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaSMEHC SherpaSMEHC)
