set(HEFT_SOURCES Effective_Higgs_Coupling.C
                 Model.C
                 SVV_LC.C
                 SVVV_LC.C
                 SVVP_LC.C
                 SVVVV_LC.C)

set(HEFT_HEADERS Effective_Higgs_Coupling.H)

sherpa_create_git_info(HEFT ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaHEFT SHARED ${HEFT_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaHEFT SherpaHEFT)
add_dependencies(SherpaHEFT HEFT_Git_Info)
