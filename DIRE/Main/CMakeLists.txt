set(DireMain_SOURCES Jet_Criterion.C
                 Dire.C)

sherpa_create_git_info(DireMain)
 
add_library(DireMain SHARED ${DireMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DireMain DireMain)
