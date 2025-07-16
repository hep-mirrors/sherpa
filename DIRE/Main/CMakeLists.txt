set(DireMain_SOURCES Jet_Criterion.C
                 Dire.C)

sherpa_create_git_info(DireMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(DireMain SHARED ${DireMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DireMain DireMain)
add_dependencies(DireMain DireMain_Git_Info)
