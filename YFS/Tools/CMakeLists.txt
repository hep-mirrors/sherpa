set(YFSTools_SOURCES Debug.C
                 Dipole.C)

set(YFSTools_HEADERS  Debug.H
                 Dipole.H)


sherpa_create_git_info(YFSTools ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(YFSTools SHARED ${YFSTools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(YFSTools YFSTools)
add_dependencies(YFSTools YFSTools_Git_Info)
