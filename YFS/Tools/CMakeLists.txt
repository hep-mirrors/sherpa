set(Tools_SOURCES Debug.C
                 Dipole.C)

set(Tools_HEADERS  Debug.H
                 Dipole.H)


sherpa_create_git_info(Tools)
 
add_library(YFSTools SHARED ${Tools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(YFSTools YFSTools)
