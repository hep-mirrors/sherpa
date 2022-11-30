set(Tools_SOURCES Parton.C
                  Singlet.C)

set(Tools_HEADERS Parton.H
                  Singlet.H)

sherpa_create_git_info(Tools)
 
add_library(MCatNLOTools SHARED ${Tools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(MCatNLOTools MCatNLOTools)
