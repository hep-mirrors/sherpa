set(MCatNLOTools_SOURCES Parton.C
                  Singlet.C)

set(MCatNLOTools_HEADERS Parton.H
                  Singlet.H)

sherpa_create_git_info(MCatNLOTools ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(MCatNLOTools SHARED ${MCatNLOTools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(MCatNLOTools MCatNLOTools)
