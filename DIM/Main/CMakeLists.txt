set(DIMMain_SOURCES Gamma.C
                 MCatNLO.C)

set(DIMMain_HEADERS Gamma.H
                 MCatNLO.H)

sherpa_create_git_info(DIMMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(DIMMain SHARED ${DIMMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DIMMain DIMMain)
add_dependencies(DIMMain DIMMain_Git_Info)
