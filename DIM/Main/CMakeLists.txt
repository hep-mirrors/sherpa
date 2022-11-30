set(Main_SOURCES Gamma.C
                 MCatNLO.C)

set(Main_HEADERS Gamma.H
                 MCatNLO.H)

sherpa_create_git_info(Main)
 
add_library(DIMMain SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DIMMain DIMMain)
