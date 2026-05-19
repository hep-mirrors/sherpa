set(SherpaMain_SOURCES Sherpa.C
                 Filter.C)
set(SherpaMain_HEADERS Sherpa.H
                 Filter.H)
sherpa_create_git_info(SherpaMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaMain SHARED ${SherpaMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaMain SherpaMain)
add_dependencies(SherpaMain SherpaMain_Git_Info)
