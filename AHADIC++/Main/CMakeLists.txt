set(AhadicMain_SOURCES  Ahadic.C)

set(AhadicMain_HEADERS Ahadic.H)

sherpa_create_git_info(AhadicMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(AhadicMain SHARED ${AhadicMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(AhadicMain AhadicMain)
add_dependencies(AhadicMain AhadicMain_Git_Info)
