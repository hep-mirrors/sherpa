set(Main_SOURCES  Ahadic.C)

set(Main_HEADERS Ahadic.H)

sherpa_create_git_info(Main)
 
add_library(AhadicMain SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(AhadicMain AhadicMain)
