set(Glauber_SOURCES Glauber.C)

set(Glauber_HEADERS Glauber.H)

sherpa_create_git_info(Glauber)
 
add_library(ShrimpsGlauber SHARED ${Glauber_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ShrimpsGlauber ShrimpsGlauber)
