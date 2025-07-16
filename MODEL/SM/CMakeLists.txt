set(SM_SOURCES  Model.C)

set(SM_HEADERS Model.H)

sherpa_create_git_info(SM ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaSM SHARED ${SM_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaSM SherpaSM)
add_dependencies(SherpaSM SM_Git_Info)
