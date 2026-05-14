set(SQED_SOURCES  Model.C Form_Factor.C)

set(SQED_HEADERS Model.H)

sherpa_create_git_info(SQED ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaSQED SHARED ${SQED_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaSQED SherpaSQED)
