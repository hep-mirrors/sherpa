set(SM_SOURCES  Model.C)

set(SM_HEADERS Model.H)

sherpa_create_git_info(ScalarQED ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaScalarQED SHARED ${SM_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaScalarQED SherpaScalarQED)
