set(CEEX_SOURCES Ceex_Base.C)

set(CEEX_HEADERS Ceex_Base.H)


sherpa_create_git_info(CEEX ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(YFSCEEX SHARED ${CEEX_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(YFSCEEX YFSCEEX)
add_dependencies(YFSCEEX CEEX_Git_Info)
