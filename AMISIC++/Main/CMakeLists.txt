set(AmisicMain_SOURCES  Amisic.C)

set(AmisicMain_HEADERS Amisic.H)

sherpa_create_git_info(AmisicMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(Amisic SHARED ${AmisicMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(Amisic Amisic)
add_dependencies(Amisic AmisicMain_Git_Info)
