set(AmisicMain_SOURCES  Amisic.C)

set(AmisicMain_HEADERS Amisic.H)

sherpa_create_git_info(AmisicMain)
 
add_library(Amisic SHARED ${AmisicMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(Amisic Amisic)
