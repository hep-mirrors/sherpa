set(Tools_SOURCES P2P_Translator.C)

set(Tools_HEADERS P2P_Translator.H)

sherpa_create_git_info(Tools)
 
add_library(AlpacaTools SHARED ${Tools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(AlpacaTools AlpacaTools)
