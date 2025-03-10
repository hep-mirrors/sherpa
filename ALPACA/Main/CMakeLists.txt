set(Main_SOURCES Alpaca.C)

set(Main_HEADERS Alpaca.H)

sherpa_create_git_info(Main)
 
add_library(AlpacaMain SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(AlpacaMain AlpacaMain)
