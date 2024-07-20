 
set(Initialization_SOURCES  Initialization_Handler.C)
set(Initialization_HEADERS Initialization_Handler.H)
sherpa_create_git_info(Initialization)
 
add_library(SherpaInitialization SHARED ${Initialization_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
target_link_libraries(SherpaInitialization PRIVATE YFSMain)
sherpa_handle_shared_library(SherpaInitialization SherpaInitialization)
