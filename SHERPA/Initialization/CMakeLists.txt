 
set(Initialization_SOURCES  Initialization_Handler.C)
set(Initialization_HEADERS Initialization_Handler.H)
sherpa_create_git_info(Initialization ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaInitialization SHARED ${Initialization_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
target_link_libraries(SherpaInitialization PRIVATE YFSMain)
sherpa_handle_shared_library(SherpaInitialization SherpaInitialization)
add_dependencies(SherpaInitialization Initialization_Git_Info)
