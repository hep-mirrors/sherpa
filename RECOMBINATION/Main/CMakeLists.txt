set(RecombinationMain_SOURCES Recombination_Handler.C)

set(RecombinationMain_HEADERS Recombination_Handler.H)

sherpa_create_git_info(RecombinationMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(Recombination SHARED ${RecombinationMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(Recombination Recombination)
add_dependencies(Recombination RecombinationMain_Git_Info)
