 
set(Weights_SOURCES  Output_Weights.C)
sherpa_create_git_info(Weights ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
add_library(SherpaWeights SHARED ${Weights_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
sherpa_handle_shared_library(SherpaWeights SherpaWeights)
add_dependencies(SherpaWeights Weights_Git_Info)
