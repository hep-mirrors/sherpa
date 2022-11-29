##### Automatically or semiautomaticaly generated, do not edit.#########
#Weights/Makefile.am
 
set(GITTAG  AddOns/Weights)
set(Weights_SOURCES  Output_Weights.C)
sherpa_create_git_info(Weights)
add_library(SherpaWeights SHARED ${Weights_SOURCES} {CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
sherpa_handle_shared_library(SherpaWeights SherpaWeights)
