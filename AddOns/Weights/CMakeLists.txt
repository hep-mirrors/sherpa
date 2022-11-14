##### Automatically or semiautomaticaly generated, do not edit.#########
#Weights/Makefile.am
 
 set(GITTAG  AddOns/Weights)
set(Weights_SOURCES  Output_Weights.C)
sherpa_create_git_info(Weights ../Weights/)
list(TRANSFORM Weights_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Weights/")
list(APPEND Weights_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaWeights SHARED ${Weights_SOURCES} )
sherpa_handle_shared_library(SherpaWeights SherpaWeights)
