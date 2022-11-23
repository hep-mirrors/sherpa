##### Automatically or semiautomaticaly generated, do not edit.#########
#UFO/Makefile.am
 
set(GITTAG  MODEL/UFO)
set(UFO_SOURCES  UFO_Param_Reader.C UFO_Model.C UFO_Color_Functions.C)
set(UFO_HEADERS UFO_Param_Reader.H UFO_Model.H UFO_Color_Functions.H)
sherpa_create_git_info(UFO ../UFO/)
list(TRANSFORM UFO_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../UFO/")
list(APPEND UFO_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ModelUFO SHARED ${UFO_SOURCES} )
sherpa_handle_shared_library(ModelUFO ModelUFO)
