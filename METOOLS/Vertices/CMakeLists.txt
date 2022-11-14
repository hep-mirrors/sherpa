##### Automatically or semiautomaticaly generated, do not edit.#########
#Vertices/Makefile.am
 
 set(GITTAG  METOOLS/Vertices)
set(Vertices_SOURCES  FFV_DLC.C VVV_DLC.C SSS_LC.C SSV_LC.C SFF_LC.C SVV_LC.C FFV_LC.C VVV_LC.C VVP_LC.C SSSS_LC.C SSVV_LC.C VVVV_LC.C)
sherpa_create_git_info(Vertices ../Vertices/)
list(TRANSFORM Vertices_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Vertices/")
list(APPEND Vertices_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(METoolsVertices SHARED ${Vertices_SOURCES} )
sherpa_handle_shared_library(METoolsVertices METoolsVertices)
