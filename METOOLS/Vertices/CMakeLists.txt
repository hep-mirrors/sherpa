set(Vertices_SOURCES FFV_DLC.C
                     VVV_DLC.C
                     SSS_LC.C
                     SSV_LC.C
                     SFF_LC.C
                     SVV_LC.C
                     FFV_LC.C
                     VVV_LC.C
                     VVP_LC.C
                     SSSS_LC.C
                     SSVV_LC.C
                     VVVV_LC.C)

sherpa_create_git_info(Vertices ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(METoolsVertices SHARED ${Vertices_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(METoolsVertices METoolsVertices)
add_dependencies(METoolsVertices Vertices_Git_Info)
