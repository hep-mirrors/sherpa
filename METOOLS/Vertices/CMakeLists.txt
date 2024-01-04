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
                     VVVV_LC.C
                     RRV_LC.C
                     RFS_LC.C)

sherpa_create_git_info(Vertices)
 
add_library(METoolsVertices SHARED ${Vertices_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(METoolsVertices METoolsVertices)
