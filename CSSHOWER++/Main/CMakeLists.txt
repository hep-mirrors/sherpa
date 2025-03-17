set(CSSMain_SOURCES CS_Cluster_Definitions.C
                 CS_Shower.C)

set(CSSMain_HEADERS CS_Cluster_Definitions.H
                 CS_Shower.H)

sherpa_create_git_info(CSSMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(CSMain SHARED ${CSSMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(CSMain CSMain)
add_dependencies(CSMain CSSMain_Git_Info)
