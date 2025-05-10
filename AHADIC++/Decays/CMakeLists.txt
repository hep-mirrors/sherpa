set(AhadicDecays_SOURCES Cluster_Decayer.C
                   Cluster_Splitter.C)

set(AhadicDecays_HEADERS Cluster_Decayer.H
                   Cluster_Splitter.H)

sherpa_create_git_info(AhadicDecays ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(AhadicDecays SHARED ${AhadicDecays_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(AhadicDecays AhadicDecays)
add_dependencies(AhadicDecays AhadicDecays_Git_Info)
