set(Decays_SOURCES Cluster_Decayer.C
                   Cluster_Splitter.C)

set(Decays_HEADERS Cluster_Decayer.H
                   Cluster_Splitter.H)

sherpa_create_git_info(Decays)
 
add_library(AhadicDecays SHARED ${Decays_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(AhadicDecays AhadicDecays)
