set(Main_SOURCES Shrimps.C
                 Cluster_Algorithm.C
                 Soft_Jet_Criterion.C
                 Hadron_Init.C)

set(Main_HEADERS Shrimps.H
                 Cluster_Algorithm.H
                 Soft_Jet_Criterion.H
                 Hadron_Init.H)

sherpa_create_git_info(Main)
 
add_library(ShrimpsMain SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ShrimpsMain ShrimpsMain)
