set(ShrimpsMain_SOURCES Shrimps.C
                 Cluster_Algorithm.C
                 Soft_Jet_Criterion.C
                 Hadron_Init.C)

set(ShrimpsMain_HEADERS Shrimps.H
                 Cluster_Algorithm.H
                 Soft_Jet_Criterion.H
                 Hadron_Init.H)

sherpa_create_git_info(ShrimpsMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(ShrimpsMain SHARED ${ShrimpsMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ShrimpsMain ShrimpsMain)
add_dependencies(ShrimpsMain ShrimpsMain_Git_Info)
