set(MCatNLOMain_SOURCES CS_Cluster_Definitions.C
                 CS_Gamma.C
                 CS_MCatNLO.C)

set(MCatNLOMain_HEADERS CS_Cluster_Definitions.H
                 CS_Gamma.H
                 CS_MCatNLO.H)

sherpa_create_git_info(MCatNLOMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(MCatNLOMain SHARED ${MCatNLOMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(MCatNLOMain MCatNLOMain)
add_dependencies(MCatNLOMain MCatNLOMain_Git_Info)
