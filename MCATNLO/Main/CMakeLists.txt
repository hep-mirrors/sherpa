##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
set(GITTAG  MCATNLO/Main)

set(Main_SOURCES CS_Cluster_Definitions.C
                 CS_Gamma.C
                 CS_MCatNLO.C)

set(Main_HEADERS CS_Cluster_Definitions.H
                 CS_Gamma.H
                 CS_MCatNLO.H)

sherpa_create_git_info(Main)
 
add_library(MCatNLOMain SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(MCatNLOMain MCatNLOMain)
