##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
set(GITTAG  CSSHOWER++/Main)

set(Main_SOURCES CS_Cluster_Definitions.C
                 CS_Shower.C)

set(Main_HEADERS CS_Cluster_Definitions.H
                 CS_Shower.H)

sherpa_create_git_info(Main)
 
add_library(CSMain SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(CSMain CSMain)
