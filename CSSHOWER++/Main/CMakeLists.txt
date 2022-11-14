##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
 set(GITTAG  CSSHOWER++/Main)
set(Main_SOURCES  CS_Cluster_Definitions.C CS_Shower.C)
set(Main_HEADERS CS_Cluster_Definitions.H CS_Shower.H)
sherpa_create_git_info(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(CSMain SHARED ${Main_SOURCES} )
sherpa_handle_shared_library(CSMain CSMain)
