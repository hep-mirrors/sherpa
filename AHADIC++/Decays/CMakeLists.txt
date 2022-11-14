##### Automatically or semiautomaticaly generated, do not edit.#########
#Decays/Makefile.am
 
 set(GITTAG  AHADIC++/Decays)
set(Decays_SOURCES  Cluster_Decayer.C Cluster_Splitter.C)
 #localincdir = $(pkg#includedir)/AHADIC++/Decays)
set(Decays_HEADERS Cluster_Decayer.H Cluster_Splitter.H)
sherpa_create_git_info(Decays ../Decays/)
list(TRANSFORM Decays_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Decays/")
list(APPEND Decays_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(AhadicDecays SHARED ${Decays_SOURCES} )
sherpa_handle_shared_library(AhadicDecays AhadicDecays)
