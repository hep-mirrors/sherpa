##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
 set(GITTAG  SHRiMPS/Main)
set(Main_SOURCES  Shrimps.C Cluster_Algorithm.C Soft_Jet_Criterion.C Hadron_Init.C)
 #localincdir = $(pkg#includedir)/SHRiMPS/Main)
set(Main_HEADERS Shrimps.H Cluster_Algorithm.H Soft_Jet_Criterion.H Hadron_Init.H)
sherpa_create_git_info(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ShrimpsMain SHARED ${Main_SOURCES} )
sherpa_handle_shared_library(ShrimpsMain ShrimpsMain)
