##### Automatically or semiautomaticaly generated, do not edit.#########
#SMDM/Makefile.am
 
 set(GITTAG  MODEL/SMDM)
set(SMDM_SOURCES  Model.C)
sherpa_create_git_info(SMDM ../SMDM/)
list(TRANSFORM SMDM_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../SMDM/")
list(APPEND SMDM_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaSMDM SHARED ${SMDM_SOURCES} )
sherpa_handle_shared_library(SherpaSMDM SherpaSMDM)
