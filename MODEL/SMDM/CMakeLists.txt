set(SMDM_SOURCES  Model.C)

sherpa_create_git_info(SMDM)
 
add_library(SherpaSMDM SHARED ${SMDM_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaSMDM SherpaSMDM)
