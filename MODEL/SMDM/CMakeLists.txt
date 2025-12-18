set(SMDM_SOURCES  Model.C)

sherpa_create_git_info(SMDM)
 
add_library(SherpaSMDM SHARED ${SMDM_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
target_link_libraries(SherpaSMDM ModelMain ToolsPhys ToolsMath ToolsYaml ToolsOrg)
sherpa_handle_shared_library(SherpaSMDM SherpaSMDM)
