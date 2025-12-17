set(SM_SOURCES  Model.C)

set(SM_HEADERS Model.H)

sherpa_create_git_info(SM)
 
add_library(SherpaSM SHARED ${SM_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
target_link_libraries(SherpaSM ModelMain ToolsPhys ToolsMath ToolsYaml ToolsOrg)
sherpa_handle_shared_library(SherpaSM SherpaSM)
