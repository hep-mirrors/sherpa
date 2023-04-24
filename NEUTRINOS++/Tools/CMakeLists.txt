set(NeutrinosTools_SOURCES
	Form_Factor_Library.C
	Form_Factor_Parameter_Maps.C
	)

set(NeutrinosTools_HEADERS
	Form_Factor_Library.H
	Form_Factor_Parameter_Maps.H
	)

sherpa_create_git_info(Tools)
 
add_library(NeutrinosTools SHARED ${NeutrinosTools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(NeutrinosTools NeutrinosTools)
