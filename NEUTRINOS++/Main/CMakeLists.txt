set(NeutrinosMain_SOURCES
	Tests.C
	)

set(NeutrinosMain_HEADERS
	Tests.H
	)

sherpa_create_git_info(Main)
 
add_library(NeutrinosMain SHARED ${NeutrinosMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(NeutrinosMain NeutrinosMain)
