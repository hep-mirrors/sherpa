set(Current_Library_SOURCES
	Scatter_Current_Base.C
	Current_ME.C
	Lepton_Lepton.C
	Nucleon_Nucleon.C
	)

set(Current_Library_HEADERS
	Scatter_Current_Base.H
	Current_ME.H
	Lepton_Lepton.H
	Nucleon_Nucleon.H
	)

sherpa_create_git_info(Current_Library)
 
add_library(NeutrinosCurrent_Library SHARED ${Current_Library_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(NeutrinosCurrent_Library NeutrinosCurrent_Library)
