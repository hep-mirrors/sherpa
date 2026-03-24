set(HadronCurrents_SOURCES Current_Base.C
			   VA_0_PiPi.C
			   VA_0_PiPiPi.C
			   Tools.C)
set(HadronCurrents_HEADERS Current_Base.H
			   VA_0_PiPi.H
			   VA_0_PiPiPi.H
			   Tools.H )


sherpa_create_git_info(HadronCurrents)
 
add_library(METoolsHadronCurrents SHARED ${HadronCurrents_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(METoolsHadronCurrents METoolsHadronCurrents)
