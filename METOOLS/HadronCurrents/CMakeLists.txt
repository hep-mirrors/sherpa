set(HadronCurrents_SOURCES Current_Base.C
			   VA_0_PiPi.C
			   VA_0_PiPiPi.C
			   VA_0_KKPi.C
			   VA_0_KPiK.C
			   VA_0_KPiPi.C
			   Tools.C )
set(HadronCurrents_HEADERS Current_Base.H
			   VA_0_PiPi.H
			   VA_0_PiPiPi.H
			   VA_0_KKPi.H
			   VA_0_KPiK.H
			   VA_0_KPiPi.H
			   Tools.H )


sherpa_create_git_info(HadronCurrents)
 
add_library(METoolsHadronCurrents SHARED ${HadronCurrents_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(METoolsHadronCurrents METoolsHadronCurrents)
