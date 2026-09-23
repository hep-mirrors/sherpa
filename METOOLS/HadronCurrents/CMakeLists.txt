set(HadronCurrents_SOURCES Current_Base.C
			   VA_0_PiPi.C
			   V_0_Isoscalar3Pi.C VA_0_PiPiPi.C
			   VA_0_KKPi.C
			   VA_0_KPiK.C
			   VA_0_KPiPi.C
			   VA_0_EtaPiPi.C
			   VA_0_KOmega3Pi.C
			   VA_0_Novosibirsk4Pi.C
			   VA_0_FF.C
			   VA_P_0.C
			   VA_P_X_Base.C
			   VA_P_PP.C
			   VA_P_P.C
			   VA_P_S.C
			   VA_P_V.C
			   VA_P_A.C
			   VA_P_T.C
			   Tools.C )
set(HadronCurrents_HEADERS Current_Base.H
			   VA_0_PiPi.H
			   V_0_Isoscalar3Pi.H VA_0_PiPiPi.H
			   VA_0_KKPi.H
			   VA_0_KPiK.H
			   VA_0_KPiPi.H
			   VA_0_EtaPiPi.H
			   VA_0_KOmega3Pi.H
			   VA_0_Novosibirsk4Pi.H
			   VA_0_FF.H
			   VA_P_0.H
			   VA_P_X_Base.H
			   VA_P_PP.H
			   VA_P_P.H
			   VA_P_S.H
			   VA_P_V.H
			   VA_P_A.H
			   VA_P_T.H
			   Tools.H )


sherpa_create_git_info(HadronCurrents)
 
add_library(METoolsHadronCurrents SHARED ${HadronCurrents_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(METoolsHadronCurrents METoolsHadronCurrents)
