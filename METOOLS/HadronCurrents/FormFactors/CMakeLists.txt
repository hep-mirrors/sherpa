set(FormFactors_SOURCES Line_Shapes.C
                        Resonance_Base.C
                        Propagator.C
                        Vector_Decays.C
                        Scalar_Decays.C
                        F0_Decays.C
                        Rho_Decays.C
                        Omega_Decays.C
                        A1_Decays.C
                        Kstar_Decays.C
                        K1_Decays.C
                        A0_Decays.C
                        K_0_800_Decays.C
                        FormFactor_Base.C
			FF_0_PP.C
			FF_0_PPP.C
			FF_0_DD.C
			RChL_Functions.C
                     Novosibirsk4pi_GTables.C
		      )
set(FormFactors_HEADERS Line_Shapes.H
                        Resonance_Base.H
                        Propagator.H
                        Vector_Decays.H
                        Scalar_Decays.H
                        F0_Decays.H
                        Rho_Decays.H
                        Omega_Decays.H
                        A1_Decays.H
                        Kstar_Decays.H
                        K1_Decays.H
                        A0_Decays.H
                        K_0_800_Decays.H
                        FormFactor_Base.H
			FF_0_PP.H
			FF_0_PPP.H
			FF_0_DD.H
			RChL_Functions.H
			Novosibirsk4pi_GTables.H
		      )


sherpa_create_git_info(FormFactors)
 
add_library(METoolsFormFactors SHARED ${FormFactors_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(METoolsFormFactors METoolsFormFactors)
