set(HadronCurrents_SOURCES Line_Shapes.C
                           Resonance_Base.C
                           Propagator.C
                           Vector_Decays.C
                           Scalar_Decays.C
                           F0_Decays.C
                           Rho_Decays.C
                           Omega_Decays.C
                           A1_Decays.C)
set(HadronCurrents_HEADERS Line_Shapes.H
			   Resonance_Base.H
                           Propagator.H
                           Vector_Decays.H
                           Scalar_Decays.H
                           F0_Decays.H
                           Rho_Decays.H
                           Omega_Decays.H
                           A1_Decays.H )


sherpa_create_git_info(HadronCurrents)
 
add_library(METoolsHadronCurrents SHARED ${HadronCurrents_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(METoolsHadronCurrents METoolsHadronCurrents)
