set(METoolsFormFactors_SOURCES Line_Shapes.C
                               Propagator.C
                               Resonance_Base.C
			       Scalar_Decays.C
			       Vector_Decays.C
			       A1_Decays.C
			       F0_Decays.C
			       Kstar_Decays.C
			       Omega_Decays.C
			       Rho_Decays.C
			     )
set(METoolsFormFactors_HEADERS Line_Shapes.C
                               Propagator.H
                               Resonance_Base.H
			       Scalar_Decays.H
			       Vector_Decays.H
			       A1_Decays.H
			       F0_Decays.H
			       Kstar_Decays.H
			       Omega_Decays.H
			       Rho_Decays.H
			     )
			       
sherpa_create_git_info(METoolsFormFactors ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(METoolsFormFactors SHARED ${METoolsFormFactors_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(METoolsFormFactors METoolsFormFactors)
add_dependencies(METoolsFormFactors METoolsFormFactors_Git_Info)
