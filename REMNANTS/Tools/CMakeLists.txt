set(RemnantsTools_SOURCES Kinematics_Generator.C
                  Colour_Generator.C
                  Beam_Decorrelator.C
                  Remnants_Parameters.C
                  Primordial_KPerp.C
		  Continued_PDF.C
		  Pseudo_TMD.C
                  Form_Factor.C)

set(RemnantsTools_HEADERS Kinematics_Generator.H
                  Colour_Generator.H
                  Beam_Decorrelator.H
                  Remnants_Parameters.H
                  Primordial_KPerp.H
		  Continued_PDF.H
		  Pseudo_TMD.H
                  Form_Factor.H)

sherpa_create_git_info(RemnantsTools ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(RemnantsTools SHARED ${RemnantsTools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(RemnantsTools RemnantsTools)
add_dependencies(RemnantsTools RemnantsTools_Git_Info)
