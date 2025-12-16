set(Beam_Remnants_SOURCES Continued_PDF.C
                          Pseudo_TMD.C
                          Hadron_Dissociation.C
                          Remnant_Handler.C
                          Colour_Generator.C)

set(Beam_Remnants_HEADERS Continued_PDF.H
                          Pseudo_TMD.H
                          Hadron_Dissociation.H
                          Remnant_Handler.H
                          Colour_Generator.H)

sherpa_create_git_info(Beam_Remnants ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(ShrimpsBeamRemnants SHARED ${Beam_Remnants_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ShrimpsBeamRemnants ShrimpsBeamRemnants)
add_dependencies(ShrimpsBeamRemnants Beam_Remnants_Git_Info)
