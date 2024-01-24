set(Beam_Remnants_SOURCES Hadron_Dissociation.C
                          Remnant_Handler.C)

set(Beam_Remnants_HEADERS Hadron_Dissociation.H
                          Remnant_Handler.H)

sherpa_create_git_info(Beam_Remnants)
 
add_library(ShrimpsBeamRemnants SHARED ${Beam_Remnants_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ShrimpsBeamRemnants ShrimpsBeamRemnants)
