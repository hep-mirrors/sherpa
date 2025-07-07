set(SoftPhysics_SOURCES Hadron_Init.C
                        Beam_Remnant_Handler.C
                        Colour_Reconnection_Handler.C
                        Soft_Collision_Handler.C
                        Hadron_Decay_Handler.C
                        Hadron_Rescattering_Handler.C
                        Soft_Photon_Handler.C
                        Singlet_Sorter.C
                        Resonance_Finder.C)

set(SoftPhysics_HEADERS Hadron_Init.H
                        Beam_Remnant_Handler.H
                        Colour_Reconnection_Handler.H
                        Soft_Collision_Handler.H
                        Hadron_Decay_Handler.H
                        Hadron_Rescattering_Handler.H
                        Soft_Photon_Handler.H
                        Singlet_Sorter.H
                        Resonance_Finder.H)

sherpa_create_git_info(SoftPhysics ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaSoftPhysics SHARED ${SoftPhysics_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaSoftPhysics SherpaSoftPhysics)
add_dependencies(SherpaSoftPhysics SoftPhysics_Git_Info)
