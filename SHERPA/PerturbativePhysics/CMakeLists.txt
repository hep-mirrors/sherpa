set(PerturbativePhysics_SOURCES Matrix_Element_Handler.C
                                MI_Handler.C
                                Shower_Handler.C
                                Perturbative_Interface.C
                                Hard_Decay_Handler.C)

set(PerturbativePhysics_HEADERS Matrix_Element_Handler.H
                                MI_Handler.H
                                Shower_Handler.H
                                Perturbative_Interface.H
                                Hard_Decay_Handler.H)

sherpa_create_git_info(PerturbativePhysics ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaPerturbativePhysics SHARED ${PerturbativePhysics_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaPerturbativePhysics SherpaPerturbativePhysics)
add_dependencies(SherpaPerturbativePhysics PerturbativePhysics_Git_Info)
