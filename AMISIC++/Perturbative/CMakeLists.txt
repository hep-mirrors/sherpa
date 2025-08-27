set(Perturbative_SOURCES Single_Collision_Handler.C
                         MI_Processes.C
                         MI_Process_Group.C
                         MI_Process.C
                         MI_Integrator.C
                         QCD_Processes.C
                         QED_Processes.C)

set(Perturbative_HEADERS Single_Collision_Handler.H
                         MI_Processes.H
                         MI_Process_Group.H
                         MI_Process.H
                         MI_Integrator.H
                         QCD_Processes.H
                         QED_Processes.H)

sherpa_create_git_info(Perturbative ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(AmisicPerturbative SHARED ${Perturbative_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(AmisicPerturbative AmisicPerturbative)
add_dependencies(AmisicPerturbative Perturbative_Git_Info)
