set(Eikonals_SOURCES Eikonal_Weights.C
                     Rapidity_Density.C
                     Omega_ik.C
                     Single_Channel_Eikonal.C
                     Eikonal_Contributor.C
                     Grid.C
                     Analytic_Eikonal.C
                     Analytic_Contributor.C
                     Eikonal_Creator.C
                     Form_Factors.C)

set(Eikonals_HEADERS Eikonal_Weights.H
                     Rapidity_Density.H
                     Omega_ik.H
                     Single_Channel_Eikonal.H
                     Eikonal_Contributor.H
                     Grid.H
                     Analytic_Eikonal.H
                     Analytic_Contributor.H
                     Eikonal_Creator.H
                     Form_Factors.H)

sherpa_create_git_info(Eikonals ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(ShrimpsEikonals SHARED ${Eikonals_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ShrimpsEikonals ShrimpsEikonals)
add_dependencies(ShrimpsEikonals Eikonals_Git_Info)
