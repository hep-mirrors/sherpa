set(Main_SOURCES Color_Integrator.C
                 Helicity_Integrator.C
                 Process_Integrator.C
                 Channel_Creator.C
                 Phase_Space_Point.C
                 Phase_Space_Enhance.C
                 Phase_Space_Handler.C
                 Phase_Space_Integrator.C)

set(Main_HEADERS Color_Integrator.H
                 Helicity_Integrator.H
                 Process_Integrator.H
                 Channel_Creator.H
                 Phase_Space_Point.H
                 Phase_Space_Enhance.H
                 Phase_Space_Variation_Weights.H
                 Phase_Space_Handler.H
                 Phase_Space_Integrator.H)

sherpa_create_git_info(Main)
 
add_library(PhasicMain SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(PhasicMain PhasicMain)
