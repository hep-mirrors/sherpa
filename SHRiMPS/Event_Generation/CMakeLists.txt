set(Event_Generation_SOURCES Event_Generator.C
                             Inelastic_Event_Generator.C
                             Elastic_Event_Generator.C
                             Soft_Diffractive_Event_Generator.C
                             Collinear_Emission_Generator.C)

set(Event_Generation_HEADERS Event_Generator_Base.H
                             Event_Generator.H
                             Inelastic_Event_Generator.H
                             Elastic_Event_Generator.H
                             Soft_Diffractive_Event_Generator.H
                             Collinear_Emission_Generator.H)

sherpa_create_git_info(Event_Generation)
 
add_library(ShrimpsEvents SHARED ${Event_Generation_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ShrimpsEvents ShrimpsEvents)
