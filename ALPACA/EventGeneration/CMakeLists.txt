set(EventGeneration_SOURCES Analysis_Handler.C
                  Colour_Handler.C
                  Dynamic_Quantities_Handler.C
                  Elastic_Kernel.C
                  Event_Generator.C
                  Inelastic_Kernel.C
                  Scatter_Merge_Handler.C
                  Simple_MEs.C
                  Split_Handler.C)

set(EventGeneration_HEADERS Analysis_Handler.H
                  Colour_Handler.H
                  Dynamic_Quantities_Handler.H
                  Elastic_Kernel.H
                  Event_Generator.H
                  Inelastic_Kernel.H
                  Scatter_Merge_Handler.H
                  Simple_MEs.H
                  Split_Handler.H)

sherpa_create_git_info(EventGeneration)
 
add_library(AlpacaEventGeneration SHARED ${EventGeneration_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )

sherpa_handle_shared_library(AlpacaEventGeneration AlpacaEventGeneration)

message(STATUS "EventGeneration_SOURCES = ${EventGeneration_SOURCES}")

message(STATUS "Looking for Event_Generator.C = ${EventGeneration_SOURCES}")
