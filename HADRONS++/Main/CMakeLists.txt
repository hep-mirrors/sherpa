set(Main_SOURCES Hadron_Decay_Map.C
                 Hadron_Decay_Table.C
                 Hadron_Decay_Channel.C
                 Mixing_Handler.C
                 Tools.C)

set(Main_HEADERS Hadron_Decay_Map.H
                 Hadron_Decay_Table.H
                 Hadron_Decay_Channel.H
                 Mixing_Handler.H
                 Tools.H)

sherpa_create_git_info(Main)
 
add_library(HadronsMain SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(HadronsMain HadronsMain)
