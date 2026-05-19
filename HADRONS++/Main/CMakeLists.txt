set(HadronsMain_SOURCES Hadron_Decay_Map.C
                 Hadron_Decay_Table.C
                 Hadron_Decay_Channel.C
                 Mixing_Handler.C
                 Tools.C)

set(HadronsMain_HEADERS Hadron_Decay_Map.H
                 Hadron_Decay_Table.H
                 Hadron_Decay_Channel.H
                 Mixing_Handler.H
                 Tools.H)

sherpa_create_git_info(HadronsMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(HadronsMain SHARED ${HadronsMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(HadronsMain HadronsMain)
add_dependencies(HadronsMain HadronsMain_Git_Info)
