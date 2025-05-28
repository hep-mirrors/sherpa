set(PhasicDecays_SOURCES Decay_Channel.C
                   Decay_Table.C
                   Decay_Map.C)

set(PhasicDecays_HEADERS Decay_Channel.H
                   Decay_Table.H
                   Decay_Map.H)

sherpa_create_git_info(PhasicDecays ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(PhasicDecays SHARED ${PhasicDecays_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(PhasicDecays PhasicDecays)
add_dependencies(PhasicDecays PhasicDecays_Git_Info)
