set(ReconnectionsMain_SOURCES Reconnection_Handler.C
                 Reconnection_Base.C
                 Reconnect_By_Singlet.C
                 Reconnect_Statistical.C
                 Reconnection_Weights.C)

set(ReconnectionsMain_HEADERS Reconnection_Handler.H
                 Reconnection_Base.H
                 Reconnect_By_Singlet.H
                 Reconnect_Statistical.H
                 Reconnection_Weights.H)

sherpa_create_git_info(ReconnectionsMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(Reconnections SHARED ${ReconnectionsMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(Reconnections Reconnections)
add_dependencies(Reconnections ReconnectionsMain_Git_Info)
