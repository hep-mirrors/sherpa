set(Main_SOURCES Reconnection_Handler.C
                 Reconnection_Base.C
                 Reconnect_By_Singlet.C
                 Reconnect_Statistical.C
                 Reconnection_Weights.C)

set(Main_HEADERS Reconnection_Handler.H
                 Reconnection_Base.H
                 Reconnect_By_Singlet.H
                 Reconnect_Statistical.H
                 Reconnection_Weights.H)

sherpa_create_git_info(Main)
 
add_library(Reconnections SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(Reconnections Reconnections)
