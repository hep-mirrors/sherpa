set(EXTAMP_SOURCES External_ME_Interface.C
                   Process.C
                   Born_Process.C
                   RS_Process.C
                   BVI_Process.C
                   CS_Dipole.C
                   CS_Dipoles.C
                   Dipole_Wrapper_Process.C)

set(EXTAMP_HEADERS External_ME_Interface.H
                   Process_Group.H
                   Process.H
                   Born_Process.H
                   RS_Process.H
                   BVI_Process.H
                   CS_Dipole.H
                   CS_Dipoles.H
                   Dipole_Wrapper_Process.H)

sherpa_create_git_info(EXTAMP ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(ExtAmp SHARED ${EXTAMP_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ExtAmp ExtAmp)
