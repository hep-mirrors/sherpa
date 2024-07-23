set(Main_SOURCES Remnant_Handler.C
                 Remnant_Base.C
                 Hadron_Remnant.C
                 Electron_Remnant.C
                 Photon_Remnant.C
                 Pomeron_Remnant.C
                 Reggeon_Remnant.C
                 No_Remnant.C)

set(Main_HEADERS Remnant_Handler.H
                 Remnant_Base.H
                 Hadron_Remnant.H
                 Electron_Remnant.H
                 Photon_Remnant.H
                 Pomeron_Remnant.H
                 Reggeon_Remnant.H
                 No_Remnant.H)

sherpa_create_git_info(Main)

add_library(Remnants SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(Remnants Remnants)
