set(RemnantsMain_SOURCES Remnant_Handler.C
                 Remnant_Base.C
                 Hadron_Remnant.C
                 Electron_Remnant.C
                 Photon_Remnant.C
                 Pomeron_Remnant.C
                 Reggeon_Remnant.C
                 No_Remnant.C)

set(RemnantsMain_HEADERS Remnant_Handler.H
                 Remnant_Base.H
                 Hadron_Remnant.H
                 Electron_Remnant.H
                 Photon_Remnant.H
                 Pomeron_Remnant.H
                 Reggeon_Remnant.H
                 No_Remnant.H)

sherpa_create_git_info(RemnantsMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(Remnants SHARED ${RemnantsMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(Remnants Remnants)
add_dependencies(Remnants RemnantsMain_Git_Info)
