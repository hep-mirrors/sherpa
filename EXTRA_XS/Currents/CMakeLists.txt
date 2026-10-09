set(Currents_SOURCES XS_Currents.C)

sherpa_create_git_info(Currents)

add_library(ExtraXSCurrents SHARED ${Currents_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ExtraXSCurrents ExtraXSCurrents)
