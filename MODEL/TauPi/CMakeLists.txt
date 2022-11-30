set(TauPi_SOURCES TauPi_LC.C
                  Model.C)

sherpa_create_git_info(TauPi)
 
add_library(SherpaTauPi SHARED ${TauPi_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaTauPi SherpaTauPi)
