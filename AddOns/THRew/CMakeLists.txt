set(THRew_SOURCES KFactor.C
                  Variation_Generator.C)

set(THRew_HEADERS KFactor.H
                  Variation_Generator.H)

sherpa_create_git_info(THRew ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(SherpaTHRew SHARED ${THRew_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
sherpa_handle_shared_library(SherpaTHRew SherpaTHRew)
add_dependencies(SherpaTHRew THRew_Git_Info)
