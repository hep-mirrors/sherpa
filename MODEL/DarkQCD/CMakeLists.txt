set(DarkQCD_SOURCES  Model.C
                     Calculators/CF_DarkQCD.C
                     Scales/DarkQCD_Scale_Setter.C)

set(DarkQCD_HEADERS Model.H)

sherpa_create_git_info(DarkQCD ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaDarkQCD SHARED ${DarkQCD_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaDarkQCD SherpaDarkQCD)
add_dependencies(SherpaDarkQCD DarkQCD_Git_Info)
