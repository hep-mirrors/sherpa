set(Main_SOURCES Jet_Criterion.C
                 Alaric.C)

sherpa_create_git_info(Main)
 
add_library(AlaricMain SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(AlaricMain AlaricMain)
