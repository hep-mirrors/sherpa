set(AlaricMain_SOURCES Jet_Criterion.C
                 Alaric.C)

sherpa_create_git_info(AlaricMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(AlaricMain SHARED ${AlaricMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(AlaricMain AlaricMain)
