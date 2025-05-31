set(AlaricGauge_SOURCES GGG.C
                  GQQ.C
                  QQG.C)

sherpa_create_git_info(AlaricGauge ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(AlaricGauge SHARED ${AlaricGauge_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(AlaricGauge AlaricGauge)
