set(Gauge_SOURCES GGG.C
                  GQQ.C
                  QQG.C)

sherpa_create_git_info(Gauge)
 
add_library(DIMGauge SHARED ${Gauge_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DIMGauge DIMGauge)
