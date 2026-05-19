set(DIMGauge_SOURCES GGG.C
                  GQQ.C
                  QQG.C)

sherpa_create_git_info(DIMGauge ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(DIMGauge SHARED ${DIMGauge_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DIMGauge DIMGauge)
add_dependencies(DIMGauge DIMGauge_Git_Info)
