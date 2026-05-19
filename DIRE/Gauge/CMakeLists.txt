set(DireGauge_SOURCES GGG.C
                  GQQ.C
                  QQG.C
                  QQQQ.C)

sherpa_create_git_info(DireGauge ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(DireGauge SHARED ${DireGauge_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DireGauge DireGauge)
add_dependencies(DireGauge DireGauge_Git_Info)
