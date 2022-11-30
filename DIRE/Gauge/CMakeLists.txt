##### Automatically or semiautomaticaly generated, do not edit.#########
#Gauge/Makefile.am
 
set(GITTAG  DIRE/Gauge)

set(Gauge_SOURCES GGG.C
                  GQQ.C
                  QQG.C
                  QQQQ.C)

sherpa_create_git_info(Gauge)
 
add_library(DireGauge SHARED ${Gauge_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DireGauge DireGauge)
