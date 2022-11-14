##### Automatically or semiautomaticaly generated, do not edit.#########
#Gauge/Makefile.am
 
 #pkglib_LTLIBRARIES = libDIMGauge.la)
 set(GITTAG  DIM/Gauge)
set(Gauge_SOURCES  GGG.C GQQ.C QQG.C)
sherpa_create_git_info(Gauge ../Gauge/)
list(TRANSFORM Gauge_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Gauge/")
list(APPEND Gauge_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(DIMGauge SHARED ${Gauge_SOURCES} )
sherpa_handle_shared_library(DIMGauge DIMGauge)
