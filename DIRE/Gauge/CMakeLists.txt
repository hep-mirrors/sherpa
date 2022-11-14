##### Automatically or semiautomaticaly generated, do not edit.#########
#Gauge/Makefile.am
 
 #pkglib_LTLIBRARIES = libDireGauge.la)
 set(GITTAG  DIRE/Gauge)
set(Gauge_SOURCES  GGG.C GQQ.C QQG.C QQQQ.C)
creategitinfo(Gauge ../Gauge/)
list(TRANSFORM Gauge_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Gauge/")
list(APPEND Gauge_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(DireGauge SHARED ${Gauge_SOURCES} )
target_include_directories(DireGauge PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(DireGauge)
install(TARGETS DireGauge DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DireGauge        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DireGauge SOVERSION 0.0.0)
set_target_properties(DireGauge        PROPERTIES DEFINE_SYMBOL "")
