##### Automatically or semiautomaticaly generated, do not edit.#########
#Gauge/Makefile.am
 
 #pkglib_LTLIBRARIES = libDIMGauge.la)
 set(GITTAG  DIM/Gauge)
set(Gauge_SOURCES  GGG.C GQQ.C QQG.C)
creategitinfo(Gauge ../Gauge/)
list(TRANSFORM Gauge_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Gauge/")
list(APPEND Gauge_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(DIMGauge SHARED ${Gauge_SOURCES} )
target_include_directories(DIMGauge PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(DIMGauge)
install(TARGETS DIMGauge DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DIMGauge        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DIMGauge SOVERSION 0.0.0)
set_target_properties(DIMGauge        PROPERTIES DEFINE_SYMBOL "")
