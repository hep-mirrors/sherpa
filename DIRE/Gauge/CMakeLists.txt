
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Gauge/Makefile.am
 
 #pkglib_LTLIBRARIES = libDireGauge.la)
 set(GITTAG  DIRE/Gauge)
 #include ../../git.make)
set(Gauge_SOURCES  GGG.C GQQ.C QQG.C QQQQ.C)
creategitinfo(Gauge ../Gauge/)
list(TRANSFORM Gauge_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Gauge/")
list(APPEND Gauge_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(Gauge_esources )
if (SHERPA_BUILD_SHARED)
add_library(DireGauge SHARED ${DireGauge_esources}
                             ${Gauge_SOURCES} 
)
target_include_directories(DireGauge PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(DireGauge PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS DireGauge DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DireGauge        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DireGauge SOVERSION 0.0.0)
set_target_properties(DireGauge        PROPERTIES DEFINE_SYMBOL "")
endif()
