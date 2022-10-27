
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Gauge/Makefile.am
 
 #pkglib_LTLIBRARIES = libDIMGauge.la)
 set(GITTAG  DIM/Gauge)
 #include ../../git.make)
set(Gauge_SOURCES  GGG.C GQQ.C QQG.C)
creategitinfo(Gauge ../Gauge/)
  list(TRANSFORM Gauge_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Gauge/")
 
set(Gauge_esources )
if (SHERPA_BUILD_STATIC)
add_library(DIMGauge_static STATIC ${DIMGauge_esources}
                             ${Gauge_SOURCES} 
)
target_include_directories(DIMGauge_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(DIMGauge_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS DIMGauge_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DIMGauge_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME DIMGauge)
endif()
if (SHERPA_BUILD_SHARED)
add_library(DIMGauge SHARED ${DIMGauge_esources}
                             ${Gauge_SOURCES} 
)
target_include_directories(DIMGauge PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(DIMGauge PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS DIMGauge DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DIMGauge        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DIMGauge SOVERSION 0.0.0)
endif()
