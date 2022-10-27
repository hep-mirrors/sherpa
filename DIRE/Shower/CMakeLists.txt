
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Shower/Makefile.am
 
 #pkglib_LTLIBRARIES = libDireShower.la)
 set(GITTAG  DIRE/Shower)
 #include ../../git.make)
set(Shower_SOURCES  Gauge.C Lorentz.C Kernel.C Alpha_QCD.C Lorentz_FF.C Lorentz_FI.C Lorentz_IF.C Lorentz_II.C Cluster_Definitions.C Shower.C)
set(Shower_HEADERS Gauge.H Lorentz.H Kernel.H Alpha_QCD.H Lorentz_FF.H Lorentz_FI.H Lorentz_IF.H Lorentz_II.H Cluster_Definitions.H Shower.H)
creategitinfo(Shower ../Shower/)
  list(TRANSFORM Shower_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Shower/")
 
set(Shower_esources )
if (SHERPA_BUILD_STATIC)
add_library(DireShower_static STATIC ${DireShower_esources}
                             ${Shower_SOURCES} 
)
target_include_directories(DireShower_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(DireShower_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS DireShower_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DireShower_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME DireShower)
endif()
if (SHERPA_BUILD_SHARED)
add_library(DireShower SHARED ${DireShower_esources}
                             ${Shower_SOURCES} 
)
target_include_directories(DireShower PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(DireShower PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS DireShower DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DireShower        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DireShower SOVERSION 0.0.0)
endif()
