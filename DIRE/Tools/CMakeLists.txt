
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Tools/Makefile.am
 
 #pkglib_LTLIBRARIES = libDireTools.la)
 set(GITTAG  DIRE/Tools)
 #include ../../git.make)
set(Tools_SOURCES  Color.C Parton.C Splitting.C Amplitude.C)
set(Tools_HEADERS Color.H Parton.H Weight.H Splitting.H Kernel_Key.H Amplitude.H)
creategitinfo(Tools ../Tools/)
  list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
 
set(Tools_esources )
if (SHERPA_BUILD_STATIC)
add_library(DireTools_static STATIC ${DireTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(DireTools_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(DireTools_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS DireTools_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DireTools_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME DireTools)
endif()
if (SHERPA_BUILD_SHARED)
add_library(DireTools SHARED ${DireTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(DireTools PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(DireTools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS DireTools DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DireTools        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DireTools SOVERSION 0.0.0)
endif()
