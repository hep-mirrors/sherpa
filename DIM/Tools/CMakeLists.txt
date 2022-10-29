
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Tools/Makefile.am
 
 #pkglib_LTLIBRARIES = libDIMTools.la)
 set(GITTAG  DIM/Tools)
 #include ../../git.make)
set(Tools_SOURCES  Color.C Parton.C Splitting.C Amplitude.C)
set(Tools_HEADERS Color.H Parton.H Weight.H Splitting.H Kernel_Key.H Amplitude.H)
creategitinfo(Tools ../Tools/)
  list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
 
set(Tools_esources )
if (SHERPA_BUILD_STATIC)
add_library(DIMTools_static STATIC ${DIMTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(DIMTools_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(DIMTools_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS DIMTools_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DIMTools_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME DIMTools)
set_target_properties(DIMTools_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(DIMTools SHARED ${DIMTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(DIMTools PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(DIMTools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS DIMTools DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DIMTools        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DIMTools SOVERSION 0.0.0)
set_target_properties(DIMTools        PROPERTIES DEFINE_SYMBOL "")
endif()
