
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Tools/Makefile.am
 
 #pkglib_LTLIBRARIES = libCSTools.la)
 set(GITTAG  CSSHOWER++/Tools)
 #include ../../git.make)
set(Tools_SOURCES  Parton.C Singlet.C)
set(Tools_HEADERS Parton.H Singlet.H)
creategitinfo(Tools ../Tools/)
  list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
 
set(Tools_esources )
if (SHERPA_BUILD_STATIC)
add_library(CSTools_static STATIC ${CSTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(CSTools_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(CSTools_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS CSTools_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(CSTools_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME CSTools)
endif()
if (SHERPA_BUILD_SHARED)
add_library(CSTools SHARED ${CSTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(CSTools PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(CSTools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS CSTools DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(CSTools        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME CSTools SOVERSION 0.0.0)
endif()
