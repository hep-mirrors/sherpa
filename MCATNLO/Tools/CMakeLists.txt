
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Tools/Makefile.am
 
 #pkglib_LTLIBRARIES = libMCatNLOTools.la)
 set(GITTAG  MCATNLO/Tools)
 #include ../../git.make)
set(Tools_SOURCES  Parton.C Singlet.C)
set(Tools_HEADERS Parton.H Singlet.H)
creategitinfo(Tools ../Tools/)
  list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
 
set(Tools_esources )
if (SHERPA_BUILD_STATIC)
add_library(MCatNLOTools_static STATIC ${MCatNLOTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(MCatNLOTools_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(MCatNLOTools_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS MCatNLOTools_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(MCatNLOTools_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME MCatNLOTools)
endif()
if (SHERPA_BUILD_SHARED)
add_library(MCatNLOTools SHARED ${MCatNLOTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(MCatNLOTools PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(MCatNLOTools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS MCatNLOTools DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(MCatNLOTools        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME MCatNLOTools SOVERSION 0.0.0)
endif()
