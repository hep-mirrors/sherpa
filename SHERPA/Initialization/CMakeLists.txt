
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Initialization/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaInitialization.la)
 set(GITTAG  SHERPA/Initialization)
 #include ../../git.make)
set(Initialization_SOURCES  Initialization_Handler.C)
set(Initialization_HEADERS Initialization_Handler.H)
creategitinfo(Initialization ../Initialization/)
  list(TRANSFORM Initialization_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Initialization/")
 
set(Initialization_esources )
if (SHERPA_BUILD_STATIC)
add_library(SherpaInitialization_static STATIC ${SherpaInitialization_esources}
                             ${Initialization_SOURCES} 
)
target_include_directories(SherpaInitialization_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaInitialization_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaInitialization_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaInitialization_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME SherpaInitialization)
endif()
if (SHERPA_BUILD_SHARED)
add_library(SherpaInitialization SHARED ${SherpaInitialization_esources}
                             ${Initialization_SOURCES} 
)
target_include_directories(SherpaInitialization PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaInitialization PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaInitialization DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaInitialization        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaInitialization SOVERSION 0.0.0)
endif()
