
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libCSMain.la)
 set(GITTAG  CSSHOWER++/Main)
 #include ../../git.make)
set(Main_SOURCES  CS_Cluster_Definitions.C CS_Shower.C)
set(Main_HEADERS CS_Cluster_Definitions.H CS_Shower.H)
creategitinfo(Main ../Main/)
  list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
 
set(Main_esources )
if (SHERPA_BUILD_STATIC)
add_library(CSMain_static STATIC ${CSMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(CSMain_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(CSMain_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS CSMain_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(CSMain_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME CSMain)
set_target_properties(CSMain_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(CSMain SHARED ${CSMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(CSMain PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(CSMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS CSMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(CSMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME CSMain SOVERSION 0.0.0)
set_target_properties(CSMain        PROPERTIES DEFINE_SYMBOL "")
endif()
