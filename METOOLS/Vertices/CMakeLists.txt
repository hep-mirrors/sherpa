
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Vertices/Makefile.am
 
 #pkglib_LTLIBRARIES = libMEToolsVertices.la)
 set(GITTAG  METOOLS/Vertices)
 #include ../../git.make)
set(Vertices_SOURCES  FFV_DLC.C VVV_DLC.C SSS_LC.C SSV_LC.C SFF_LC.C SVV_LC.C FFV_LC.C VVV_LC.C VVP_LC.C SSSS_LC.C SSVV_LC.C VVVV_LC.C)
creategitinfo(Vertices ../Vertices/)
  list(TRANSFORM Vertices_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Vertices/")
 
set(Vertices_esources )
if (SHERPA_BUILD_STATIC)
add_library(METoolsVertices_static STATIC ${METoolsVertices_esources}
                             ${Vertices_SOURCES} 
)
target_include_directories(METoolsVertices_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(METoolsVertices_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS METoolsVertices_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(METoolsVertices_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME METoolsVertices)
set_target_properties(METoolsVertices_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(METoolsVertices SHARED ${METoolsVertices_esources}
                             ${Vertices_SOURCES} 
)
target_include_directories(METoolsVertices PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(METoolsVertices PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS METoolsVertices DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(METoolsVertices        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME METoolsVertices SOVERSION 0.0.0)
set_target_properties(METoolsVertices        PROPERTIES DEFINE_SYMBOL "")
endif()
