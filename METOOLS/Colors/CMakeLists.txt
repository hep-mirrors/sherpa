
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Colors/Makefile.am
 
 #pkglib_LTLIBRARIES = libMEToolsColors.la)
 set(GITTAG  METOOLS/Colors)
 #include ../../git.make)
set(Colors_SOURCES  None_CC.C D_CC.C T_CC.C F_CC.C FF_CC.C TT_CC.C ST_CC.C SF_CC.C)
set(Colors_HEADERS)
creategitinfo(Colors ../Colors/)
  list(TRANSFORM Colors_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Colors/")
 
set(Colors_esources )
if (SHERPA_BUILD_STATIC)
add_library(METoolsColors_static STATIC ${METoolsColors_esources}
                             ${Colors_SOURCES} 
)
target_include_directories(METoolsColors_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(METoolsColors_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS METoolsColors_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(METoolsColors_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME METoolsColors)
set_target_properties(METoolsColors_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(METoolsColors SHARED ${METoolsColors_esources}
                             ${Colors_SOURCES} 
)
target_include_directories(METoolsColors PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(METoolsColors PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS METoolsColors DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(METoolsColors        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME METoolsColors SOVERSION 0.0.0)
set_target_properties(METoolsColors        PROPERTIES DEFINE_SYMBOL "")
endif()
