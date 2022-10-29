
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Special/Makefile.am
 
 #pkglib_LTLIBRARIES = libExtraXSSpecial.la)
 set(GITTAG  EXTRA_XS/Special)
 #include ../../git.make)
set(Special_SOURCES  Instanton.C)
creategitinfo(Special ../Special/)
  list(TRANSFORM Special_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Special/")
 
set(Special_esources )
if (SHERPA_BUILD_STATIC)
add_library(ExtraXSSpecial_static STATIC ${ExtraXSSpecial_esources}
                             ${Special_SOURCES} 
)
target_include_directories(ExtraXSSpecial_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ExtraXSSpecial_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ExtraXSSpecial_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ExtraXSSpecial_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME ExtraXSSpecial)
set_target_properties(ExtraXSSpecial_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(ExtraXSSpecial SHARED ${ExtraXSSpecial_esources}
                             ${Special_SOURCES} 
)
target_include_directories(ExtraXSSpecial PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ExtraXSSpecial PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ExtraXSSpecial DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ExtraXSSpecial        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ExtraXSSpecial SOVERSION 0.0.0)
set_target_properties(ExtraXSSpecial        PROPERTIES DEFINE_SYMBOL "")
endif()
