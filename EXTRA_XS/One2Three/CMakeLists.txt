
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#One2Three/Makefile.am
 
 #pkglib_LTLIBRARIES = libExtraXS1_3.la)
 set(GITTAG  EXTRA_XS/One2Three)
 #include ../../git.make)
set(One2Three_SOURCES  Comix1to3.C)
set(One2Three_HEADERS Comix1to3.H)
creategitinfo(One2Three ../One2Three/)
  list(TRANSFORM One2Three_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../One2Three/")
 
set(One2Three_esources )
if (SHERPA_BUILD_STATIC)
add_library(ExtraXS1_3_static STATIC ${ExtraXS1_3_esources}
                             ${One2Three_SOURCES} 
)
target_include_directories(ExtraXS1_3_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ExtraXS1_3_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ExtraXS1_3_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ExtraXS1_3_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME ExtraXS1_3)
endif()
if (SHERPA_BUILD_SHARED)
add_library(ExtraXS1_3 SHARED ${ExtraXS1_3_esources}
                             ${One2Three_SOURCES} 
)
target_include_directories(ExtraXS1_3 PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ExtraXS1_3 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ExtraXS1_3 DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ExtraXS1_3        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ExtraXS1_3 SOVERSION 0.0.0)
endif()
