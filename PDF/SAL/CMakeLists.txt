
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#SAL/Makefile.am
 
 #pkglib_LTLIBRARIES = libSALSherpa.la)
 set(GITTAG  PDF/SAL)
 #include ../../git.make)
set(SAL_SOURCES  sal.cpp SALph_CPP_Interface.C)
set(SAL_HEADERS SALph_CPP_Interface.H sal.h math_ut.h)
#nobase_#dist_pkgdata_DATA = SALGrid/SAL.dat
creategitinfo(SAL ../SAL/)
  list(TRANSFORM SAL_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../SAL/")
 
set(SAL_esources )
if (SHERPA_BUILD_STATIC)
add_library(SALSherpa_static STATIC ${SALSherpa_esources}
                             ${SAL_SOURCES} 
)
target_include_directories(SALSherpa_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SALSherpa_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SALSherpa_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SALSherpa_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME SALSherpa)
set_target_properties(SALSherpa_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(SALSherpa SHARED ${SALSherpa_esources}
                             ${SAL_SOURCES} 
)
target_include_directories(SALSherpa PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SALSherpa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SALSherpa DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SALSherpa        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SALSherpa SOVERSION 0.0.0)
set_target_properties(SALSherpa        PROPERTIES DEFINE_SYMBOL "")
endif()
