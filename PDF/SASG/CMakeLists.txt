
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#SASG/Makefile.am
 
 #pkglib_LTLIBRARIES = libSASGSherpa.la)
 set(GITTAG  PDF/SASG)
 #include ../../git.make)
set(SASG_SOURCES  sasgam2.f SASGph_Fortran_Interface.C)
set(SASG_HEADERS SASGph_Fortran_Interface.H)
creategitinfo(SASG ../SASG/)
  list(TRANSFORM SASG_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../SASG/")
 
set(SASG_esources )
if (SHERPA_BUILD_STATIC)
add_library(SASGSherpa_static STATIC ${SASGSherpa_esources}
                             ${SASG_SOURCES} 
)
target_include_directories(SASGSherpa_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SASGSherpa_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SASGSherpa_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SASGSherpa_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME SASGSherpa)
endif()
if (SHERPA_BUILD_SHARED)
add_library(SASGSherpa SHARED ${SASGSherpa_esources}
                             ${SASG_SOURCES} 
)
target_include_directories(SASGSherpa PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SASGSherpa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SASGSherpa DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SASGSherpa        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SASGSherpa SOVERSION 0.0.0)
endif()
