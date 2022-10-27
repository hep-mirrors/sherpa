
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#LHAPDF/Makefile.am
 
 set(GITTAG  PDF/LHAPDF)
 #include ../../git.make)
 #pkglib_LTLIBRARIES = libLHAPDFSherpa.la)
set(LHAPDF_SOURCES  LHAPDF_CPP_Interface.C)
creategitinfo(LHAPDF ../LHAPDF/)
  list(TRANSFORM LHAPDF_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../LHAPDF/")
 
set(LHAPDF_esources )
if (SHERPA_BUILD_STATIC)
add_library(LHAPDFSherpa_static STATIC ${LHAPDFSherpa_esources}
                             ${LHAPDF_SOURCES} 
)
target_include_directories(LHAPDFSherpa_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(LHAPDFSherpa_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS LHAPDFSherpa_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(LHAPDFSherpa_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME LHAPDFSherpa)
endif()
if (SHERPA_BUILD_SHARED)
add_library(LHAPDFSherpa SHARED ${LHAPDFSherpa_esources}
                             ${LHAPDF_SOURCES} 
)
target_include_directories(LHAPDFSherpa PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(LHAPDFSherpa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS LHAPDFSherpa DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(LHAPDFSherpa        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME LHAPDFSherpa SOVERSION 0.0.0)
endif()
