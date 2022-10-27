
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#GRV/Makefile.am
 
 #pkglib_LTLIBRARIES = libGRVSherpa.la)
 set(GITTAG  PDF/GRV)
 #include ../../git.make)
set(GRV_SOURCES  grvphoton.f GRVph_Fortran_Interface.C)
set(GRV_HEADERS GRVph_Fortran_Interface.H)
creategitinfo(GRV ../GRV/)
  list(TRANSFORM GRV_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../GRV/")
 
set(GRV_esources )
if (SHERPA_BUILD_STATIC)
add_library(GRVSherpa_static STATIC ${GRVSherpa_esources}
                             ${GRV_SOURCES} 
)
target_include_directories(GRVSherpa_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(GRVSherpa_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS GRVSherpa_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(GRVSherpa_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME GRVSherpa)
endif()
if (SHERPA_BUILD_SHARED)
add_library(GRVSherpa SHARED ${GRVSherpa_esources}
                             ${GRV_SOURCES} 
)
target_include_directories(GRVSherpa PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(GRVSherpa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS GRVSherpa DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(GRVSherpa        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME GRVSherpa SOVERSION 0.0.0)
endif()
