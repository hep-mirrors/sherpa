##### Automatically or semiautomaticaly generated, do not edit.#########
#SAL/Makefile.am
 
 #pkglib_LTLIBRARIES = libSALSherpa.la)
 set(GITTAG  PDF/SAL)
set(SAL_SOURCES  sal.cpp SALph_CPP_Interface.C)
set(SAL_HEADERS SALph_CPP_Interface.H sal.h math_ut.h)
#nobase_#dist_pkgdata_DATA = SALGrid/SAL.dat
creategitinfo(SAL ../SAL/)
list(TRANSFORM SAL_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../SAL/")
list(APPEND SAL_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SALSherpa SHARED ${SAL_SOURCES} )
target_include_directories(SALSherpa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(SALSherpa)
install(TARGETS SALSherpa DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SALSherpa        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SALSherpa SOVERSION 0.0.0)
set_target_properties(SALSherpa        PROPERTIES DEFINE_SYMBOL "")
