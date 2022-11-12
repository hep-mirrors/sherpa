
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
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
list(APPEND LHAPDF_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(LHAPDFSherpa SHARED 
                             ${LHAPDF_SOURCES} 
)
target_include_directories(LHAPDFSherpa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(LHAPDFSherpa)
install(TARGETS LHAPDFSherpa DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(LHAPDFSherpa        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME LHAPDFSherpa SOVERSION 0.0.0)
set_target_properties(LHAPDFSherpa        PROPERTIES DEFINE_SYMBOL "")


target_link_libraries(LHAPDFSherpa  PRIVATE ${LHAPDF_LIBRARIES})
target_include_directories(LHAPDFSherpa  PRIVATE ${LHAPDF_INCLUDE_DIRS})
