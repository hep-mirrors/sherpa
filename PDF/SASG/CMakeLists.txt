
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
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
list(APPEND SASG_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SASGSherpa SHARED 
                             ${SASG_SOURCES} 
)
target_include_directories(SASGSherpa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(SASGSherpa)
install(TARGETS SASGSherpa DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SASGSherpa        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SASGSherpa SOVERSION 0.0.0)
set_target_properties(SASGSherpa        PROPERTIES DEFINE_SYMBOL "")
