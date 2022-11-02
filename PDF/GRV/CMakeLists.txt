
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
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
list(APPEND GRV_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(GRVSherpa SHARED 
                             ${GRV_SOURCES} 
)
target_include_directories(GRVSherpa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS GRVSherpa DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(GRVSherpa        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME GRVSherpa SOVERSION 0.0.0)
set_target_properties(GRVSherpa        PROPERTIES DEFINE_SYMBOL "")
