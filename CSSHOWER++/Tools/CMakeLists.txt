
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Tools/Makefile.am
 
 #pkglib_LTLIBRARIES = libCSTools.la)
 set(GITTAG  CSSHOWER++/Tools)
 #include ../../git.make)
set(Tools_SOURCES  Parton.C Singlet.C)
set(Tools_HEADERS Parton.H Singlet.H)
creategitinfo(Tools ../Tools/)
list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
list(APPEND Tools_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(CSTools SHARED 
                             ${Tools_SOURCES} 
)
target_include_directories(CSTools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(CSTools)
install(TARGETS CSTools DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(CSTools        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME CSTools SOVERSION 0.0.0)
set_target_properties(CSTools        PROPERTIES DEFINE_SYMBOL "")
