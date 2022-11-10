
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Tools/Makefile.am
 
 #pkglib_LTLIBRARIES = libDIMTools.la)
 set(GITTAG  DIM/Tools)
 #include ../../git.make)
set(Tools_SOURCES  Color.C Parton.C Splitting.C Amplitude.C)
set(Tools_HEADERS Color.H Parton.H Weight.H Splitting.H Kernel_Key.H Amplitude.H)
creategitinfo(Tools ../Tools/)
list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
list(APPEND Tools_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(DIMTools SHARED 
                             ${Tools_SOURCES} 
)
target_include_directories(DIMTools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(DIMTools)
install(TARGETS DIMTools DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DIMTools        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DIMTools SOVERSION 0.0.0)
set_target_properties(DIMTools        PROPERTIES DEFINE_SYMBOL "")
