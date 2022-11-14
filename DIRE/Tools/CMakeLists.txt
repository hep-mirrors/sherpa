##### Automatically or semiautomaticaly generated, do not edit.#########
#Tools/Makefile.am
 
 #pkglib_LTLIBRARIES = libDireTools.la)
 set(GITTAG  DIRE/Tools)
set(Tools_SOURCES  Color.C Parton.C Splitting.C Amplitude.C)
set(Tools_HEADERS Color.H Parton.H Weight.H Splitting.H Kernel_Key.H Amplitude.H)
creategitinfo(Tools ../Tools/)
list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
list(APPEND Tools_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(DireTools SHARED ${Tools_SOURCES} )
target_include_directories(DireTools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(DireTools)
install(TARGETS DireTools DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DireTools        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DireTools SOVERSION 0.0.0)
set_target_properties(DireTools        PROPERTIES DEFINE_SYMBOL "")
