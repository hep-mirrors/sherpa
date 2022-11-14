##### Automatically or semiautomaticaly generated, do not edit.#########
#Tools/Makefile.am
 
 #pkglib_LTLIBRARIES = libMCatNLOTools.la)
 set(GITTAG  MCATNLO/Tools)
set(Tools_SOURCES  Parton.C Singlet.C)
set(Tools_HEADERS Parton.H Singlet.H)
creategitinfo(Tools ../Tools/)
list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
list(APPEND Tools_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(MCatNLOTools SHARED ${Tools_SOURCES} )
target_include_directories(MCatNLOTools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(MCatNLOTools)
install(TARGETS MCatNLOTools DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(MCatNLOTools        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME MCatNLOTools SOVERSION 0.0.0)
set_target_properties(MCatNLOTools        PROPERTIES DEFINE_SYMBOL "")
