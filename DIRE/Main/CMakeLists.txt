##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libDireMain.la)
 set(GITTAG  DIRE/Main)
set(Main_SOURCES  Jet_Criterion.C Dire.C)
creategitinfo(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(DireMain SHARED ${Main_SOURCES} )
target_include_directories(DireMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(DireMain)
install(TARGETS DireMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DireMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DireMain SOVERSION 0.0.0)
set_target_properties(DireMain        PROPERTIES DEFINE_SYMBOL "")
