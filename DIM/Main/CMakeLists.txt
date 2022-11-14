##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libDIMMain.la)
 set(GITTAG  DIM/Main)
set(Main_SOURCES  Gamma.C MCatNLO.C)
set(Main_HEADERS Gamma.H MCatNLO.H)
creategitinfo(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(DIMMain SHARED ${Main_SOURCES} )
target_include_directories(DIMMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(DIMMain)
install(TARGETS DIMMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DIMMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DIMMain SOVERSION 0.0.0)
set_target_properties(DIMMain        PROPERTIES DEFINE_SYMBOL "")
