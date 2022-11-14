##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libRemnants.la)
 set(GITTAG  REMNANTS/Main)
set(Main_SOURCES  Remnant_Handler.C Remnant_Base.C Hadron_Remnant.C Electron_Remnant.C Photon_Remnant.C No_Remnant.C)
set(Main_HEADERS Remnant_Handler.H Remnant_Base.H Hadron_Remnant.H Electron_Remnant.H Photon_Remnant.H No_Remnant.H)
creategitinfo(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(Remnants SHARED ${Main_SOURCES} )
target_include_directories(Remnants PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(Remnants)
install(TARGETS Remnants DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(Remnants        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME Remnants SOVERSION 0.0.0)
set_target_properties(Remnants        PROPERTIES DEFINE_SYMBOL "")
