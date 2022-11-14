##### Automatically or semiautomaticaly generated, do not edit.#########
#TauPi/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaTauPi.la)
 set(GITTAG  MODEL/TauPi)
set(TauPi_SOURCES  TauPi_LC.C Model.C)
creategitinfo(TauPi ../TauPi/)
list(TRANSFORM TauPi_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../TauPi/")
list(APPEND TauPi_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaTauPi SHARED ${TauPi_SOURCES} )
target_include_directories(SherpaTauPi PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(SherpaTauPi)
install(TARGETS SherpaTauPi DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaTauPi        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaTauPi SOVERSION 0.0.0)
set_target_properties(SherpaTauPi        PROPERTIES DEFINE_SYMBOL "")
