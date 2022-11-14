##### Automatically or semiautomaticaly generated, do not edit.#########
#Lorentz/Makefile.am
 
 #pkglib_LTLIBRARIES = libDireLorentz.la)
 set(GITTAG  DIRE/Lorentz)
set(Lorentz_SOURCES  FFV_FF.C FFV_FI.C FFV_IF.C FFV_II.C VVV_FF.C VVV_FI.C VVV_IF.C VVV_II.C FFFF_FF.C FFFF_FI.C FFFF_IF.C FFFF_II.C)
creategitinfo(Lorentz ../Lorentz/)
list(TRANSFORM Lorentz_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Lorentz/")
list(APPEND Lorentz_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(DireLorentz SHARED ${Lorentz_SOURCES} )
target_include_directories(DireLorentz PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(DireLorentz)
install(TARGETS DireLorentz DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DireLorentz        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DireLorentz SOVERSION 0.0.0)
set_target_properties(DireLorentz        PROPERTIES DEFINE_SYMBOL "")
