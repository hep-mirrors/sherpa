##### Automatically or semiautomaticaly generated, do not edit.#########
#Beam_Remnants/Makefile.am
 
 #pkglib_LTLIBRARIES = libShrimpsBeamRemnants.la)
 set(GITTAG  SHRiMPS/Beam_Remnants)
set(Beam_Remnants_SOURCES  Continued_PDF.C Hadron_Dissociation.C Remnant_Handler.C Colour_Generator.C)
 #localincdir = $(pkg#includedir)/SHRiMPS/Beam_Remnants)
set(Beam_Remnants_HEADERS Continued_PDF.H Hadron_Dissociation.H Remnant_Handler.H Colour_Generator.H)
creategitinfo(Beam_Remnants ../Beam_Remnants/)
list(TRANSFORM Beam_Remnants_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Beam_Remnants/")
list(APPEND Beam_Remnants_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ShrimpsBeamRemnants SHARED ${Beam_Remnants_SOURCES} )
target_include_directories(ShrimpsBeamRemnants PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(ShrimpsBeamRemnants)
install(TARGETS ShrimpsBeamRemnants DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ShrimpsBeamRemnants        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ShrimpsBeamRemnants SOVERSION 0.0.0)
set_target_properties(ShrimpsBeamRemnants        PROPERTIES DEFINE_SYMBOL "")
