##### Automatically or semiautomaticaly generated, do not edit.#########
#Beam_Remnants/Makefile.am
 
set(GITTAG  SHRiMPS/Beam_Remnants)
set(Beam_Remnants_SOURCES  Continued_PDF.C Hadron_Dissociation.C Remnant_Handler.C Colour_Generator.C)
set(Beam_Remnants_HEADERS Continued_PDF.H Hadron_Dissociation.H Remnant_Handler.H Colour_Generator.H)
sherpa_create_git_info(Beam_Remnants ../Beam_Remnants/)
list(TRANSFORM Beam_Remnants_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Beam_Remnants/")
list(APPEND Beam_Remnants_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ShrimpsBeamRemnants SHARED ${Beam_Remnants_SOURCES} )
sherpa_handle_shared_library(ShrimpsBeamRemnants ShrimpsBeamRemnants)
