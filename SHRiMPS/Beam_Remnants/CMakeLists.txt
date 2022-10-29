
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Beam_Remnants/Makefile.am
 
 #pkglib_LTLIBRARIES = libShrimpsBeamRemnants.la)
 set(GITTAG  SHRiMPS/Beam_Remnants)
 #include ../../git.make)
set(Beam_Remnants_SOURCES  Continued_PDF.C Hadron_Dissociation.C Remnant_Handler.C Colour_Generator.C)
 #localincdir = $(pkg#includedir)/SHRiMPS/Beam_Remnants)
set(Beam_Remnants_HEADERS Continued_PDF.H Hadron_Dissociation.H Remnant_Handler.H Colour_Generator.H)
creategitinfo(Beam_Remnants ../Beam_Remnants/)
  list(TRANSFORM Beam_Remnants_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Beam_Remnants/")
 
set(Beam_Remnants_esources )
if (SHERPA_BUILD_STATIC)
add_library(ShrimpsBeamRemnants_static STATIC ${ShrimpsBeamRemnants_esources}
                             ${Beam_Remnants_SOURCES} 
)
target_include_directories(ShrimpsBeamRemnants_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ShrimpsBeamRemnants_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ShrimpsBeamRemnants_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ShrimpsBeamRemnants_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME ShrimpsBeamRemnants)
set_target_properties(ShrimpsBeamRemnants_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(ShrimpsBeamRemnants SHARED ${ShrimpsBeamRemnants_esources}
                             ${Beam_Remnants_SOURCES} 
)
target_include_directories(ShrimpsBeamRemnants PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ShrimpsBeamRemnants PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ShrimpsBeamRemnants DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ShrimpsBeamRemnants        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ShrimpsBeamRemnants SOVERSION 0.0.0)
set_target_properties(ShrimpsBeamRemnants        PROPERTIES DEFINE_SYMBOL "")
endif()
