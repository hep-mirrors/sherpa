
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#SoftPhysics/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaSoftPhysics.la)
 set(GITTAG  SHERPA/SoftPhysics)
 #include ../../git.make)
set(SoftPhysics_SOURCES  Hadron_Init.C Beam_Remnant_Handler.C Colour_Reconnection_Handler.C Soft_Collision_Handler.C Hadron_Decay_Handler.C Lund_Decay_Handler.C Soft_Photon_Handler.C Singlet_Sorter.C Resonance_Finder.C)
set(SoftPhysics_HEADERS Hadron_Init.H Beam_Remnant_Handler.H Colour_Reconnection_Handler.H Soft_Collision_Handler.H Hadron_Decay_Handler.H Lund_Decay_Handler.H Soft_Photon_Handler.H Singlet_Sorter.H Resonance_Finder.H)
creategitinfo(SoftPhysics ../SoftPhysics/)
  list(TRANSFORM SoftPhysics_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../SoftPhysics/")
 
set(SoftPhysics_esources )
if (SHERPA_BUILD_STATIC)
add_library(SherpaSoftPhysics_static STATIC ${SherpaSoftPhysics_esources}
                             ${SoftPhysics_SOURCES} 
)
target_include_directories(SherpaSoftPhysics_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaSoftPhysics_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaSoftPhysics_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaSoftPhysics_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME SherpaSoftPhysics)
set_target_properties(SherpaSoftPhysics_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(SherpaSoftPhysics SHARED ${SherpaSoftPhysics_esources}
                             ${SoftPhysics_SOURCES} 
)
target_include_directories(SherpaSoftPhysics PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaSoftPhysics PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaSoftPhysics DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaSoftPhysics        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaSoftPhysics SOVERSION 0.0.0)
set_target_properties(SherpaSoftPhysics        PROPERTIES DEFINE_SYMBOL "")
endif()
