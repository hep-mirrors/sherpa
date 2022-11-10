
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
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
list(APPEND SoftPhysics_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaSoftPhysics SHARED 
                             ${SoftPhysics_SOURCES} 
)
target_include_directories(SherpaSoftPhysics PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(SherpaSoftPhysics)
install(TARGETS SherpaSoftPhysics DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaSoftPhysics        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaSoftPhysics SOVERSION 0.0.0)
set_target_properties(SherpaSoftPhysics        PROPERTIES DEFINE_SYMBOL "")
