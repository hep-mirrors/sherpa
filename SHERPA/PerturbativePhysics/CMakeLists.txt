
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#PerturbativePhysics/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaPerturbativePhysics.la)
 set(GITTAG  SHERPA/PerturbativePhysics)
 #include ../../git.make)
set(PerturbativePhysics_SOURCES  Matrix_Element_Handler.C MI_Handler.C Shower_Handler.C Perturbative_Interface.C Hard_Decay_Handler.C)
set(PerturbativePhysics_HEADERS Matrix_Element_Handler.H MI_Handler.H Shower_Handler.H Perturbative_Interface.H Hard_Decay_Handler.H)
creategitinfo(PerturbativePhysics ../PerturbativePhysics/)
list(TRANSFORM PerturbativePhysics_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../PerturbativePhysics/")
list(APPEND PerturbativePhysics_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaPerturbativePhysics SHARED 
                             ${PerturbativePhysics_SOURCES} 
)
target_include_directories(SherpaPerturbativePhysics PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaPerturbativePhysics DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaPerturbativePhysics        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaPerturbativePhysics SOVERSION 0.0.0)
set_target_properties(SherpaPerturbativePhysics        PROPERTIES DEFINE_SYMBOL "")
