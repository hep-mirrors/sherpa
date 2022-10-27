
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
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
 
set(PerturbativePhysics_esources )
if (SHERPA_BUILD_STATIC)
add_library(SherpaPerturbativePhysics_static STATIC ${SherpaPerturbativePhysics_esources}
                             ${PerturbativePhysics_SOURCES} 
)
target_include_directories(SherpaPerturbativePhysics_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaPerturbativePhysics_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaPerturbativePhysics_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaPerturbativePhysics_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME SherpaPerturbativePhysics)
endif()
if (SHERPA_BUILD_SHARED)
add_library(SherpaPerturbativePhysics SHARED ${SherpaPerturbativePhysics_esources}
                             ${PerturbativePhysics_SOURCES} 
)
target_include_directories(SherpaPerturbativePhysics PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaPerturbativePhysics PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaPerturbativePhysics DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaPerturbativePhysics        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaPerturbativePhysics SOVERSION 0.0.0)
endif()
