
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Formation/Makefile.am
 
 #pkglib_LTLIBRARIES = libAhadicFormation.la)
 set(GITTAG  AHADIC++/Formation)
 #include ../../git.make)
set(Formation_SOURCES  Beam_Particles_Shifter.C Singlet_Former.C Singlet_Checker.C Gluon_Decayer.C Gluon_Splitter.C Trivial_Splitter.C)
 #localincdir = $(pkg#includedir)/AHADIC++/Formation)
set(Formation_HEADERS Beam_Particles_Shifter.H Singlet_Former.H Singlet_Checker.H Gluon_Decayer.H Gluon_Splitter.H Trivial_Splitter.H)
creategitinfo(Formation ../Formation/)
  list(TRANSFORM Formation_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Formation/")
 
set(Formation_esources )
if (SHERPA_BUILD_STATIC)
add_library(AhadicFormation_static STATIC ${AhadicFormation_esources}
                             ${Formation_SOURCES} 
)
target_include_directories(AhadicFormation_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(AhadicFormation_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS AhadicFormation_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(AhadicFormation_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME AhadicFormation)
set_target_properties(AhadicFormation_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(AhadicFormation SHARED ${AhadicFormation_esources}
                             ${Formation_SOURCES} 
)
target_include_directories(AhadicFormation PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(AhadicFormation PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS AhadicFormation DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(AhadicFormation        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME AhadicFormation SOVERSION 0.0.0)
set_target_properties(AhadicFormation        PROPERTIES DEFINE_SYMBOL "")
endif()
