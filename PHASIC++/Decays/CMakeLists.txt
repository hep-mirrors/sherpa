
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Decays/Makefile.am
 
 #pkglib_LTLIBRARIES = libPhasicDecays.la)
 set(GITTAG  PHASIC++/Decays)
 #include ../../git.make)
set(Decays_SOURCES  Decay_Channel.C Decay_Table.C Decay_Map.C)
set(Decays_HEADERS Decay_Channel.H Decay_Table.H Decay_Map.H)
creategitinfo(Decays ../Decays/)
  list(TRANSFORM Decays_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Decays/")
 
set(Decays_esources )
if (SHERPA_BUILD_STATIC)
add_library(PhasicDecays_static STATIC ${PhasicDecays_esources}
                             ${Decays_SOURCES} 
)
target_include_directories(PhasicDecays_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PhasicDecays_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhasicDecays_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhasicDecays_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME PhasicDecays)
endif()
if (SHERPA_BUILD_SHARED)
add_library(PhasicDecays SHARED ${PhasicDecays_esources}
                             ${Decays_SOURCES} 
)
target_include_directories(PhasicDecays PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PhasicDecays PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhasicDecays DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhasicDecays        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME PhasicDecays SOVERSION 0.0.0)
endif()
