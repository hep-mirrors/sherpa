
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libHadronsMain.la)
 set(GITTAG  HADRONS++/Main)
 #include ../../git.make)
set(Main_SOURCES  Hadron_Decay_Map.C Hadron_Decay_Table.C Hadron_Decay_Channel.C Mixing_Handler.C Tools.C)
set(Main_HEADERS Hadron_Decay_Map.H Hadron_Decay_Table.H Hadron_Decay_Channel.H Mixing_Handler.H Tools.H)
creategitinfo(Main ../Main/)
  list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
 
set(Main_esources )
if (SHERPA_BUILD_STATIC)
add_library(HadronsMain_static STATIC ${HadronsMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(HadronsMain_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(HadronsMain_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS HadronsMain_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(HadronsMain_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME HadronsMain)
endif()
if (SHERPA_BUILD_SHARED)
add_library(HadronsMain SHARED ${HadronsMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(HadronsMain PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(HadronsMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS HadronsMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(HadronsMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME HadronsMain SOVERSION 0.0.0)
endif()
