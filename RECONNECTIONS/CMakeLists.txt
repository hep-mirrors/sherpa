
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Main/Makefile.am
# Makefile.am

########################################################################
#RECONNECTIONS/Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libReconnections.la)
 set(GITTAG  RECONNECTIONS/Main)
 #include ../../git.make)
set(RECONNECTIONS_Main_SOURCES Reconnection_Handler.C Reconnection_Base.C Reconnect_By_Singlet.C Reconnect_Statistical.C Reconnection_Weights.C)
set(RECONNECTIONS_Main_HEADERS Reconnection_Handler.H Reconnection_Base.H Reconnect_By_Singlet.H Reconnect_Statistical.H Reconnection_Weights.H)
creategitinfo(RECONNECTIONS_Main ../RECONNECTIONS/Main/)
  list(TRANSFORM RECONNECTIONS_Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../RECONNECTIONS/Main/")
 
#RECONNECTIONS/Makefile.am
 
 #SUBDIRS = Main)
 
set(RECONNECTIONS_esources )
if (SHERPA_BUILD_STATIC)
add_library(Reconnections_static STATIC ${Reconnections_esources}
                             ${RECONNECTIONS_Main_SOURCES} 
)
target_include_directories(Reconnections_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(Reconnections_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS Reconnections_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(Reconnections_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME Reconnections)
endif()
if (SHERPA_BUILD_SHARED)
add_library(Reconnections SHARED ${Reconnections_esources}
                             ${RECONNECTIONS_Main_SOURCES} 
)
target_include_directories(Reconnections PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(Reconnections PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS Reconnections DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(Reconnections        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME Reconnections SOVERSION 0.0.0)
endif()
