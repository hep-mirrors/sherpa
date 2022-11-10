
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Main/Makefile.am
# Makefile.am

########################################################################
#RECONNECTIONS/Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libReconnections.la)
 set(GITTAG  RECONNECTIONS/Main)
 #include ../../git.make)
set(RECONNECTIONS_Main_SOURCES  Reconnection_Handler.C Reconnection_Base.C Reconnect_By_Singlet.C Reconnect_Statistical.C Reconnection_Weights.C)
set(RECONNECTIONS_Main_HEADERS Reconnection_Handler.H Reconnection_Base.H Reconnect_By_Singlet.H Reconnect_Statistical.H Reconnection_Weights.H)
creategitinfo(RECONNECTIONS_Main ../RECONNECTIONS/Main/)
list(TRANSFORM RECONNECTIONS_Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../RECONNECTIONS/Main/")
list(APPEND RECONNECTIONS_Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
#RECONNECTIONS/Makefile.am
 
 #SUBDIRS = Main)
 
add_library(Reconnections SHARED 
                             ${RECONNECTIONS_Main_SOURCES} 
)
target_include_directories(Reconnections PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(Reconnections)
install(TARGETS Reconnections DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(Reconnections        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME Reconnections SOVERSION 0.0.0)
set_target_properties(Reconnections        PROPERTIES DEFINE_SYMBOL "")
