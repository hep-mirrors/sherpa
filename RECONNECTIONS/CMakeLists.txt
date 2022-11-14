##### Automatically or semiautomaticaly generated, do not edit.#########
#RECONNECTIONS/Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libReconnections.la)
 set(GITTAG  RECONNECTIONS/Main)
set(RECONNECTIONS_Main_SOURCES  Reconnection_Handler.C Reconnection_Base.C Reconnect_By_Singlet.C Reconnect_Statistical.C Reconnection_Weights.C)
set(RECONNECTIONS_Main_HEADERS Reconnection_Handler.H Reconnection_Base.H Reconnect_By_Singlet.H Reconnect_Statistical.H Reconnection_Weights.H)
sherpa_create_git_info(RECONNECTIONS_Main ../RECONNECTIONS/Main/)
list(TRANSFORM RECONNECTIONS_Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../RECONNECTIONS/Main/")
list(APPEND RECONNECTIONS_Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
#RECONNECTIONS/Makefile.am
 
 #SUBDIRS = Main)
 
add_library(Reconnections SHARED ${RECONNECTIONS_Main_SOURCES} )
sherpa_handle_shared_library(Reconnections Reconnections)
