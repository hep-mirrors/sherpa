##### Automatically or semiautomaticaly generated, do not edit.#########
#Decays/Makefile.am
 
 #pkglib_LTLIBRARIES = libPhasicDecays.la)
 set(GITTAG  PHASIC++/Decays)
set(Decays_SOURCES  Decay_Channel.C Decay_Table.C Decay_Map.C)
set(Decays_HEADERS Decay_Channel.H Decay_Table.H Decay_Map.H)
sherpa_create_git_info(Decays ../Decays/)
list(TRANSFORM Decays_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Decays/")
list(APPEND Decays_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(PhasicDecays SHARED ${Decays_SOURCES} )
sherpa_handle_shared_library(PhasicDecays PhasicDecays)
