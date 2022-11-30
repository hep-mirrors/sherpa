##### Automatically or semiautomaticaly generated, do not edit.#########
#Decays/Makefile.am
 
set(GITTAG  PHASIC++/Decays)

set(Decays_SOURCES Decay_Channel.C
                   Decay_Table.C
                   Decay_Map.C)

set(Decays_HEADERS Decay_Channel.H
                   Decay_Table.H
                   Decay_Map.H)

sherpa_create_git_info(Decays)
 
add_library(PhasicDecays SHARED ${Decays_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(PhasicDecays PhasicDecays)
