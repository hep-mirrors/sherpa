##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
set(GITTAG  DIRE/Main)

set(Main_SOURCES Jet_Criterion.C
                 Dire.C)

sherpa_create_git_info(Main)
 
add_library(DireMain SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DireMain DireMain)
