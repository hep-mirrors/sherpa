##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
set(GITTAG  PHOTONS++/Main)

set(Main_SOURCES Photons.C
                 Define_Dipole.C
                 Dipole_Type.C)

set(Main_HEADERS Photons.H
                 Dipole_Type.H
                 Define_Dipole.H)

sherpa_create_git_info(Main)
 
add_library(PhotonsMain SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(PhotonsMain PhotonsMain)
