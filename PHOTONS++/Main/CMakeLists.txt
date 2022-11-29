##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
set(GITTAG  PHOTONS++/Main)

set(Main_SOURCES Photons.C
                 Define_Dipole.C
                 Dipole_Type.C)

set(Main_HEADERS Photons.H
                 Dipole_Type.H
                 Define_Dipole.H)

sherpa_create_git_info(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(PhotonsMain SHARED ${Main_SOURCES} )
sherpa_handle_shared_library(PhotonsMain PhotonsMain)
