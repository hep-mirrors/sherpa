##### Automatically or semiautomaticaly generated, do not edit.#########
#Tools/Makefile.am
 
set(GITTAG  CSSHOWER++/Tools)

set(Tools_SOURCES Parton.C
                  Singlet.C)

set(Tools_HEADERS Parton.H
                  Singlet.H)

sherpa_create_git_info(Tools)
 
add_library(CSTools SHARED ${Tools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(CSTools CSTools)
