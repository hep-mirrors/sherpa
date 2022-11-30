##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
set(GITTAG  AMISIC/Main)

set(Main_SOURCES  Amisic.C)

set(Main_HEADERS Amisic.H)

sherpa_create_git_info(Main)
 
add_library(Amisic SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(Amisic Amisic)
