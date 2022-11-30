##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
set(GITTAG  SHERPA/Main)

set(Main_SOURCES Sherpa.C
                 Filter.C)
set(Main_HEADERS Sherpa.H
                 Filter.H)
sherpa_create_git_info(Main)
 
add_library(SherpaMain SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaMain SherpaMain)
