##### Automatically or semiautomaticaly generated, do not edit.#########
#Tools/Makefile.am
 
set(GITTAG  DIRE/Tools)

set(Tools_SOURCES Color.C
                  Parton.C
                  Splitting.C
                  Amplitude.C)

set(Tools_HEADERS Color.H
                  Parton.H
                  Weight.H
                  Splitting.H
                  Kernel_Key.H
                  Amplitude.H)

sherpa_create_git_info(Tools)
 
add_library(DireTools SHARED ${Tools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DireTools DireTools)
