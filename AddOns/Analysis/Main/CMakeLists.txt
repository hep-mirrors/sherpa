##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
set(GITTAG  AddOns/Analysis/Main)

set(Main_SOURCES Analysis_Object.C
                 Primitive_Analysis.C
                 Analysis_Handler.C)

set(Main_HEADERS Analysis_Object.H
                 Primitive_Analysis.H
                 Analysis_Handler.H)

sherpa_create_git_info(Main)
 
add_library(SherpaAnalysis SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaAnalysis SherpaAnalysis)
