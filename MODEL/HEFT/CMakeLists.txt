##### Automatically or semiautomaticaly generated, do not edit.#########
#HEFT/Makefile.am
 
set(GITTAG  MODEL/HEFT)

set(HEFT_SOURCES Effective_Higgs_Coupling.C
                 Model.C
                 SVV_LC.C
                 SVVV_LC.C
                 SVVP_LC.C
                 SVVVV_LC.C)

set(HEFT_HEADERS Effective_Higgs_Coupling.H)

sherpa_create_git_info(HEFT)
 
add_library(SherpaHEFT SHARED ${HEFT_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaHEFT SherpaHEFT)
