##### Automatically or semiautomaticaly generated, do not edit.#########
#SMEHC/Makefile.am
 
set(GITTAG  MODEL/SMEHC)
set(SMEHC_SOURCES  Effective_Higgs_Coupling.C Model.C SVV_LC.C SVVV_LC.C SVVP_LC.C SVVVV_LC.C)
set(SMEHC_HEADERS Effective_Higgs_Coupling.H)
sherpa_create_git_info(SMEHC ../SMEHC/)
list(TRANSFORM SMEHC_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../SMEHC/")
list(APPEND SMEHC_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaSMEHC SHARED ${SMEHC_SOURCES} )
sherpa_handle_shared_library(SherpaSMEHC SherpaSMEHC)
