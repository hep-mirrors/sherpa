##### Automatically or semiautomaticaly generated, do not edit.#########
#SMGold/Makefile.am
 
set(GITTAG  MODEL/SMGold)
set(SMGold_SOURCES  Model.C Identity_2_1.C T_3_2_1.C VSS1.C VSS2.C VVVV2.C VVVV3.C VVVV4.C)
sherpa_create_git_info(SMGold ../SMGold/)
list(TRANSFORM SMGold_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../SMGold/")
list(APPEND SMGold_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaSMGold SHARED ${SMGold_SOURCES} )
sherpa_handle_shared_library(SherpaSMGold SherpaSMGold)
