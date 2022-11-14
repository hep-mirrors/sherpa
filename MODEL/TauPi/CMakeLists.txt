##### Automatically or semiautomaticaly generated, do not edit.#########
#TauPi/Makefile.am
 
set(GITTAG  MODEL/TauPi)
set(TauPi_SOURCES  TauPi_LC.C Model.C)
sherpa_create_git_info(TauPi ../TauPi/)
list(TRANSFORM TauPi_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../TauPi/")
list(APPEND TauPi_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaTauPi SHARED ${TauPi_SOURCES} )
sherpa_handle_shared_library(SherpaTauPi SherpaTauPi)
