##### Automatically or semiautomaticaly generated, do not edit.#########
#Special/Makefile.am
 
set(GITTAG  EXTRA_XS/Special)
set(Special_SOURCES  Instanton.C)
sherpa_create_git_info(Special ../Special/)
list(TRANSFORM Special_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Special/")
list(APPEND Special_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ExtraXSSpecial SHARED ${Special_SOURCES} )
sherpa_handle_shared_library(ExtraXSSpecial ExtraXSSpecial)
