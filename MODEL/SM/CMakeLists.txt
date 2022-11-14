##### Automatically or semiautomaticaly generated, do not edit.#########
#SM/Makefile.am
 
 set(GITTAG  MODEL/SM)
set(SM_SOURCES  Model.C)
sherpa_create_git_info(SM ../SM/)
list(TRANSFORM SM_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../SM/")
list(APPEND SM_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaSM SHARED ${SM_SOURCES} )
sherpa_handle_shared_library(SherpaSM SherpaSM)
