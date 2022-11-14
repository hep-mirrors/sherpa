##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
set(GITTAG  AMISIC/Main)
set(Main_SOURCES  Amisic.C)
set(Main_HEADERS Amisic.H)
sherpa_create_git_info(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(Amisic SHARED ${Main_SOURCES} )
sherpa_handle_shared_library(Amisic Amisic)
