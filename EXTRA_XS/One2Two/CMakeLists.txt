##### Automatically or semiautomaticaly generated, do not edit.#########
#One2Two/Makefile.am
 
set(GITTAG  EXTRA_XS/One2Two)
set(One2Two_SOURCES  Comix1to2.C)
set(One2Two_HEADERS Comix1to2.H)
sherpa_create_git_info(One2Two ../One2Two/)
list(TRANSFORM One2Two_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../One2Two/")
list(APPEND One2Two_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ExtraXS1_2 SHARED ${One2Two_SOURCES} )
sherpa_handle_shared_library(ExtraXS1_2 ExtraXS1_2)
