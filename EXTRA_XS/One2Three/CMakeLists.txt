##### Automatically or semiautomaticaly generated, do not edit.#########
#One2Three/Makefile.am
 
 #pkglib_LTLIBRARIES = libExtraXS1_3.la)
 set(GITTAG  EXTRA_XS/One2Three)
set(One2Three_SOURCES  Comix1to3.C)
set(One2Three_HEADERS Comix1to3.H)
sherpa_create_git_info(One2Three ../One2Three/)
list(TRANSFORM One2Three_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../One2Three/")
list(APPEND One2Three_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ExtraXS1_3 SHARED ${One2Three_SOURCES} )
sherpa_handle_shared_library(ExtraXS1_3 ExtraXS1_3)
