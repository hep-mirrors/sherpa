##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libDIMMain.la)
 set(GITTAG  DIM/Main)
set(Main_SOURCES  Gamma.C MCatNLO.C)
set(Main_HEADERS Gamma.H MCatNLO.H)
sherpa_create_git_info(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(DIMMain SHARED ${Main_SOURCES} )
sherpa_handle_shared_library(DIMMain DIMMain)
