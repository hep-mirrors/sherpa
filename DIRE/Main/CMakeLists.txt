##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libDireMain.la)
 set(GITTAG  DIRE/Main)
set(Main_SOURCES  Jet_Criterion.C Dire.C)
sherpa_create_git_info(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(DireMain SHARED ${Main_SOURCES} )
sherpa_handle_shared_library(DireMain DireMain)
