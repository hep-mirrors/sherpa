##### Automatically or semiautomaticaly generated, do not edit.#########
#Tools/Makefile.am
 
 #pkglib_LTLIBRARIES = libMCatNLOTools.la)
 set(GITTAG  MCATNLO/Tools)
set(Tools_SOURCES  Parton.C Singlet.C)
set(Tools_HEADERS Parton.H Singlet.H)
sherpa_create_git_info(Tools ../Tools/)
list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
list(APPEND Tools_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(MCatNLOTools SHARED ${Tools_SOURCES} )
sherpa_handle_shared_library(MCatNLOTools MCatNLOTools)
