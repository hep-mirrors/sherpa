##### Automatically or semiautomaticaly generated, do not edit.#########
#Initialization/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaInitialization.la)
 set(GITTAG  SHERPA/Initialization)
set(Initialization_SOURCES  Initialization_Handler.C)
set(Initialization_HEADERS Initialization_Handler.H)
sherpa_create_git_info(Initialization ../Initialization/)
list(TRANSFORM Initialization_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Initialization/")
list(APPEND Initialization_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaInitialization SHARED ${Initialization_SOURCES} )
sherpa_handle_shared_library(SherpaInitialization SherpaInitialization)
