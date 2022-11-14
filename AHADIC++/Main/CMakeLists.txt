##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libAhadicMain.la)
 set(GITTAG  AHADIC++/Main)
set(Main_SOURCES  Ahadic.C)
 #localincdir = $(pkg#includedir)/AHADIC++/Main)
set(Main_HEADERS Ahadic.H)
sherpa_create_git_info(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(AhadicMain SHARED ${Main_SOURCES} )
sherpa_handle_shared_library(AhadicMain AhadicMain)
