##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libMCatNLOMain.la)
 set(GITTAG  MCATNLO/Main)
set(Main_SOURCES  CS_Cluster_Definitions.C CS_Gamma.C CS_MCatNLO.C)
set(Main_HEADERS CS_Cluster_Definitions.H CS_Gamma.H CS_MCatNLO.H)
sherpa_create_git_info(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(MCatNLOMain SHARED ${Main_SOURCES} )
sherpa_handle_shared_library(MCatNLOMain MCatNLOMain)
