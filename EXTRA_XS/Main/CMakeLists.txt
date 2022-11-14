##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libExtraXS.la)
 set(GITTAG  EXTRA_XS/Main)
set(Main_SOURCES  Single_Process.C Process_Group.C Simple_XS.C ME_Tools.C ME2_Base.C)
set(Main_HEADERS Single_Process.H Process_Group.H ME_Tools.H ME2_Base.H)
sherpa_create_git_info(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ExtraXS SHARED ${Main_SOURCES} )
sherpa_handle_shared_library(ExtraXS ExtraXS)
