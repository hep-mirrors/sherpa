##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
 set(GITTAG  COMIX/Main)
set(Main_SOURCES  Process_Base.C Process_Group.C Single_Process.C Single_Dipole_Term.C Comix.C)
set(Main_HEADERS Process_Base.H Process_Group.H Single_Process.H Single_Dipole_Term.H)
 #dist_bin_SCRIPTS = init_nlo.sh)
sherpa_create_git_info(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(Comix SHARED ${Main_SOURCES} )
sherpa_handle_shared_library(Comix Comix)
