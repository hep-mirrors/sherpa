##### Automatically or semiautomaticaly generated, do not edit.#########
#PS_Library/Makefile.am
 
 #pkglib_LTLIBRARIES = libHadronsPSs.la)
 set(GITTAG  HADRONS++/PS_Library)
set(PS_Library_SOURCES  HD_PS_Base.C Three_Body_PSs.C Two_Body_PSs.C Four_Body_PSs.C ResonanceFlavour.C)
set(PS_Library_HEADERS HD_PS_Base.H Three_Body_PSs.H Two_Body_PSs.H Four_Body_PSs.H ResonanceFlavour.H)
sherpa_create_git_info(PS_Library ../PS_Library/)
list(TRANSFORM PS_Library_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../PS_Library/")
list(APPEND PS_Library_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(HadronsPSs SHARED ${PS_Library_SOURCES} )
sherpa_handle_shared_library(HadronsPSs HadronsPSs)
