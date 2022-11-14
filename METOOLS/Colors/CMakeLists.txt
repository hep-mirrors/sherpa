##### Automatically or semiautomaticaly generated, do not edit.#########
#Colors/Makefile.am
 
 set(GITTAG  METOOLS/Colors)
set(Colors_SOURCES  None_CC.C D_CC.C T_CC.C F_CC.C FF_CC.C TT_CC.C ST_CC.C SF_CC.C)
set(Colors_HEADERS)
sherpa_create_git_info(Colors ../Colors/)
list(TRANSFORM Colors_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Colors/")
list(APPEND Colors_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(METoolsColors SHARED ${Colors_SOURCES} )
sherpa_handle_shared_library(METoolsColors METoolsColors)
