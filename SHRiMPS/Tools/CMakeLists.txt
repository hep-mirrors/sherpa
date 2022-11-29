##### Automatically or semiautomaticaly generated, do not edit.#########
#Tools/Makefile.am
 
set(GITTAG  SHRiMPS/Tools)

set(Tools_SOURCES DEQ_Solver.C
                  Kernels.C
                  MinBias_Parameters.C)

set(Tools_HEADERS Parameter_Structures.H
                  Special_Functions.H
                  DEQ_Solver.H
                  Kernels.H
                  MinBias_Parameters.H)

sherpa_create_git_info(Tools ../Tools/)
list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
list(APPEND Tools_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ShrimpsTools SHARED ${Tools_SOURCES} )
sherpa_handle_shared_library(ShrimpsTools ShrimpsTools)
