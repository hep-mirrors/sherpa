set(Tools_SOURCES DEQ_Solver.C
                  Kernels.C
                  MinBias_Parameters.C)

set(Tools_HEADERS Parameter_Structures.H
                  Special_Functions.H
                  DEQ_Solver.H
                  Kernels.H
                  MinBias_Parameters.H)

sherpa_create_git_info(Tools)
 
add_library(ShrimpsTools SHARED ${Tools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ShrimpsTools ShrimpsTools)
