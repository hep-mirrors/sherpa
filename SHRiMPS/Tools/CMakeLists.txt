set(ShrimpsTools_SOURCES DEQ_Solver.C
                  Kernels.C
                  MinBias_Parameters.C)

set(ShrimpsTools_HEADERS Parameter_Structures.H
                  DEQ_Solver.H
                  Kernels.H
                  MinBias_Parameters.H)

sherpa_create_git_info(ShrimpsTools ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(ShrimpsTools SHARED ${ShrimpsTools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ShrimpsTools ShrimpsTools)
add_dependencies(ShrimpsTools ShrimpsTools_Git_Info)
