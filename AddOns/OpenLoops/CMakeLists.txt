set(OpenLoops_SOURCES OpenLoops_Virtual.C
                      OpenLoops_Born.C
                      OpenLoops_Interface.C
                      Color_Correlated_ME2.C
                      Spin_Color_Correlated_ME2.C
                      GGH_Process_Manager.C
                      GGH_KFactor_Setter.C
                      GGH_Decay_KFactor_Setter.C)

set(OpenLoops_HEADERS OpenLoops_Virtual.H
                      OpenLoops_Born.H
                      OpenLoops_Interface.H
                      Color_Correlated_ME2.H
                      Spin_Color_Correlated_ME2.H
                      GGH_Process_Manager.H
                      GGH_KFactor_Setter.H)

sherpa_create_git_info(OpenLoops ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaOpenLoops SHARED ${OpenLoops_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaOpenLoops SherpaOpenLoops)
target_link_libraries(SherpaOpenLoops PRIVATE OpenLoops::OpenLoops )
add_dependencies(SherpaOpenLoops OpenLoops_Git_Info)
