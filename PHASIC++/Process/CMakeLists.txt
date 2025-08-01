set(Process_SOURCES Subprocess_Info.C
                    Process_Info.C
                    Process_Base.C
                    Process_Group.C
                    Single_Process.C
                    MCatNLO_Process.C
                    Virtual_ME2_Base.C
                    Tree_ME2_Base.C
                    Color_Correlated_ME2.C
                    Spin_Color_Correlated_ME2.C
                    Massive_Kernels.C
                    KP_Terms.C
                    ME_Generators.C
                    ME_Generator_Base.C
                    YFS_Process.C)

set(Process_HEADERS Subprocess_Info.H
                    Process_Info.H
                    Process_Base.H
                    Process_Group.H
                    Single_Process.H
                    MCatNLO_Process.H
                    Virtual_ME2_Base.H
                    Tree_ME2_Base.H
                    Color_Correlated_ME2.H
                    Spin_Color_Correlated_ME2.H
                    External_ME_Args.H
                    Massive_Kernels.H
                    KP_Terms.H
                    ME_Generators.H
                    ME_Generator_Base.H
                    YFS_Process.H)

sherpa_create_git_info(Process ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(PhasicProcess SHARED ${Process_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(PhasicProcess PhasicProcess)
add_dependencies(PhasicProcess Process_Git_Info)
if (SHERPA_ENABLE_THREADING)
  target_link_libraries(PhasicProcess PRIVATE Threads::Threads)
endif()
