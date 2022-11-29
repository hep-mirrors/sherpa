##### Automatically or semiautomaticaly generated, do not edit.#########
#Process/Makefile.am
 
set(GITTAG  PHASIC++/Process)

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
                    ME_Generator_Base.C)

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
                    ME_Generator_Base.H)

sherpa_create_git_info(Process ../Process/)
list(TRANSFORM Process_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Process/")
list(APPEND Process_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(PhasicProcess SHARED ${Process_SOURCES} )
sherpa_handle_shared_library(PhasicProcess PhasicProcess)
