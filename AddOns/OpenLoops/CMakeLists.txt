#Manually edited
#######  Automatically or semiautomaticaly generated, do not edit.
#OpenLoops/Makefile.am
 
set(GITTAG  AddOns/OpenLoops)
set(OpenLoops_SOURCES  OpenLoops_Virtual.C OpenLoops_Born.C OpenLoops_Interface.C Color_Correlated_ME2.C Spin_Color_Correlated_ME2.C GGH_Process_Manager.C GGH_KFactor_Setter.C GGH_Decay_KFactor_Setter.C)
set(OpenLoops_HEADERS OpenLoops_Virtual.H OpenLoops_Born.H OpenLoops_Interface.H Color_Correlated_ME2.H Spin_Color_Correlated_ME2.H GGH_Process_Manager.H GGH_KFactor_Setter.H)
 if (SHERPA_SUPPORT_OPENLOOPS)
set(OpenLoops_SOURCES ${OpenLoops_SOURCES})
 endif()
sherpa_create_git_info(OpenLoops ../OpenLoops/)
list(TRANSFORM OpenLoops_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../OpenLoops/")
list(APPEND OpenLoops_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaOpenLoops SHARED ${OpenLoops_SOURCES} )
target_link_libraries(SherpaOpenLoops PRIVATE ${OPENLOOPS_LIBRARY} )
sherpa_handle_shared_library(SherpaOpenLoops SherpaOpenLoops)
