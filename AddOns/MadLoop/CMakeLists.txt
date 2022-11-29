
set(GITTAG  AddOns/MadLoop)
set(MadLoop_SOURCES ML5_Interface.C)
sherpa_create_git_info(MadLoop ../MadLoop/)
list(TRANSFORM EWSud_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../MadLoop/")
list(APPEND MadLoop_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)

add_library(SherpaML5 SHARED ${MadLoop_SOURCES})
sherpa_handle_shared_library(SherpaML5 SherpaML5)
