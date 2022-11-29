
set(GITTAG  AddOns/MadLoop)
set(MadLoop_SOURCES ML5_Interface.C)
sherpa_create_git_info(MadLoop)
add_library(SherpaML5 SHARED ${MadLoop_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
sherpa_handle_shared_library(SherpaML5 SherpaML5)
