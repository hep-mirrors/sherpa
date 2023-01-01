set(Main_SOURCES Single_Process.C
                 Process_Group.C
                 Simple_XS.C
                 ME_Tools.C
                 ME2_Base.C)

set(Main_HEADERS Single_Process.H
                 Process_Group.H
                 ME_Tools.H
                 ME2_Base.H)

sherpa_create_git_info(Main)
 
add_library(ExtraXS SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ExtraXS ExtraXS)
