set(ExtraXSMain_SOURCES Single_Process.C
                 Process_Group.C
                 Simple_XS.C
                 ME_Tools.C
                 ME2_Base.C)

set(ExtraXSMain_HEADERS Single_Process.H
                 Process_Group.H
                 ME_Tools.H
                 ME2_Base.H)

sherpa_create_git_info(ExtraXSMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(ExtraXS SHARED ${ExtraXSMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ExtraXS ExtraXS)
add_dependencies(ExtraXS ExtraXSMain_Git_Info)
