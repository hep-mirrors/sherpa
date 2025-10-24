set(Fusing_SOURCES Fusing_Direct_Hook.C
                   Fusing_Fragmentation_Hook.C)

sherpa_create_git_info(Fusing ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(SherpaFusing SHARED ${Fusing_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaFusing SherpaFusing)
