set(SMGold_SOURCES Model.C
                   Identity_2_1.C
                   T_3_2_1.C
                   VSS1.C
                   VSS2.C
                   VVVV2.C
                   VVVV3.C
                   VVVV4.C)

sherpa_create_git_info(SMGold ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaSMGold SHARED ${SMGold_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaSMGold SherpaSMGold)
add_dependencies(SherpaSMGold SMGold_Git_Info)
