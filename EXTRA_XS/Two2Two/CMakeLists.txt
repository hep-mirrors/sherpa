set(Two2Two_SOURCES XS_Drell_Yan.C
                    XS_DiPhoton.C
                    XS_QCD.C
                    XS_HadronsQED.C
                    ee_Y4S_BB.C
                    gg_yy.C
                    DMDM_mumu.C)

sherpa_create_git_info(Two2Two ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(ExtraXS2_2 SHARED ${Two2Two_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ExtraXS2_2 ExtraXS2_2)
add_dependencies(ExtraXS2_2 Two2Two_Git_Info)
