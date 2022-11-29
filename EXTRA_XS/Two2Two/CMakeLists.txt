##### Automatically or semiautomaticaly generated, do not edit.#########
#Two2Two/Makefile.am
 
set(GITTAG  EXTRA_XS/Two2Two)

set(Two2Two_SOURCES XS_Drell_Yan.C
                    XS_DiPhoton.C
                    XS_QCD.C
                    ee_Y4S_BB.C
                    gg_yy.C
                    yy_bobo.C
                    DMDM_mumu.C)

sherpa_create_git_info(Two2Two ../Two2Two/)
list(TRANSFORM Two2Two_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Two2Two/")
list(APPEND Two2Two_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ExtraXS2_2 SHARED ${Two2Two_SOURCES} )
sherpa_handle_shared_library(ExtraXS2_2 ExtraXS2_2)
