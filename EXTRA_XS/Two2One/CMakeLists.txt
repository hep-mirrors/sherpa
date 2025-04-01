set(Two2One_SOURCES XS_vv2quarkonium.C
)

sherpa_create_git_info(Two2One)
 
add_library(ExtraXS2_1 SHARED ${Two2One_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ExtraXS2_1 ExtraXS2_1)
