set(Two2Three_SOURCES ee_3pi.C)

sherpa_create_git_info(Two2Three)

add_library(ExtraXS2_3 SHARED ${Two2Three_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ExtraXS2_3 ExtraXS2_3)
