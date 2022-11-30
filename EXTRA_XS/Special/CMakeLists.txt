set(Special_SOURCES Instanton.C)

sherpa_create_git_info(Special)
 
add_library(ExtraXSSpecial SHARED ${Special_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ExtraXSSpecial ExtraXSSpecial)
