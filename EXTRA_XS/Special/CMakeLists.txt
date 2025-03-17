set(Special_SOURCES Instanton.C)

sherpa_create_git_info(Special ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(ExtraXSSpecial SHARED ${Special_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ExtraXSSpecial ExtraXSSpecial)
add_dependencies(ExtraXSSpecial Special_Git_Info)
