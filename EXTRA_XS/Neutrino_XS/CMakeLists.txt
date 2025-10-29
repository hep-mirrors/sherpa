set(Neutrino_XS_SOURCES 
  vn_ep.C
)

sherpa_create_git_info(Neutrino_XS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(ExtraXS2_2 SHARED ${Neutrino_XS_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ExtraXS2_2 ExtraXS2_2) #TODO: what/where is ExtraXS2_2?
add_dependencies(ExtraXS2_2 Neutrino_XS_Git_Info)
