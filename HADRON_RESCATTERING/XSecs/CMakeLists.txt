set(HadronRescatteringXSecs_SOURCES
  BaryonBaryon.C
  NNfits.C
  HPR1R2.C)

set(HadronRescatteringXSecs_HEADERS
  BaryonBaryon.H
  NNfits.H
  HPR1R2.H)

sherpa_create_git_info(HadronRescatteringXSecs ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(HadronRescatteringXSecs SHARED ${HadronRescatteringXSecs_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(HadronRescatteringXSecs HadronRescatteringXSecs)
add_dependencies(HadronRescatteringXSecs HadronRescatteringXSecs_Git_Info)
