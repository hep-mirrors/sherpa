set(HadronRescatteringMain_SOURCES
  Hadron_Rescatterings.C)

set(HadronRescatteringMain_HEADERS
  Hadron_Rescatterings.H)

sherpa_create_git_info(HadronRescatteringMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(HadronRescattering SHARED ${HadronRescatteringMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(HadronRescattering HadronRescattering)
add_dependencies(HadronRescattering HadronRescatteringMain_Git_Info)
