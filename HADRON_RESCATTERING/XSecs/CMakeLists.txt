set(HadronRescatteringXSecs_SOURCES
  BaryonBaryon.C
  MesonMeson.C
  NNfits.C
  HPR1R2.C
  PiK.C
  HR_Parameters.C
  BaryonMeson.C
  Helper.C
  HeraParam.C
  ScatteringBase.C
  )

set(HadronRescatteringXSecs_HEADERS
  BaryonBaryon.H
  MesonMeson.H
  NNfits.H
  HPR1R2.H
  PiK.H
  HR_Parameters.H
  BaryonMeson.H
  Helper.H
  HeraParam.H
  ScatteringBase.H
  )

sherpa_create_git_info(HadronRescatteringXSecs ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(HadronRescatteringXSecs SHARED ${HadronRescatteringXSecs_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(HadronRescatteringXSecs HadronRescatteringXSecs)
add_dependencies(HadronRescatteringXSecs HadronRescatteringXSecs_Git_Info)
