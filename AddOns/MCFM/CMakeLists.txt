
set(MCFM_SOURCES MCFM_Interface.C)
sherpa_create_git_info(MCFM ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(SherpaMCFM SHARED ${MCFM_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
sherpa_handle_shared_library(SherpaMCFM SherpaMCFM)
target_link_libraries(SherpaMCFM  PRIVATE MCFM::MCFM)
add_dependencies(SherpaMCFM MCFM_Git_Info)

if (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
  target_compile_options(SherpaMCFM       PRIVATE -fopenmp)
endif()
if (${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
  target_compile_options(SherpaMCFM       PRIVATE -qopenmp)
endif()
