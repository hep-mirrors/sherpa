#Manually edited

set(GITTAG  AddOns/MCFM)
set(MCFM_SOURCES MCFM_Interface.C)
sherpa_create_git_info(MCFM ../MCFM/)
list(TRANSFORM MCFM_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../MCFM/")
list(APPEND MCFM_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)

add_library(SherpaMCFM SHARED MCFM_Interface.C)
sherpa_handle_shared_library(SherpaMCFM SherpaMCFM)
target_include_directories(SherpaMCFM  PRIVATE ${MCFM_INCLUDE_DIR})
target_link_libraries(SherpaMCFM  PRIVATE ${MCFM_LIBRARY})

#Not sure if this is needed
if (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
  target_compile_options(SherpaMCFM       PRIVATE -fopenmp)
endif()
if (${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
  target_compile_options(SherpaMCFM       PRIVATE -qopenmp)
endif()
