set(DiHiggsNLO_SOURCES DiHiggs_Virtual.C
                       hhgrid.c)

set(DiHiggsNLO_HEADERS DiHiggs_Virtual.H
                       hhgrid.h)

sherpa_create_git_info(DiHiggsNLO ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaDiHiggsNLO SHARED ${DiHiggsNLO_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaDiHiggsNLO SherpaDiHiggsNLO)
add_dependencies(SherpaDiHiggsNLO DiHiggsNLO_Git_Info)

target_link_libraries(SherpaDiHiggsNLO PRIVATE Python::Module)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/creategrid.py DESTINATION ${SHERPA_Python_SITEARCH}/  COMPONENT python)
