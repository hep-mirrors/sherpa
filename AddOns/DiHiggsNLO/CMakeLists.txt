##### Automatically or semiautomaticaly generated, do not edit.#########
#DiHiggsNLO/Makefile.am
 
set(GITTAG  AddOns/DiHiggsNLO)

set(DiHiggsNLO_SOURCES DiHiggs_Virtual.C
                       hhgrid.c)

set(DiHiggsNLO_HEADERS DiHiggs_Virtual.H
                       hhgrid.h)

sherpa_create_git_info(DiHiggsNLO)
 
add_library(SherpaDiHiggsNLO SHARED ${DiHiggsNLO_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaDiHiggsNLO SherpaDiHiggsNLO)

target_include_directories(SherpaDiHiggsNLO PRIVATE ${Python_INCLUDE_DIRS})
target_link_libraries(SherpaDiHiggsNLO PRIVATE ${Python_LIBRARIES})
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/creategrid.py DESTINATION ${SHERPA_Python_SITEARCH}/  COMPONENT python)
