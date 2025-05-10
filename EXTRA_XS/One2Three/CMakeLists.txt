set(One2Three_SOURCES Comix1to3.C)

set(One2Three_HEADERS Comix1to3.H)
list(APPEND One2Three_SOURCES Comix1to3.C H_to_bbg_at_NLO.C)
list(APPEND One2Three_HEADERS Comix1to3.H H_to_bbg_at_NLO.H)

sherpa_create_git_info(One2Three ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(ExtraXS1_3 SHARED ${One2Three_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ExtraXS1_3 ExtraXS1_3)
add_dependencies(ExtraXS1_3 One2Three_Git_Info)
