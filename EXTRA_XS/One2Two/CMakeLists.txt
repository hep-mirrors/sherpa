set(One2Two_SOURCES Comix1to2.C)

set(One2Two_HEADERS Comix1to2.H)

sherpa_create_git_info(One2Two ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(ExtraXS1_2 SHARED ${One2Two_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ExtraXS1_2 ExtraXS1_2)
add_dependencies(ExtraXS1_2 One2Two_Git_Info)
