set(CEEX_SOURCES Ceex_Base.C
                 Ceex_Spinors.C
                 Ceex_Born.C
                 Ceex_Virtual.C
                 Ceex_Real.C
                 Ceex_Partitions.C
                 Ceex_Checks.C
                 Ceex_Beta2.C)

set(CEEX_HEADERS Ceex_Base.H)


sherpa_create_git_info(CEEX ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(YFSCEEX ${SHERPA_SUBLIB_TYPE} ${CEEX_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(YFSCEEX YFSCEEX)
add_dependencies(YFSCEEX CEEX_Git_Info)
