
set(Recola_SOURCES Recola_Interface.C Recola_Born.C  Recola_Virtual.C)

sherpa_create_git_info(Recola)

add_library(SherpaRecola SHARED ${Recola_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
sherpa_handle_shared_library(SherpaRecola SherpaRecola)
target_include_directories(SherpaRecola  PRIVATE ${RECOLA_INCLUDE_DIR})
target_link_libraries(SherpaRecola  PRIVATE ${RECOLA_LIBRARY})
