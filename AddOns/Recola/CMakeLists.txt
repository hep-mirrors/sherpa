
set(Recola_SOURCES Recola_Interface.C Recola_Born.C  Recola_Virtual.C Recola_Virtual_VBS_WpWp.C Recola_Virtual_SChannel_WpWp.C)
set(Recola_HEADERS Recola_Interface.H Recola_Born.H  Recola_Virtual.H Recola_Virtual_VBS_WpWp.H Recola_Virtual_SChannel_WpWp.H)

sherpa_create_git_info(Recola ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(SherpaRecola SHARED ${Recola_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
sherpa_handle_shared_library(SherpaRecola SherpaRecola)
target_link_libraries(SherpaRecola  PRIVATE recola::recola)
add_dependencies(SherpaRecola Recola_Git_Info)
