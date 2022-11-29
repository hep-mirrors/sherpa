#Manually edited

set(GITTAG  AddOns/Recola)
set(Recola_SOURCES Recola_Interface.C Recola_Born.C  Recola_Virtual.C)

sherpa_create_git_info(Recola ../Recola/)
list(TRANSFORM Recola_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Recola/")
list(APPEND Recola_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)

add_library(SherpaRecola SHARED ${Recola_SOURCES})
sherpa_handle_shared_library(SherpaRecola SherpaRecola)
target_include_directories(SherpaRecola  PRIVATE ${RECOLA_INCLUDE_DIR})
target_link_libraries(SherpaRecola  PRIVATE ${RECOLA_LIBRARY})
