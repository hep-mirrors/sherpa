set(Griffin_SOURCES Griffin_Interface.C Griffin_Virtual.C Griffin_Born.C)
set(Griffin_HEADERS Griffin_Interface.H Griffin_Virtual.H Griffin_Born.H)

sherpa_create_git_info(Griffin ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(SherpaGriffin SHARED ${Griffin_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
sherpa_handle_shared_library(SherpaGriffin SherpaGriffin)
target_include_directories(SherpaGriffin  PRIVATE ${Griffin_INCLUDE_DIR})
if(NOT (CMAKE_SYSTEM_NAME MATCHES Darwin))
	target_link_libraries(SherpaGriffin  PRIVATE ${GRIFFIN_PREFIX}/libgriffin.so)
else()
	target_link_libraries(SherpaGriffin  PRIVATE ${GRIFFIN_PREFIX}/libgriffin.dylib)
endif()