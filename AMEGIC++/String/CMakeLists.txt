set(String_SOURCES String_Handler.C
                   String_Generator.C
                   String_Tree.C
                   String_Output.C
                   String_Library.C)

set(String_HEADERS String_Handler.H
                   String_Generator.H
                   String_Tree.H
                   String_Output.H
                   MyString.H
                   Values.H
                   String_Library.H)

sherpa_create_git_info(String ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(String SHARED ${String_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(String String)
target_include_directories(String PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
add_dependencies(String String_Git_Info)
