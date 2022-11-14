##### Automatically or semiautomaticaly generated, do not edit.#########
#String/Makefile.am
 
 set(GITTAG  AMEGIC++/String)
set(String_SOURCES  String_Handler.C String_Generator.C String_Tree.C String_Output.C MyString.C String_Library.C)
set(String_HEADERS String_Handler.H String_Generator.H String_Tree.H String_Output.H MyString.H Values.H String_Library.H)
sherpa_create_git_info(String ../String/)
list(TRANSFORM String_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../String/")
list(APPEND String_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(String SHARED ${String_SOURCES} )
sherpa_handle_shared_library(String String)
target_include_directories(String PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
