set(CSSTools_SOURCES Parton.C
                  Singlet.C)

set(CSSTools_HEADERS Parton.H
                  Singlet.H)

sherpa_create_git_info(CSSTools)
 
add_library(CSTools SHARED ${CSSTools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(CSTools CSTools)
