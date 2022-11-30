set(GoSam_SOURCES GoSam_Virtual.C
                  GoSam_Interface.C)

set(GoSam_HEADERS GoSam_Virtual.H
                  GoSam_Interface.H)

sherpa_create_git_info(GoSam)
 
add_library(SherpaGoSam SHARED ${GoSam_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaGoSam SherpaGoSam)
