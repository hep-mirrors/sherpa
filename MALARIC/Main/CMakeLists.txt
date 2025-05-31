set(MalaricMain_SOURCES Gamma.C
                 MCatNLO.C)

set(MalaricMain_HEADERS Gamma.H
                 MCatNLO.H)

sherpa_create_git_info(MalaricMain)
 
add_library(MalaricMain SHARED ${MalaricMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(MalaricMain MalaricMain)
