set(DIMTools_SOURCES Color.C
                  Parton.C
                  Splitting.C
                  Amplitude.C)

set(DIMTools_HEADERS Color.H
                  Parton.H
                  Weight.H
                  Splitting.H
                  Kernel_Key.H
                  Amplitude.H)

sherpa_create_git_info(DIMTools ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(DIMTools SHARED ${DIMTools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DIMTools DIMTools)
add_dependencies(DIMTools DIMTools_Git_Info)
