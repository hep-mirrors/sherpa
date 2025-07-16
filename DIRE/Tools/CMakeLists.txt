set(DireTools_SOURCES Color.C
                  Parton.C
                  Splitting.C
                  Amplitude.C)

set(DireTools_HEADERS Color.H
                  Parton.H
                  Weight.H
                  Splitting.H
                  Kernel_Key.H
                  Amplitude.H)

sherpa_create_git_info(DireTools ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(DireTools SHARED ${DireTools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DireTools DireTools)
add_dependencies(DireTools DireTools_Git_Info)
