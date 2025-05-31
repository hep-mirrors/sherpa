set(AlaricTools_SOURCES Color.C
                  Parton.C
                  Splitting.C
                  Amplitude.C)

set(AlaricTools_HEADERS Color.H
                  Parton.H
                  Weight.H
                  Splitting.H
                  Kernel_Key.H
                  Amplitude.H)

sherpa_create_git_info(AlaricTools ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(AlaricTools SHARED ${AlaricTools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(AlaricTools AlaricTools)
