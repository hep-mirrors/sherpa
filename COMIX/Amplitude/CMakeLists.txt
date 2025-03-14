set(ComixAmplitude_SOURCES Amplitude.C)

set(ComixAmplitude_HEADERS Amplitude.H)

sherpa_create_git_info(ComixAmplitude)
 
add_library(ComixAmplitude SHARED ${ComixAmplitude_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ComixAmplitude ComixAmplitude)
