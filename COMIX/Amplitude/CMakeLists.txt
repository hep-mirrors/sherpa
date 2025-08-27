set(ComixAmplitude_SOURCES Amplitude.C)

set(ComixAmplitude_HEADERS Amplitude.H)

sherpa_create_git_info(ComixAmplitude ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(ComixAmplitude SHARED ${ComixAmplitude_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ComixAmplitude ComixAmplitude)
add_dependencies(ComixAmplitude ComixAmplitude_Git_Info)
if (SHERPA_ENABLE_THREADING)
  target_link_libraries(ComixAmplitude PRIVATE Threads::Threads)
endif()
