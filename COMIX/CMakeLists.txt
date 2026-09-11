add_subdirectory(Main)
add_subdirectory(Amplitude)
add_subdirectory(Phasespace)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(ComixCore SHARED
    $<TARGET_OBJECTS:Comix> $<TARGET_OBJECTS:ComixAmplitude> $<TARGET_OBJECTS:ComixPhasespace>)
  target_link_libraries(ComixCore PRIVATE Comix ComixAmplitude ComixPhasespace)
  sherpa_handle_shared_library(ComixCore ComixCore)
  add_dependencies(ComixCore ComixMain_Git_Info ComixAmplitude_Git_Info ComixPhasespace_Git_Info)
endif()
