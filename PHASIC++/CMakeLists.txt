add_subdirectory(Channels)
add_subdirectory(Decays)
add_subdirectory(Main)
add_subdirectory(Enhance)
add_subdirectory(Process)
add_subdirectory(Scales)
add_subdirectory(Selectors)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(Phasic SHARED
    $<TARGET_OBJECTS:PhasicChannels> $<TARGET_OBJECTS:PhasicDecays>
    $<TARGET_OBJECTS:PhasicEnhance> $<TARGET_OBJECTS:PhasicMain>
    $<TARGET_OBJECTS:PhasicProcess> $<TARGET_OBJECTS:PhasicScales>
    $<TARGET_OBJECTS:PhasicSelectors>)
  target_link_libraries(Phasic PRIVATE PhasicChannels PhasicDecays PhasicEnhance PhasicMain PhasicProcess PhasicScales PhasicSelectors)
  sherpa_handle_shared_library(Phasic Phasic)
  add_dependencies(Phasic Channels_Git_Info PhasicDecays_Git_Info Enhance_Git_Info PhasicMain_Git_Info Process_Git_Info Scales_Git_Info Selectors_Git_Info)
endif()
