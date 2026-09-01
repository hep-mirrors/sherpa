add_subdirectory(Gauge)
add_subdirectory(Lorentz)
add_subdirectory(Main)
add_subdirectory(Shower)
add_subdirectory(Tools)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(Dire SHARED
    $<TARGET_OBJECTS:DireGauge> $<TARGET_OBJECTS:DireLorentz>
    $<TARGET_OBJECTS:DireMain> $<TARGET_OBJECTS:DireShower> $<TARGET_OBJECTS:DireTools>)
  target_link_libraries(Dire PRIVATE DireGauge DireLorentz DireMain DireShower DireTools)
  sherpa_handle_shared_library(Dire Dire)
  add_dependencies(Dire DireGauge_Git_Info DireLorentz_Git_Info DireMain_Git_Info DireShower_Git_Info DireTools_Git_Info)
endif()
