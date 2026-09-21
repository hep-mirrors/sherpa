add_subdirectory(Gauge)
add_subdirectory(Lorentz)
add_subdirectory(Main)
add_subdirectory(Shower)
add_subdirectory(Tools)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(DIM SHARED
    $<TARGET_OBJECTS:DIMGauge> $<TARGET_OBJECTS:DIMLorentz>
    $<TARGET_OBJECTS:DIMMain> $<TARGET_OBJECTS:DIMShower> $<TARGET_OBJECTS:DIMTools>)
  target_link_libraries(DIM PRIVATE DIMGauge DIMLorentz DIMMain DIMShower DIMTools)
  sherpa_handle_shared_library(DIM DIM)
  add_dependencies(DIM DIMGauge_Git_Info DIMLorentz_Git_Info DIMMain_Git_Info DIMShower_Git_Info DIMTools_Git_Info)
endif()
