add_subdirectory(Calculators)
add_subdirectory(Showers)
add_subdirectory(Main)
add_subdirectory(Tools)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(MCatNLO SHARED
    $<TARGET_OBJECTS:MCatNLOCalculators> $<TARGET_OBJECTS:MCatNLOShowers>
    $<TARGET_OBJECTS:MCatNLOMain> $<TARGET_OBJECTS:MCatNLOTools>)
  target_link_libraries(MCatNLO PRIVATE MCatNLOCalculators MCatNLOShowers MCatNLOMain MCatNLOTools)
  sherpa_handle_shared_library(MCatNLO MCatNLO)
  add_dependencies(MCatNLO MCatNLOCalculators_Git_Info MCatNLOShowers_Git_Info MCatNLOMain_Git_Info MCatNLOTools_Git_Info)
endif()
