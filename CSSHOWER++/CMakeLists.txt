add_subdirectory(Main)
add_subdirectory(Calculators)
add_subdirectory(Showers)
add_subdirectory(Tools)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(CSShower SHARED
    $<TARGET_OBJECTS:CSMain> $<TARGET_OBJECTS:CSCalculators>
    $<TARGET_OBJECTS:CSShowers> $<TARGET_OBJECTS:CSTools>)
  target_link_libraries(CSShower PRIVATE CSMain CSCalculators CSShowers CSTools)
  sherpa_handle_shared_library(CSShower CSShower)
  add_dependencies(CSShower CSSMain_Git_Info CSSCalculators_Git_Info CSSShowers_Git_Info CSSTools_Git_Info)
endif()
