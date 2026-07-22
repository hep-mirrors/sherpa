add_subdirectory(Colors)
add_subdirectory(Currents)
add_subdirectory(Explicit)
add_subdirectory(Loops)
add_subdirectory(Main)
add_subdirectory(SpinCorrelations)
add_subdirectory(Vertices)
add_subdirectory(FormFactors)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(METools SHARED
    $<TARGET_OBJECTS:METoolsColors> $<TARGET_OBJECTS:METoolsCurrents>
    $<TARGET_OBJECTS:METoolsExplicit> $<TARGET_OBJECTS:METoolsLoops>
    $<TARGET_OBJECTS:METoolsMain> $<TARGET_OBJECTS:METoolsSpinCorrelations>
    $<TARGET_OBJECTS:METoolsVertices> $<TARGET_OBJECTS:METoolsFormFactors>)
  target_link_libraries(METools PRIVATE METoolsColors METoolsCurrents METoolsExplicit METoolsLoops METoolsMain METoolsSpinCorrelations METoolsVertices METoolsFormFactors)
  sherpa_handle_shared_library(METools METools)
  add_dependencies(METools Colors_Git_Info Currents_Git_Info Explicit_Git_Info Loops_Git_Info METoolsMain_Git_Info SpinCorrelations_Git_Info Vertices_Git_Info)
endif()
