add_subdirectory(Main)
add_subdirectory(Perturbative)
add_subdirectory(Tools)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(AmisicCore SHARED
    $<TARGET_OBJECTS:Amisic> $<TARGET_OBJECTS:AmisicPerturbative> $<TARGET_OBJECTS:AmisicTools>)
  target_link_libraries(AmisicCore PRIVATE Amisic AmisicPerturbative AmisicTools)
  sherpa_handle_shared_library(AmisicCore AmisicCore)
  add_dependencies(AmisicCore AmisicMain_Git_Info Perturbative_Git_Info AmisicTools_Git_Info)
endif()
