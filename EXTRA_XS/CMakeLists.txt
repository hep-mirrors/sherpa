add_subdirectory(Main)
add_subdirectory(NLO)
add_subdirectory(One2Three)
add_subdirectory(One2Two)
add_subdirectory(Two2Two)
add_subdirectory(Special)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(ExtraXSCore SHARED
    $<TARGET_OBJECTS:ExtraXS> $<TARGET_OBJECTS:ExtraXSNLO>
    $<TARGET_OBJECTS:ExtraXS1_3> $<TARGET_OBJECTS:ExtraXS1_2>
    $<TARGET_OBJECTS:ExtraXS2_2> $<TARGET_OBJECTS:ExtraXSSpecial>)
  target_link_libraries(ExtraXSCore PRIVATE ExtraXS ExtraXSNLO ExtraXS1_3 ExtraXS1_2 ExtraXS2_2 ExtraXSSpecial)
  sherpa_handle_shared_library(ExtraXSCore ExtraXSCore)
  add_dependencies(ExtraXSCore ExtraXSMain_Git_Info ExtraXSNLO_Git_Info One2Three_Git_Info One2Two_Git_Info Two2Two_Git_Info Special_Git_Info)
endif()
