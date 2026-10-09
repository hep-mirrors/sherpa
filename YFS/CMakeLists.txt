add_subdirectory(Main)
add_subdirectory(NLO)
add_subdirectory(ACRAIC)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(YFSCore SHARED
    $<TARGET_OBJECTS:YFSMain> $<TARGET_OBJECTS:YFSNLO>
    $<TARGET_OBJECTS:YFSACRAIC>)
  target_link_libraries(YFSCore PRIVATE YFSMain YFSNLO YFSACRAIC)
  sherpa_handle_shared_library(YFSCore YFSCore)
  add_dependencies(YFSCore YFSMain_Git_Info YFSNLO_Git_Info ACRAIC_Git_Info)
endif()
