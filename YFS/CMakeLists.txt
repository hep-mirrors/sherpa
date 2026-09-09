add_subdirectory(Main)
add_subdirectory(NLO)
add_subdirectory(CEEX)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(YFSCore SHARED
    $<TARGET_OBJECTS:YFSMain> $<TARGET_OBJECTS:YFSNLO>
    $<TARGET_OBJECTS:YFSCEEX>)
  target_link_libraries(YFSCore PRIVATE YFSMain YFSNLO YFSCEEX)
  sherpa_handle_shared_library(YFSCore YFSCore)
  add_dependencies(YFSCore YFSMain_Git_Info YFSNLO_Git_Info CEEX_Git_Info)
endif()
