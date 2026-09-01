add_subdirectory(Main)
add_subdirectory(NLO)
add_subdirectory(CEEX)
add_subdirectory(Tools)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(YFS SHARED
    $<TARGET_OBJECTS:YFSMain> $<TARGET_OBJECTS:YFSNLO>
    $<TARGET_OBJECTS:YFSCEEX> $<TARGET_OBJECTS:YFSTools>)
  target_link_libraries(YFS PRIVATE YFSMain YFSNLO YFSCEEX YFSTools)
  sherpa_handle_shared_library(YFS YFS)
  add_dependencies(YFS YFSMain_Git_Info YFSNLO_Git_Info CEEX_Git_Info YFSTools_Git_Info)
endif()
