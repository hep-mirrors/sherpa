add_subdirectory(Decays)
add_subdirectory(Formation)
add_subdirectory(Main)
add_subdirectory(Tools)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(Ahadic SHARED
    $<TARGET_OBJECTS:AhadicMain> $<TARGET_OBJECTS:AhadicTools>
    $<TARGET_OBJECTS:AhadicFormation> $<TARGET_OBJECTS:AhadicDecays>)
  target_link_libraries(Ahadic PRIVATE AhadicMain AhadicTools AhadicFormation AhadicDecays)
  sherpa_handle_shared_library(Ahadic Ahadic)
  add_dependencies(Ahadic AhadicMain_Git_Info AhadicTools_Git_Info Formation_Git_Info AhadicDecays_Git_Info)
endif()
