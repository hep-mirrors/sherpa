add_subdirectory(Main)
add_subdirectory(Tools)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(RemnantsCore SHARED
    $<TARGET_OBJECTS:Remnants> $<TARGET_OBJECTS:RemnantsTools>)
  target_link_libraries(RemnantsCore PRIVATE Remnants RemnantsTools)
  sherpa_handle_shared_library(RemnantsCore RemnantsCore)
  add_dependencies(RemnantsCore RemnantsMain_Git_Info RemnantsTools_Git_Info)
endif()
