add_subdirectory(Math)
add_subdirectory(Org)
add_subdirectory(Phys)
add_subdirectory(YAML)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(ATOOLSCore SHARED
    $<TARGET_OBJECTS:ToolsMath> $<TARGET_OBJECTS:ToolsOrg>
    $<TARGET_OBJECTS:ToolsPhys> $<TARGET_OBJECTS:ToolsYaml>)
  target_link_libraries(ATOOLSCore PRIVATE ToolsMath ToolsOrg ToolsPhys ToolsYaml)
  sherpa_handle_shared_library(ATOOLSCore ATOOLSCore)
  add_dependencies(ATOOLSCore Math_Git_Info Org_Git_Info Phys_Git_Info YAML_Git_Info)
endif()
