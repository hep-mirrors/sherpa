add_subdirectory(Ladders)
add_subdirectory(Beam_Remnants)
add_subdirectory(Cross_Sections)
add_subdirectory(Eikonals)
add_subdirectory(Event_Generation)
add_subdirectory(Main)
add_subdirectory(Tools)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(Shrimps SHARED
    $<TARGET_OBJECTS:ShrimpsLadders> $<TARGET_OBJECTS:ShrimpsBeamRemnants>
    $<TARGET_OBJECTS:ShrimpsXsecs> $<TARGET_OBJECTS:ShrimpsEikonals>
    $<TARGET_OBJECTS:ShrimpsEvents> $<TARGET_OBJECTS:ShrimpsMain>
    $<TARGET_OBJECTS:ShrimpsTools>)
  target_link_libraries(Shrimps PRIVATE ShrimpsLadders ShrimpsBeamRemnants ShrimpsXsecs ShrimpsEikonals ShrimpsEvents ShrimpsMain ShrimpsTools)
  sherpa_handle_shared_library(Shrimps Shrimps)
  add_dependencies(Shrimps Ladders_Git_Info Beam_Remnants_Git_Info Cross_Sections_Git_Info Eikonals_Git_Info Event_Generation_Git_Info ShrimpsMain_Git_Info ShrimpsTools_Git_Info)
endif()
