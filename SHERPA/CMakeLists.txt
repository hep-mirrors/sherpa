add_subdirectory(Main)
add_subdirectory(Initialization)
add_subdirectory(PerturbativePhysics)
add_subdirectory(Single_Events)
add_subdirectory(SoftPhysics)
add_subdirectory(Tools)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(SherpaCore SHARED
    $<TARGET_OBJECTS:SherpaMain> $<TARGET_OBJECTS:SherpaInitialization>
    $<TARGET_OBJECTS:SherpaPerturbativePhysics> $<TARGET_OBJECTS:SherpaSingleEvents>
    $<TARGET_OBJECTS:SherpaSoftPhysics> $<TARGET_OBJECTS:SherpaTools>)
  target_link_libraries(SherpaCore PRIVATE SherpaMain SherpaInitialization SherpaPerturbativePhysics SherpaSingleEvents SherpaSoftPhysics SherpaTools)
  target_link_libraries(SherpaCore PRIVATE YFS)
  sherpa_handle_shared_library(SherpaCore SherpaCore)
  add_dependencies(SherpaCore SherpaMain_Git_Info Initialization_Git_Info PerturbativePhysics_Git_Info Single_Events_Git_Info SoftPhysics_Git_Info SherpaTools_Git_Info)
endif()
