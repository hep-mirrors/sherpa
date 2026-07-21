add_subdirectory(MEs)
add_subdirectory(Main)
add_subdirectory(PhaseSpace)
add_subdirectory(Tools)
add_subdirectory(PhotonSplitter)
if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(Photons SHARED
    $<TARGET_OBJECTS:PhotonsMEs> $<TARGET_OBJECTS:PhotonsMain>
    $<TARGET_OBJECTS:PhotonsPhaseSpace> $<TARGET_OBJECTS:PhotonsTools>
    $<TARGET_OBJECTS:PhotonsPhotonSplitter>)
  target_link_libraries(Photons PRIVATE PhotonsMEs PhotonsMain PhotonsPhaseSpace PhotonsTools PhotonsPhotonSplitter)
  sherpa_handle_shared_library(Photons Photons)
  add_dependencies(Photons MEs_Git_Info PhotonsMain_Git_Info PhaseSpace_Git_Info PhotonsTools_Git_Info PhotonSplitter_Git_Info)
endif()
