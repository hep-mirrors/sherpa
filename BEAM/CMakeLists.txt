add_subdirectory(Main)
add_subdirectory(Spectra)

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(Beam SHARED
    $<TARGET_OBJECTS:BeamMain> $<TARGET_OBJECTS:BeamSpectra>)
  target_link_libraries(Beam PRIVATE BeamMain BeamSpectra)
  sherpa_handle_shared_library(Beam Beam)
  add_dependencies(Beam BeamMain_Git_Info Spectra_Git_Info)
endif()
