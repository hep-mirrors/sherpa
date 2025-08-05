 
 if (SHERPA_ENABLE_PYTHIA8)
set(Pythia_SOURCES Pythia_Jet_Criterion.C
                   Pythia8_Hadronisation.C)

 else()
set(Pythia_SOURCES Pythia_Jet_Criterion.C)

 endif()

sherpa_create_git_info(Pythia)
 
add_library(SherpaPythia8 SHARED ${Pythia_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaPythia8 SherpaPythia8)

if (SHERPA_ENABLE_PYTHIA8)   
  target_link_libraries(SherpaPythia8 PRIVATE Pythia8::Pythia8)
endif()
