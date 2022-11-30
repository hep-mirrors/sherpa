##### Automatically or semiautomaticaly generated, do not edit.#########
#Pythia/Makefile.am
 
 if (SHERPA_ENABLE_PYTHIA8)
set(Pythia_SOURCES Pythia_Jet_Criterion.C
                   Pythia8_Hadronisation.C)

 else()
set(Pythia_SOURCES Pythia_Jet_Criterion.C)

 endif()
set(GITTAG  AddOns/Pythia)

sherpa_create_git_info(Pythia)
 
add_library(SherpaPythia SHARED ${Pythia_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaPythia SherpaPythia)

if (SHERPA_ENABLE_PYTHIA8)   
  target_link_libraries(SherpaPythia PRIVATE ${PYTHIA8_LIBRARIES})
  target_include_directories(SherpaPythia PRIVATE ${PYTHIA8_INCLUDE_DIRS})
endif()
