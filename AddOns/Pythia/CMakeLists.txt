##### Automatically or semiautomaticaly generated, do not edit.#########
#Pythia/Makefile.am
 
 if (SHERPA_ENABLE_PYTHIA8)
set(Pythia_SOURCES Pythia_Jet_Criterion.C
                   Pythia8_Hadronisation.C)

 else()
set(Pythia_SOURCES Pythia_Jet_Criterion.C)

 endif()
set(GITTAG  AddOns/Pythia)
sherpa_create_git_info(Pythia ../Pythia/)
list(TRANSFORM Pythia_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Pythia/")
list(APPEND Pythia_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaPythia SHARED ${Pythia_SOURCES} )
sherpa_handle_shared_library(SherpaPythia SherpaPythia)

if (SHERPA_ENABLE_PYTHIA8)   
  target_link_libraries(SherpaPythia PRIVATE ${PYTHIA8_LIBRARIES})
  target_include_directories(SherpaPythia PRIVATE ${PYTHIA8_INCLUDE_DIRS})
endif()
