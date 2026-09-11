add_subdirectory(Current_Library)
add_subdirectory(ME_Library)
add_subdirectory(Main)
add_subdirectory(PS_Library)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/Decaydata.yaml DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} COMPONENT doc  )

if (SHERPA_ENABLE_MERGED_LIBS)
  add_library(Hadrons SHARED
    $<TARGET_OBJECTS:HadronsCurrents> $<TARGET_OBJECTS:HadronsMEs>
    $<TARGET_OBJECTS:HadronsMain> $<TARGET_OBJECTS:HadronsPSs>)
  target_link_libraries(Hadrons PRIVATE HadronsCurrents HadronsMEs HadronsMain HadronsPSs)
  sherpa_handle_shared_library(Hadrons Hadrons)
  add_dependencies(Hadrons Current_Library_Git_Info ME_Library_Git_Info HadronsMain_Git_Info PS_Library_Git_Info)
endif()
