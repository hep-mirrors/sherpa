add_subdirectory(Current_Library)
add_subdirectory(ME_Library)
add_subdirectory(Main)
add_subdirectory(PS_Library)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/Decaydata.yaml DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} COMPONENT doc  )
