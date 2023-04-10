add_subdirectory(Main)
add_subdirectory(Current_Library)
add_subdirectory(Tools)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/Transitiondata.yaml DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} COMPONENT doc  )
