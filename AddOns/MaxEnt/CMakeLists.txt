if (SHERPA_ENABLE_JSON)
  set(JSONSOURCES DY_Weight_Hook.C)
endif()

set(MaxEnt_SOURCES Thrust_Weight_Hook.C
		   ${JSONSOURCES})

sherpa_create_git_info(MaxEnt ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaMaxEnt SHARED ${MaxEnt_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaMaxEnt SherpaMaxEnt)
if (SHERPA_ENABLE_JSON)
  target_link_libraries(SherpaMaxEnt PRIVATE nlohmann_json)
endif()
