set(Currents_SOURCES C_Scalar.C
                     C_Spinor.C
                     C_Vector.C
                     C_Pseudo.C
                     S_C.C
                     F_C.C
                     V_C.C
                     P_C.C)

set(Currents_HEADERS C_Scalar.H
                     C_Spinor.H
                     C_Vector.H
                     C_Pseudo.H)

sherpa_create_git_info(Currents ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(METoolsCurrents SHARED ${Currents_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(METoolsCurrents METoolsCurrents)
add_dependencies(METoolsCurrents Currents_Git_Info)
if (SHERPA_ENABLE_THREADING)
  target_link_libraries(METoolsCurrents PRIVATE Threads::Threads)
endif()
