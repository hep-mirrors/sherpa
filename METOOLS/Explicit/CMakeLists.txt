set(Explicit_SOURCES C_Object.C
                     Current.C
                     Vertex_Key.C
                     Color_Calculator.C
                     Lorentz_Calculator.C
                     NLO_Counter_Terms.C
                     Dipole_Kinematics.C
                     Dipole_Color.C
                     Dipole_Terms.C
                     Form_Factor.C
                     Vertex.C)

set(Explicit_HEADERS C_Object.H
                     Current.H
                     Vertex_Key.H
                     Color_Calculator.H
                     Lorentz_Calculator.H
                     NLO_Counter_Terms.H
                     Dipole_Info.H
                     Dipole_Kinematics.H
                     Dipole_Color.H
                     Dipole_Terms.H
                     Form_Factor.H
                     Vertex.H)

sherpa_create_git_info(Explicit ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(METoolsExplicit SHARED ${Explicit_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(METoolsExplicit METoolsExplicit)
add_dependencies(METoolsExplicit Explicit_Git_Info)
if (SHERPA_ENABLE_THREADING)
  target_link_libraries(METoolsExplicit PRIVATE Threads::Threads)
endif()
