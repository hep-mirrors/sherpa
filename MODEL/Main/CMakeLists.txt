set(ModelMain_SOURCES Color_Function.C
                 Single_Vertex.C
                 Coupling_Data.C
                 Model_Base.C
                 Running_AlphaS.C
                 Strong_Coupling.C
                 Running_AlphaQED.C
                 Running_Fermion_Mass.C)

set(ModelMain_HEADERS Color_Function.H
                 Single_Vertex.H
                 Coupling_Data.H
                 Model_Base.H
                 Running_AlphaS.H
                 Strong_Coupling.H
                 Running_AlphaQED.H
                 Running_Fermion_Mass.H)

sherpa_create_git_info(ModelMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
if (SHERPA_ENABLE_ALPHAQEDHAD) 
    list(APPEND ModelMain_SOURCES hadr5x19.f)
endif()

add_library(ModelMain SHARED ${ModelMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
add_dependencies(ModelMain ModelMain_Git_Info)

if (SHERPA_ENABLE_ALPHAQEDHAD) 
    sherpa_handle_shared_library(ModelMain ModelMain hadr5x19)
else()
    sherpa_handle_shared_library(ModelMain ModelMain)
endif()
