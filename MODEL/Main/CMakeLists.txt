##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
 set(GITTAG  MODEL/Main)
set(Main_SOURCES  Color_Function.C Single_Vertex.C Coupling_Data.C Model_Base.C Running_AlphaS.C Strong_Coupling.C Running_AlphaQED.C Running_Fermion_Mass.C)
set(Main_HEADERS Color_Function.H Single_Vertex.H Coupling_Data.H Model_Base.H Running_AlphaS.H Strong_Coupling.H Running_AlphaQED.H Running_Fermion_Mass.H)
sherpa_create_git_info(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ModelMain SHARED ${Main_SOURCES} )
sherpa_handle_shared_library(ModelMain ModelMain)
