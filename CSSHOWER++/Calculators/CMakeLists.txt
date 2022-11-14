##### Automatically or semiautomaticaly generated, do not edit.#########
#Calculators/Makefile.am
 
set(GITTAG  CSSHOWER++/Calculators)
set(Calculators_SOURCES  CF_QCD.C CF_QED.C CF_EW_FFV.C CF_Triangle.C CF_HV.C LF_SSV.C LF_FFV.C LF_VVV.C LF_Triangle.C)
sherpa_create_git_info(Calculators ../Calculators/)
list(TRANSFORM Calculators_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Calculators/")
list(APPEND Calculators_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(CSCalculators SHARED ${Calculators_SOURCES} )
sherpa_handle_shared_library(CSCalculators CSCalculators)
