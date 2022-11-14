##### Automatically or semiautomaticaly generated, do not edit.#########
#Calculators/Makefile.am
 
set(GITTAG  MCATNLO/Calculators)
set(Calculators_SOURCES  CF_QCD.C LF_FFV.C LF_VVV.C)
sherpa_create_git_info(Calculators ../Calculators/)
list(TRANSFORM Calculators_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Calculators/")
list(APPEND Calculators_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(MCatNLOCalculators SHARED ${Calculators_SOURCES} )
sherpa_handle_shared_library(MCatNLOCalculators MCatNLOCalculators)
