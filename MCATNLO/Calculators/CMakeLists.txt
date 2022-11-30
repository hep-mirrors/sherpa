set(Calculators_SOURCES CF_QCD.C
                        LF_FFV.C
                        LF_VVV.C)

sherpa_create_git_info(Calculators)
 
add_library(MCatNLOCalculators SHARED ${Calculators_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(MCatNLOCalculators MCatNLOCalculators)
