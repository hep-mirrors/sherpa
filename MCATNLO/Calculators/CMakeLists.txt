set(MCatNLOCalculators_SOURCES CF_QCD.C
                        LF_FFV.C
                        LF_VVV.C)

sherpa_create_git_info(MCatNLOCalculators ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(MCatNLOCalculators SHARED ${MCatNLOCalculators_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(MCatNLOCalculators MCatNLOCalculators)
add_dependencies(MCatNLOCalculators MCatNLOCalculators_Git_Info)
