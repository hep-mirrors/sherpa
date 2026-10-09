set(ExtraXSNNLO_SOURCES ee_aa_QED_TwoLoop.C
                        ee_aa_QED_Poles.C
                        ee_aa_QED_VV.C)

set(ExtraXSNNLO_HEADERS ee_aa_QED_TwoLoop.H
                        ee_aa_QED_Poles.H
                        ee_aa_QED_MSeries.inc
                        ee_aa_QED_VV.H)

sherpa_create_git_info(ExtraXSNNLO ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(ExtraXSNNLO ${SHERPA_SUBLIB_TYPE} ${ExtraXSNNLO_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ExtraXSNNLO ExtraXSNNLO)
add_dependencies(ExtraXSNNLO ExtraXSNNLO_Git_Info)
