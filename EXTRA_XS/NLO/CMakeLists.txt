set(ExtraXSNLO_SOURCES Logarithms.C
                Dummy_QCD_Virtual.C
                DY_QCD_Virtual.C
                DY_Jet_QCD_Virtual.C
                Higgs_QCD_Virtual.C
                Higgs_Jet_QCD_Virtual.C
                ee3jet_CSS_approx.C
                ppll1jet_CSS_approx.C
                epe2jet_CSS_approx.C)

set(ExtraXSNLO_HEADERS Logarithms.H)

sherpa_create_git_info(ExtraXSNLO ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(ExtraXSNLO SHARED ${ExtraXSNLO_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ExtraXSNLO ExtraXSNLO)
add_dependencies(ExtraXSNLO ExtraXSNLO_Git_Info)
