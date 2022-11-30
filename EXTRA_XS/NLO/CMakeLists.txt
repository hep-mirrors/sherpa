set(NLO_SOURCES Logarithms.C
                Dummy_QCD_Virtual.C
                DY_QCD_Virtual.C
                DY_Jet_QCD_Virtual.C
                Higgs_QCD_Virtual.C
                Higgs_Jet_QCD_Virtual.C
                ee3jet_CSS_approx.C
                ppll1jet_CSS_approx.C
                epe2jet_CSS_approx.C)

set(NLO_HEADERS Logarithms.H)

sherpa_create_git_info(NLO)
 
add_library(ExtraXSNLO SHARED ${NLO_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ExtraXSNLO ExtraXSNLO)
