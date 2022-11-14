##### Automatically or semiautomaticaly generated, do not edit.#########
#NLO/Makefile.am
 
set(GITTAG  EXTRA_XS/NLO)
set(NLO_SOURCES  Logarithms.C Dummy_QCD_Virtual.C DY_QCD_Virtual.C DY_Jet_QCD_Virtual.C Higgs_QCD_Virtual.C Higgs_Jet_QCD_Virtual.C ee3jet_CSS_approx.C ppll1jet_CSS_approx.C epe2jet_CSS_approx.C)
set(NLO_HEADERS Logarithms.H)
sherpa_create_git_info(NLO ../NLO/)
list(TRANSFORM NLO_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../NLO/")
list(APPEND NLO_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ExtraXSNLO SHARED ${NLO_SOURCES} )
sherpa_handle_shared_library(ExtraXSNLO ExtraXSNLO)
