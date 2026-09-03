set(ExtraXSNLO_SOURCES Logarithms.C
                Dummy_QCD_Virtual.C
                DY_QCD_Virtual.C
                DY_Jet_QCD_Virtual.C
                Higgs_QCD_Virtual.C
                Higgs_Jet_QCD_Virtual.C
                ee3jet_CSS_approx.C
                ppll1jet_CSS_approx.C
                epe2jet_CSS_approx.C
                eenunu_TwoLoop.C)

set(ExtraXSNLO_HEADERS Logarithms.H)

# Compiled unconditionally, NOT under SHERPA_ENABLE_LOOPTOOLS: the file guards
# its own LoopTools half with #ifdef USING__LOOPTOOLS so that a build without
# LoopTools still registers the getter and can tell the user why it cannot run.
# (This list previously said NLO_SOURCES -- not the variable add_library uses --
# so PionPionVirtual.C was never compiled under any setting.)
list(APPEND ExtraXSNLO_SOURCES PionPionVirtual.C)

sherpa_create_git_info(ExtraXSNLO ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(ExtraXSNLO ${SHERPA_SUBLIB_TYPE} ${ExtraXSNLO_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
if (SHERPA_ENABLE_LOOPTOOLS)
    # emu_emu_Virtual.C is NOT added here: this block runs after add_library,
    # so appending to any source list has no effect. Left as it was.
    target_link_libraries(ExtraXSNLO  PRIVATE LoopTools::LoopTools)
endif()
sherpa_handle_shared_library(ExtraXSNLO ExtraXSNLO)
add_dependencies(ExtraXSNLO ExtraXSNLO_Git_Info)