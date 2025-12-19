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

if (SHERPA_ENABLE_LOOPTOOLS) 
    list(APPEND NLO_SOURCES PionPionVirtual.C)
endif()

sherpa_create_git_info(ExtraXSNLO ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(ExtraXSNLO SHARED ${ExtraXSNLO_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
if (SHERPA_ENABLE_LOOPTOOLS) 
    list(APPEND NLO_SOURCES emu_emu_Virtual.C)
    target_include_directories(ExtraXSNLO  PRIVATE ${LOOPTOOLS_DIRS}/include)
    target_link_libraries(ExtraXSNLO  PRIVATE ${LOOPTOOLS_DIRS}/lib/libooptools.a)
endif()
sherpa_handle_shared_library(ExtraXSNLO ExtraXSNLO)
add_dependencies(ExtraXSNLO ExtraXSNLO_Git_Info)