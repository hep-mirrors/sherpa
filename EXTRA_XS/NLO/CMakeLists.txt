set(ExtraXSNLO_SOURCES Logarithms.C
                Dummy_QCD_Virtual.C
                DY_QCD_Virtual.C
                DY_Jet_QCD_Virtual.C
                Higgs_QCD_Virtual.C
                Higgs_Jet_QCD_Virtual.C
                PionPionGamma.C
                ee3jet_CSS_approx.C
                ee2nunuA.C
                # ee2FFBarGamma.C
                ppll1jet_CSS_approx.C
                epe2jet_CSS_approx.C)

set(ExtraXSNLO_HEADERS Logarithms.H)
sherpa_create_git_info(ExtraXSNLO ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

if (SHERPA_ENABLE_LOOPTOOLS) 
    list(APPEND ExtraXSNLO_SOURCES PionPionVirtual.C)
endif()

add_library(ExtraXSNLO SHARED ${ExtraXSNLO_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
if (SHERPA_ENABLE_LOOPTOOLS) 
    list(APPEND NLO_SOURCES emu_emu_Virtual.C)
    target_include_directories(ExtraXSNLO  PRIVATE ${LOOPTOOLS_DIRS}/include)
    target_link_libraries(ExtraXSNLO  PRIVATE ${LOOPTOOLS_DIRS}/lib/libooptools.a)
endif()
sherpa_handle_shared_library(ExtraXSNLO ExtraXSNLO)
if (SHERPA_ENABLE_LOOPTOOLS)
    set(LOOPTOOLS_DIRS "/Users/alanprice/software/LoopTools-2.15/arm64-Darwin")
    find_library(LOOPTOOLS_LIB NAMES ooptools PATHS ${LOOPTOOLS_DIRS}/lib)
    if (LOOPTOOLS_LIB)
        target_include_directories(ExtraXSNLO  PRIVATE ${LOOPTOOLS_DIRS}/include)
        target_link_libraries(ExtraXSNLO PRIVATE ${LOOPTOOLS_LIB} gfortran)
        # target_link_libraries(ExtraXSNLO PRIVATE ${LOOPTOOLS_LIB})
    else()
        message(FATAL_ERROR "LoopTools library not found!")
    endif()
endif()

sherpa_create_git_info(ExtraXSNLO ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
add_dependencies(ExtraXSNLO ExtraXSNLO_Git_Info)