set(Analyses_SOURCES Analysis_Base.C
                     WPolarization_Analysis.C
                     ATLAS_CONF_2011_038.C
                     xQ2_Analysis.C)

set(Analyses_HEADERS Analysis_Base.H)

sherpa_create_git_info(Analyses)
 
add_library(SherpaAnalyses SHARED ${Analyses_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaAnalyses SherpaAnalyses)
