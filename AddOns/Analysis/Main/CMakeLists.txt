set(Main_SOURCES Analysis_Object.C
                 Primitive_Analysis.C
                 Analysis_Handler.C)

set(Main_HEADERS Analysis_Object.H
                 Primitive_Analysis.H
                 Analysis_Handler.H)

sherpa_create_git_info(Main ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(SherpaAnalysis SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
target_link_libraries(SherpaAnalysis PRIVATE SherpaAnalyses SherpaObservables SherpaAnalysisTools SherpaAnalysisTrigger)
sherpa_handle_shared_library(SherpaAnalysis SherpaAnalysis)
