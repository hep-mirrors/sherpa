##### Automatically or semiautomaticaly generated, do not edit.#########
#Analyses/Makefile.am
 
set(GITTAG  AddOns/Analysis/Analyses)

set(Analyses_SOURCES Analysis_Base.C
                     WPolarization_Analysis.C
                     ATLAS_CONF_2011_038.C
                     xQ2_Analysis.C)

set(Analyses_HEADERS Analysis_Base.H)
sherpa_create_git_info(Analyses ../Analyses/)
list(TRANSFORM Analyses_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Analyses/")
list(APPEND Analyses_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaAnalyses SHARED ${Analyses_SOURCES} )
sherpa_handle_shared_library(SherpaAnalyses SherpaAnalyses)
