##### Automatically or semiautomaticaly generated, do not edit.#########
#Analyses/Makefile.am
 
 set(GITTAG  AddOns/Analysis/Analyses)
set(Analyses_SOURCES  Analysis_Base.C WPolarization_Analysis.C ATLAS_CONF_2011_038.C xQ2_Analysis.C)
set(ANAANA_ADDS  Analysis_Base.H)
 if (SHERPA_ENABLE_ANALYSIS)
 #pkglib_LTLIBRARIES = libSherpaAnalyses.la)
set(Analyses_SOURCES ${ANAANA_SOURCES})
set(Analyses_HEADERS ${ANAANA_ADDS})
 #ANAANA_#EXTRA_DIST = ${ANAANA_ADDS})
 else()
 #ANAANA_#EXTRA_DIST = ${ANAANA_SOURCES} ${ANAANA_ADDS})
 endif()
 #EXTRA_DIST = $(#ANAANA_#EXTRA_DIST)
sherpa_create_git_info(Analyses ../Analyses/)
list(TRANSFORM Analyses_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Analyses/")
list(APPEND Analyses_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaAnalyses SHARED ${Analyses_SOURCES} )
sherpa_handle_shared_library(SherpaAnalyses SherpaAnalyses)
