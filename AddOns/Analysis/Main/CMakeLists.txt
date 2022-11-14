##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
 set(GITTAG  AddOns/Analysis/Main)
set(Main_SOURCES  Analysis_Object.C Primitive_Analysis.C Analysis_Handler.C)
set(ANAMAIN_ADDS  Analysis_Object.H Primitive_Analysis.H Analysis_Handler.H)
 if (SHERPA_ENABLE_ANALYSIS)
 #pkglib_LTLIBRARIES = libSherpaAnalysis.la)
set(Main_SOURCES ${ANAMAIN_SOURCES})
set(Main_HEADERS ${ANAMAIN_ADDS})
#libSherpaAnalysis_la_LIBADD = -L../Tools -L../Triggers -L../Observables -L../Analyses -lSherpaAnalysisTools -lSherpaAnalysisTrigger -lSherpaObservables -lSherpaAnalyses)
 #ANAMAIN_#EXTRA_DIST = ${ANAMAIN_ADDS})
 else()
 #ANAMAIN_#EXTRA_DIST = ${ANAMAIN_SOURCES} ${ANAMAIN_ADDS})
 endif()
 #EXTRA_DIST = $(#ANAMAIN_#EXTRA_DIST)
sherpa_create_git_info(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaAnalysis SHARED ${Main_SOURCES} )
sherpa_handle_shared_library(SherpaAnalysis SherpaAnalysis)
