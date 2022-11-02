
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 set(GITTAG  AddOns/Analysis/Main)
 #include ../../../git.make)
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
creategitinfo(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaAnalysis SHARED 
                             ${Main_SOURCES} 
)
target_include_directories(SherpaAnalysis PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaAnalysis DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaAnalysis        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaAnalysis SOVERSION 0.0.0)
set_target_properties(SherpaAnalysis        PROPERTIES DEFINE_SYMBOL "")
