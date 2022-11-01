
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Analyses/Makefile.am
 
 set(GITTAG  AddOns/Analysis/Analyses)
 #include ../../../git.make)
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
creategitinfo(Analyses ../Analyses/)
list(TRANSFORM Analyses_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Analyses/")
list(APPEND Analyses_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(Analyses_esources )
if (SHERPA_BUILD_SHARED)
add_library(SherpaAnalyses SHARED ${SherpaAnalyses_esources}
                             ${Analyses_SOURCES} 
)
target_include_directories(SherpaAnalyses PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaAnalyses PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaAnalyses DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaAnalyses        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaAnalyses SOVERSION 0.0.0)
set_target_properties(SherpaAnalyses        PROPERTIES DEFINE_SYMBOL "")
endif()
