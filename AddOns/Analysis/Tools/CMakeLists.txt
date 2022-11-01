
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Tools/Makefile.am
 
 set(GITTAG  AddOns/Analysis/Tools)
 #include ../../../git.make)
set(Tools_SOURCES  Lorentz.C Breit_Frame.C MET.C Blob_Selector.C List_Creator.C List_Merger.C Particle_Ordering.C Particle_Qualifier.C Particle_Selector.C)
set(ANATOOL_ADDS  Particle_Qualifier.H)
 if (SHERPA_ENABLE_ANALYSIS)
 #pkglib_LTLIBRARIES = libSherpaAnalysisTools.la)
set(Tools_SOURCES ${ANATOOLS_SOURCES})
set(Tools_HEADERS ${ANATOOLS_ADDS})
 #ANATOOLS_#EXTRA_DIST = ${ANATOOLS_ADDS})
 else()
 #ANATOOLS_#EXTRA_DIST = ${ANATOOLS_SOURCES} ${ANATOOLS_ADDS})
 endif()
 #EXTRA_DIST = $(#ANATOOLS_#EXTRA_DIST)
creategitinfo(Tools ../Tools/)
list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
list(APPEND Tools_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(Tools_esources )
if (SHERPA_BUILD_SHARED)
add_library(SherpaAnalysisTools SHARED ${SherpaAnalysisTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(SherpaAnalysisTools PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaAnalysisTools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaAnalysisTools DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaAnalysisTools        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaAnalysisTools SOVERSION 0.0.0)
set_target_properties(SherpaAnalysisTools        PROPERTIES DEFINE_SYMBOL "")
endif()
