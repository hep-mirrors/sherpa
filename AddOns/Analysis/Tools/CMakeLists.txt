##### Automatically or semiautomaticaly generated, do not edit.#########
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
sherpa_create_git_info(Tools ../Tools/)
list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
list(APPEND Tools_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaAnalysisTools SHARED ${Tools_SOURCES} )
sherpa_handle_shared_library(SherpaAnalysisTools SherpaAnalysisTools)
