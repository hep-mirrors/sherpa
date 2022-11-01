
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Triggers/Makefile.am
 
 set(GITTAG  AddOns/Analysis/Triggers)
 #include ../../../git.make)
set(Triggers_SOURCES  Trigger_Base.C One_Variable_Selector.C One_Particle_Selector.C Two_Particle_Selector.C All_Particle_Selector.C Jet_Particle_Selector.C Two_Particle_X_Selector.C Custom_Selector.C
 ranlux.cpp reference.cpp geom_2d.cpp momentum.cpp hash.cpp quadtree.cpp vicinity.cpp protocones.cpp split_merge.cpp siscone.cpp siscone_error.cpp area.cpp MySISCone.C ConeMCFM.C MCFMCone.C HT_Selector.C
 Q2_Selector.C Durham_Algorithm.C Cambridge_Algorithm.C Kt_Algorithm.C DIS_Algorithm.C Final_Selector.C Midpoint_Cone.C Calorimeter_Cone.C FrixioneIsolation_Selector.C Isolation_Selector.C Fastjet_Interface.C
 Primitive_Detector.C Primitive_Calorimeter.C)
set(ANATRIGGER_ADDS  Trigger_Base.H Durham_Algorithm.H Cambridge_Algorithm.H area.h circulator.h defines.h geom_2d.h hash.h momentum.h protocones.h quadtree.h ranlux.h reference.h siscone_error.h siscone.h
 split_merge.h vicinity.h MySISCone.H ConeMCFM.H MCFMCone.H Kt_Algorithm.H DIS_Algorithm.H Final_Selector.H Midpoint_Cone.H Calorimeter_Cone.H Primitive_Detector.H Primitive_Calorimeter.H)
 if (SHERPA_ENABLE_ANALYSIS)
 #pkglib_LTLIBRARIES = libSherpaAnalysisTrigger.la)
set(Triggers_SOURCES ${ANATRIGGER_SOURCES})
set(Triggers_HEADERS ${ANATRIGGER_ADDS})
 #ANATRIGGER_#EXTRA_DIST = ${ANATRIGGER_ADDS})
 else()
 #ANATRIGGER_#EXTRA_DIST = ${ANATRIGGER_SOURCES} ${ANATRIGGER_ADDS})
 endif()
 #EXTRA_DIST = $(#ANATRIGGER_#EXTRA_DIST)
creategitinfo(Triggers ../Triggers/)
list(TRANSFORM Triggers_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Triggers/")
list(APPEND Triggers_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(Triggers_esources )
if (SHERPA_BUILD_SHARED)
add_library(SherpaAnalysisTrigger SHARED ${SherpaAnalysisTrigger_esources}
                             ${Triggers_SOURCES} 
)
target_include_directories(SherpaAnalysisTrigger PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaAnalysisTrigger PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaAnalysisTrigger DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaAnalysisTrigger        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaAnalysisTrigger SOVERSION 0.0.0)
set_target_properties(SherpaAnalysisTrigger        PROPERTIES DEFINE_SYMBOL "")
endif()