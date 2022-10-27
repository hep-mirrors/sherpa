
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Scales/Makefile.am
 
 #pkglib_LTLIBRARIES = libPhasicScales.la)
 set(GITTAG  PHASIC++/Scales)
 #include ../../git.make)
set(Scales_SOURCES  Cluster_Definitions.C Scale_Setter_Base.C Core_Scale_Setter.C Default_Core_Scale.C Higgs_Core_Scale.C KFactor_Setter_Base.C Color_Setter.C Tag_Setter.C MEPS_Scale_Setter.C VBF_Scale_Setter.C
 Variable_Scale_Setter.C Variable_Core_Scale.C Fastjet_Scale_Setter.C Democratic_Scale_Setter.C QCD_Scale_Setter.C QCD_Core_Scale.C TTBar_Core_Scale.C TZj_Core_Scale.C SingleTop_Core_Scale.C QQSinglet_Core_Scale.C
 HTPrime_Core_Scale.C MINLO_Scale_Setter.C MINLO_KFactor_Setter.C EWVirtKFactor_Setter.C Variable_KFactor_Setter.C No_KFactor_Setter.C)
set(Scales_HEADERS Cluster_Definitions.H Scale_Setter_Base.H Core_Scale_Setter.H KFactor_Setter_Base.H MINLO_Scale_Setter.H Color_Setter.H Tag_Setter.H)
creategitinfo(Scales ../Scales/)
  list(TRANSFORM Scales_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Scales/")
 
set(Scales_esources )
if (SHERPA_BUILD_STATIC)
add_library(PhasicScales_static STATIC ${PhasicScales_esources}
                             ${Scales_SOURCES} 
)
target_include_directories(PhasicScales_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PhasicScales_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhasicScales_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhasicScales_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME PhasicScales)
endif()
if (SHERPA_BUILD_SHARED)
add_library(PhasicScales SHARED ${PhasicScales_esources}
                             ${Scales_SOURCES} 
)
target_include_directories(PhasicScales PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PhasicScales PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhasicScales DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhasicScales        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME PhasicScales SOVERSION 0.0.0)
endif()
