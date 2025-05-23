set(Scales_SOURCES Cluster_Definitions.C
                   Scale_Setter_Base.C
                   Core_Scale_Setter.C
                   Default_Core_Scale.C
                   Higgs_Core_Scale.C
                   KFactor_Setter_Base.C
                   Color_Setter.C
                   Tag_Setter.C
                   MEPS_Scale_Setter.C
                   VBF_Scale_Setter.C
                   Variable_Scale_Setter.C
                   Variable_Core_Scale.C
                   Fastjet_Scale_Setter.C
                   Democratic_Scale_Setter.C
                   QCD_Core_Scale.C
                   Photons_Core_Scale.C
                   TTBar_Core_Scale.C
                   TZj_Core_Scale.C
                   SingleTop_Core_Scale.C
                   QQSinglet_Core_Scale.C
                   HTPrime_Core_Scale.C
		           Photons_Core_Scale.C
                   MINLO_Scale_Setter.C
                   MINLO_KFactor_Setter.C
                   EWVirtKFactor_Setter.C
                   Variable_KFactor_Setter.C
                   No_KFactor_Setter.C)

set(Scales_HEADERS Cluster_Definitions.H
                   Scale_Setter_Base.H
                   Core_Scale_Setter.H
                   KFactor_Setter_Base.H
                   MINLO_Scale_Setter.H
                   Color_Setter.H
                   Tag_Setter.H)

sherpa_create_git_info(Scales ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(PhasicScales SHARED ${Scales_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(PhasicScales PhasicScales)
add_dependencies(PhasicScales Scales_Git_Info)
