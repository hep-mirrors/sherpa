##### Automatically or semiautomaticaly generated, do not edit.#########
#Triggers/Makefile.am
 
set(GITTAG  AddOns/Analysis/Triggers)

set(Triggers_SOURCES Trigger_Base.C
                     One_Variable_Selector.C
                     One_Particle_Selector.C
                     Two_Particle_Selector.C
                     All_Particle_Selector.C
                     Jet_Particle_Selector.C
                     Two_Particle_X_Selector.C
                     Custom_Selector.C
                     ranlux.cpp
                     reference.cpp
                     geom_2d.cpp
                     momentum.cpp
                     hash.cpp
                     quadtree.cpp
                     vicinity.cpp
                     protocones.cpp
                     split_merge.cpp
                     siscone.cpp
                     siscone_error.cpp
                     area.cpp
                     MySISCone.C
                     ConeMCFM.C
                     MCFMCone.C
                     HT_Selector.C
                     Q2_Selector.C
                     Durham_Algorithm.C
                     Cambridge_Algorithm.C
                     Kt_Algorithm.C
                     DIS_Algorithm.C
                     Final_Selector.C
                     Midpoint_Cone.C
                     Calorimeter_Cone.C
                     FrixioneIsolation_Selector.C
                     Isolation_Selector.C
                     Fastjet_Interface.C
                     Primitive_Detector.C
                     Primitive_Calorimeter.C)

set(Trigger_HEADERS Trigger_Base.H
                    Durham_Algorithm.H
                    Cambridge_Algorithm.H
                    area.h
                    circulator.h
                    defines.h
                    geom_2d.h
                    hash.h
                    momentum.h
                    protocones.h
                    quadtree.h
                    ranlux.h
                    reference.h
                    siscone_error.h
                    siscone.h
                    split_merge.h
                    vicinity.h
                    MySISCone.H
                    ConeMCFM.H
                    MCFMCone.H
                    Kt_Algorithm.H
                    DIS_Algorithm.H
                    Final_Selector.H
                    Midpoint_Cone.H
                    Calorimeter_Cone.H
                    Primitive_Detector.H
                    Primitive_Calorimeter.H)

sherpa_create_git_info(Triggers)
 
add_library(SherpaAnalysisTrigger SHARED ${Triggers_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaAnalysisTrigger SherpaAnalysisTrigger)
