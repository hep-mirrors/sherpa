set(Tools_SOURCES Proto_Particle.C
                  Cluster.C
                  Constituents.C
                  Wave_Function.C
                  Multiplet_Constructor.C
                  Single_Transitions.C
                  Double_Transitions.C
                  Hadronisation_Parameters.C
                  Soft_Cluster_Handler.C
                  Singlet_Tools.C
                  Splitter_Base.C
                  Flavour_Selector.C
                  KT_Selector.C
                  Z_Selector.C)

set(Tools_HEADERS Proto_Particle.H
                  Cluster.H
                  Constituents.H
                  Wave_Function.H
                  Hadron_Multiplet.H
                  Multiplet_Constructor.H
                  Single_Transitions.H
                  Double_Transitions.H
                  Hadronisation_Parameters.H
                  Soft_Cluster_Handler.H
                  Singlet_Tools.H
                  Splitter_Base.H
                  Flavour_Selector.H
                  KT_Selector.H
                  Z_Selector.H)

sherpa_create_git_info(Tools)
 
add_library(AhadicTools SHARED ${Tools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(AhadicTools AhadicTools)
