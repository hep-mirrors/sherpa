set(ME_Library_SOURCES HD_ME_Base.C
                       Generic.C
                       Current_ME.C
                       Eta_Decay_MEs.C
                       Four_Body_MEs.C
                       Three_Body_MEs.C
                       Two_Body_MEs.C
                       B_K_Semileptonic.C
                       B_KStar_Semileptonic.C)

set(ME_Library_HEADERS HD_ME_Base.H
                       Generic.H
                       Current_ME.H
                       PseudoScalar_Decay_MEs.H
                       Eta_Decay_MEs.H
                       Vector_Decay_MEs.H
                       Phi_Decay_MEs.H
                       Tensor_Decay_MEs.H
                       Baryon_SimpleDecay_MEs.H
                       HeavyBaryon_SimpleTwoBody_Decay_MEs.H
                       Four_Body_MEs.H
                       Three_Body_MEs.H
                       Two_Body_MEs.H
                       B_K_Semileptonic.H
                       B_KStar_Semileptonic.H)

sherpa_create_git_info(ME_Library)
 
add_library(HadronsMEs SHARED ${ME_Library_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(HadronsMEs HadronsMEs)
