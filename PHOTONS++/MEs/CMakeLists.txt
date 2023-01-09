set(MEs_SOURCES PHOTONS_ME_Base.C
                Scalar_To_Scalar_Scalar.C
                Scalar_To_Fermion_Fermion.C
                Scalar_To_Scalar_Lepton_Neutrino.C
                Scalar_To_Vector_Lepton_Neutrino.C
                Tau_To_Lepton_Neutrinos.C
                W_To_Lepton_Neutrino.C
                Z_To_Fermion_Fermion.C
                Vector_To_Fermion_Fermion.C
                Vector_To_Scalar_Scalar.C
                Collinear_Approximation_FI.C 
                Collinear_Approximation_FF.C 
                EW_One_Loop_Functions_Base.C
                Higgs_Decay_EW_One_Loop_Functions.C 
                Higgs_To_Fermion_Fermion.C 
                W_Decay_EW_One_Loop_Functions.C 
                Z_Decay_EW_One_Loop_Functions.C )



set(NNLO_SOURCES Higgs_Decay_Two_Loop_Functions.C
                 Higgs_Decay_RV_Box_1.C
                 Higgs_Decay_RV_Box_2.C
                 Higgs_Decay_RV_Bubbles.C
                 Higgs_Decay_RV_CT.C
                 Higgs_Decay_RV_Diagrams.C
                 Higgs_Decay_RV_Vertices.C
                 Z_Decay_RV_Box_1.C 
                 Z_Decay_RV_Box_2.C 
                 Z_Decay_RV_Bubbles.C 
                 Z_Decay_RV_CT.C 
                 Z_Decay_RV_Diagrams.C 
                 Z_Decay_RV_Vertices.C 
                 Z_Decay_Two_Loop_Functions.C)


set(MEs_HEADERS PHOTONS_ME_Base.H
                Scalar_To_Scalar_Scalar.H
                Scalar_To_Fermion_Fermion.H
                Scalar_To_Scalar_Lepton_Neutrino.H
                Scalar_To_Vector_Lepton_Neutrino.H
                Tau_To_Lepton_Neutrinos.H
                W_To_Lepton_Neutrino.H
                Z_To_Fermion_Fermion.H
                Vector_To_Fermion_Fermion.H
                Vector_To_Scalar_Scalar.H
                EW_One_Loop_Functions_Base.H
                Higgs_Decay_Two_Loop_Functions.H
                Higgs_Decay_RV_Diagrams.H
                W_Decay_EW_One_Loop_Functions.H
                Z_Decay_EW_One_Loop_Functions.H   
                Z_Decay_Two_Loop_Functions.H
                Z_Decay_RV_Diagrams.H)

sherpa_create_git_info(MEs)
 
if(SHERPA_ENABLE_YFSNNLO)
    add_library(PhotonsMEs SHARED ${MEs_SOURCES} ${NNLO_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C  )
else()
    add_library(PhotonsMEs SHARED ${MEs_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
endif()
sherpa_handle_shared_library(PhotonsMEs PhotonsMEs)
