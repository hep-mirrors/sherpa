set(MEs_SOURCES PHOTONS_ME_Base.C
                Scalar_To_Scalar_Scalar.C
                Scalar_To_Fermion_Fermion.C
                Scalar_To_Scalar_Lepton_Neutrino.C
                Scalar_To_Vector_Lepton_Neutrino.C
                Tau_To_Lepton_Neutrinos.C
                W_To_Lepton_Neutrino.C
                Z_To_Fermion_Fermion.C
                Vector_To_Fermion_Fermion.C
                Vector_To_Scalar_Scalar.C)

set(MEs_HEADERS PHOTONS_ME_Base.H
                Scalar_To_Scalar_Scalar.H
                Scalar_To_Fermion_Fermion.H
                Scalar_To_Scalar_Lepton_Neutrino.H
                Scalar_To_Vector_Lepton_Neutrino.H
                Tau_To_Lepton_Neutrinos.H
                W_To_Lepton_Neutrino.H
                Z_To_Fermion_Fermion.H
                Vector_To_Fermion_Fermion.H
                Vector_To_Scalar_Scalar.H)

sherpa_create_git_info(MEs ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(PhotonsMEs SHARED ${MEs_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(PhotonsMEs PhotonsMEs)
add_dependencies(PhotonsMEs MEs_Git_Info)
