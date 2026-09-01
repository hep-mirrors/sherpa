set(BeamMain_SOURCES Beam_Base.C
                 Beam_Parameters.C
                 Kinematics_Base.C
                 Collider_Kinematics.C
                 RelicDensity_Kinematics.C
                 DM_Annihilation_Kinematics.C
                 Weight_Base.C
                 Collider_Weight.C
                 RelicDensity_Weight.C
                 DM_Annihilation_Weight.C
                 Beam_Spectra_Handler.C)

set(BeamMain_HEADERS Beam_Base.H
                 Beam_Rejection.H
                 Kinematics_Base.H
                 Collider_Kinematics.H
                 RelicDensity_Kinematics.H
                 DM_Annihilation_Kinematics.H
                 Weight_Base.H
                 Collider_Weight.H
                 RelicDensity_Weight.H
                 DM_Annihilation_Weight.H
                 Beam_Parameters.H
                 Beam_Spectra_Handler.H)

sherpa_create_git_info(BeamMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(BeamMain SHARED ${BeamMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(BeamMain BeamMain)
add_dependencies(BeamMain BeamMain_Git_Info)
