##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
set(GITTAG  BEAM/Main)

set(Main_SOURCES Beam_Base.C
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

set(Main_HEADERS Beam_Base.H
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

sherpa_create_git_info(Main)
 
add_library(BeamMain SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(BeamMain BeamMain)
