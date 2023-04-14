set(Single_Events_SOURCES Event_Handler.C
                          Event_Phase_Handler.C
                          Signal_Processes.C
                          Multiple_Interactions.C
                          Minimum_Bias.C
                          Jet_Evolution.C
                          Signal_Process_FS_QED_Correction.C
                          Hard_Decays.C
                          Beam_Remnants.C
                          Hadronization.C
                          Hadron_Decays.C
                          Userhook_Phase.C
                          Analysis_Phase.C
                          Output_Phase.C
                          EvtReadin_Phase.C
                          Soft_Photons.C
                          Decay_Handler_Base.C)

set(Single_Events_HEADERS Event_Handler.H
                          Event_Phase_Handler.H
                          Signal_Processes.H
                          Multiple_Interactions.H
                          Minimum_Bias.H
                          Jet_Evolution.H
                          Signal_Process_FS_QED_Correction.H
                          Hard_Decays.H
                          Beam_Remnants.H
                          Hadronization.H
                          Hadron_Decays.H
                          Userhook_Phase.H
                          Analysis_Phase.H
                          Output_Phase.H
                          EvtReadin_Phase.H
                          Soft_Photons.H
                          Decay_Handler_Base.H)

sherpa_create_git_info(Single_Events)
 
add_library(SherpaSingleEvents SHARED ${Single_Events_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaSingleEvents SherpaSingleEvents)
