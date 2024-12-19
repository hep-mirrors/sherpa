set(Channels_SOURCES Transverse_Kinematics.C
                     CSS_Kinematics.C
                     Antenna_Kinematics.C
                     Channel_Elements.C
                     Channel_Elements_KK.C
                     Single_Channel.C
                     FSR_Channel.C
                     Rambo.C
                     RamboKK.C
                     VHAAG.C
                     VHAAG_ND.C
                     VHAAG_res.C
                     VHAAG_Threshold.C
                     Decay_Dalitz.C
                     Vegas.C
                     Simple_Pole_Channels.C
                     Resonance_Channels.C
                     Leading_Log_Channels.C
                     LBS_Compton_Peak_Channels.C
                     Threshold_Channels.C
                     Channel_Interface.C
                     Multi_Channel.C
                     ISR_Channels.C
                     Beam_Channels.C
                     Channel_Generator.C
                     FSR_Channels.C
                     BBar_Multi_Channel.C
                     CS_Dipole.C
                     CS_Dipoles.C
                     BBar_Emission_Generator.C
                     Exponential_Channels.C)

set(Channels_HEADERS Transverse_Kinematics.H
                     CSS_Kinematics.H
                     Antenna_Kinematics.H
                     Channel_Elements.H
                     Channel_Elements_KK.H
                     Single_Channel.H
                     FSR_Channel.H
                     Channel_Interface.H
                     Rambo.H
                     RamboKK.H
                     VHAAG.H
                     VHAAG_ND.H
                     VHAAG_res.H
                     Decay_Dalitz.H
                     Vegas.H
                     ISR_Channel_Base.H
                     Simple_Pole_Channels.H
                     Resonance_Channels.H
                     Leading_Log_Channels.H
                     LBS_Compton_Peak_Channels.H
                     Threshold_Channels.H
                     Multi_Channel.H
                     ISR_Channels.H
                     Beam_Channels.H
                     Channel_Generator.H
                     FSR_Channels.H
                     BBar_Multi_Channel.H
                     CS_Dipole.H
                     CS_Dipoles.H
                     BBar_Emission_Generator.H
                     Exponential_Channels.H)

sherpa_create_git_info(Channels)
 
add_library(PhasicChannels SHARED ${Channels_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(PhasicChannels PhasicChannels)
