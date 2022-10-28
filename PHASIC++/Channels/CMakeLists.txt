
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Channels/Makefile.am
 
 #pkglib_LTLIBRARIES = libPhasicChannels.la)
 set(GITTAG  PHASIC++/Channels)
 #include ../../git.make)
set(Channels_SOURCES  CSS_Kinematics.C Channel_Elements.C Channel_Elements_KK.C Single_Channel.C FSR_Channel.C Rambo.C RamboKK.C VHAAG.C VHAAG_ND.C VHAAG_res.C VHAAG_Threshold.C Decay_Dalitz.C Vegas.C
 Simple_Pole_Channels.C Resonance_Channels.C Leading_Log_Channels.C LBS_Compton_Peak_Channels.C Threshold_Channels.C Channel_Interface.C Multi_Channel.C ISR_Channels.C Beam_Channels.C Channel_Generator.C
 FSR_Channels.C BBar_Multi_Channel.C CS_Dipole.C CS_Dipoles.C BBar_Emission_Generator.C Exponential_Channels.C)
set(Channels_HEADERS CSS_Kinematics.H Channel_Elements.H Channel_Elements_KK.H Single_Channel.H FSR_Channel.H Channel_Interface.H Rambo.H RamboKK.H VHAAG.H VHAAG_ND.H VHAAG_res.H Decay_Dalitz.H Vegas.H
 ISR_Channel_Base.H Simple_Pole_Channels.H Resonance_Channels.H Leading_Log_Channels.H LBS_Compton_Peak_Channels.H Threshold_Channels.H Multi_Channel.H ISR_Channels.H Beam_Channels.H Channel_Generator.H
 FSR_Channels.H BBar_Multi_Channel.H CS_Dipole.H CS_Dipoles.H BBar_Emission_Generator.H Exponential_Channels.H)
creategitinfo(Channels ../Channels/)
  list(TRANSFORM Channels_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Channels/")
 
set(Channels_esources )
if (SHERPA_BUILD_STATIC)
add_library(PhasicChannels_static STATIC ${PhasicChannels_esources}
                             ${Channels_SOURCES} 
)
target_include_directories(PhasicChannels_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PhasicChannels_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhasicChannels_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhasicChannels_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME PhasicChannels)
endif()
if (SHERPA_BUILD_SHARED)
add_library(PhasicChannels SHARED ${PhasicChannels_esources}
                             ${Channels_SOURCES} 
)
target_include_directories(PhasicChannels PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PhasicChannels PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhasicChannels DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhasicChannels        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME PhasicChannels SOVERSION 0.0.0)
endif()