
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Single_Events/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaSingleEvents.la)
 set(GITTAG  SHERPA/Single_Events)
 #include ../../git.make)
set(Single_Events_SOURCES  Event_Handler.C Event_Phase_Handler.C Signal_Processes.C Multiple_Interactions.C Minimum_Bias.C Jet_Evolution.C Signal_Process_FS_QED_Correction.C Hard_Decays.C Beam_Remnants.C
 Hadronization.C Hadron_Decays.C Userhook_Phase.C Analysis_Phase.C Output_Phase.C EvtReadin_Phase.C Decay_Handler_Base.C)
set(Single_Events_HEADERS Event_Handler.H Event_Phase_Handler.H Signal_Processes.H Multiple_Interactions.H Minimum_Bias.H Jet_Evolution.H Signal_Process_FS_QED_Correction.H Hard_Decays.H Beam_Remnants.H
 Hadronization.H Hadron_Decays.H Userhook_Phase.H Analysis_Phase.H Output_Phase.H EvtReadin_Phase.H Decay_Handler_Base.H)
creategitinfo(Single_Events ../Single_Events/)
list(TRANSFORM Single_Events_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Single_Events/")
list(APPEND Single_Events_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaSingleEvents SHARED 
                             ${Single_Events_SOURCES} 
)
target_include_directories(SherpaSingleEvents PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaSingleEvents DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaSingleEvents        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaSingleEvents SOVERSION 0.0.0)
set_target_properties(SherpaSingleEvents        PROPERTIES DEFINE_SYMBOL "")
