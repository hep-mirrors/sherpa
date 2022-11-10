
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libBeamMain.la)
 set(GITTAG  BEAM/Main)
 #include ../../git.make)
set(Main_SOURCES  Beam_Base.C Beam_Parameters.C Kinematics_Base.C Collider_Kinematics.C RelicDensity_Kinematics.C DM_Annihilation_Kinematics.C Weight_Base.C Collider_Weight.C RelicDensity_Weight.C DM_Annihilation_Weight.C
 Beam_Spectra_Handler.C)
set(Main_HEADERS Beam_Base.H Kinematics_Base.H Collider_Kinematics.H RelicDensity_Kinematics.H DM_Annihilation_Kinematics.H Weight_Base.H Collider_Weight.H RelicDensity_Weight.H DM_Annihilation_Weight.H
 Beam_Parameters.H Beam_Spectra_Handler.H)
creategitinfo(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(BeamMain SHARED 
                             ${Main_SOURCES} 
)
target_include_directories(BeamMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(BeamMain)
install(TARGETS BeamMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(BeamMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME BeamMain SOVERSION 0.0.0)
set_target_properties(BeamMain        PROPERTIES DEFINE_SYMBOL "")
