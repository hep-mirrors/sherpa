
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libMEToolsMain.la)
 set(GITTAG  METOOLS/Main)
 #include ../../git.make)
set(Main_SOURCES  Polarization_Tools.C XYZFuncs.C Polarization_Index.C Spin_Structure.C Partial_Amplitude_Base.C Three_Particle_Amplitudes.C Four_Particle_Amplitudes.C)
set(Main_HEADERS Polarization_Tools.H XYZFuncs.H Polarization_Index.H Spin_Structure.H Partial_Amplitude_Base.H Three_Particle_Amplitudes.H Four_Particle_Amplitudes.H)
creategitinfo(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(METoolsMain SHARED 
                             ${Main_SOURCES} 
)
target_include_directories(METoolsMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(METoolsMain)
install(TARGETS METoolsMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(METoolsMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME METoolsMain SOVERSION 0.0.0)
set_target_properties(METoolsMain        PROPERTIES DEFINE_SYMBOL "")
