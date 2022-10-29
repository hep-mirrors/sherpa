
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
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
 
set(Main_esources )
if (SHERPA_BUILD_STATIC)
add_library(METoolsMain_static STATIC ${METoolsMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(METoolsMain_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(METoolsMain_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS METoolsMain_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(METoolsMain_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME METoolsMain)
set_target_properties(METoolsMain_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(METoolsMain SHARED ${METoolsMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(METoolsMain PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(METoolsMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS METoolsMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(METoolsMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME METoolsMain SOVERSION 0.0.0)
set_target_properties(METoolsMain        PROPERTIES DEFINE_SYMBOL "")
endif()
