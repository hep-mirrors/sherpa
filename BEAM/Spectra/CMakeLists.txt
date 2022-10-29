
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Spectra/Makefile.am
 
 #pkglib_LTLIBRARIES = libBeamSpectra.la)
 set(GITTAG  BEAM/Spectra)
 #include ../../git.make)
set(Spectra_SOURCES  Monochromatic.C Laser_Backscattering.C EPA.C DM_beam.C Spectrum_Reader.C)
set(Spectra_HEADERS Monochromatic.H Laser_Backscattering.H EPA.H DM_beam.H Spectrum_Reader.H)
creategitinfo(Spectra ../Spectra/)
  list(TRANSFORM Spectra_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Spectra/")
 
set(Spectra_esources )
if (SHERPA_BUILD_STATIC)
add_library(BeamSpectra_static STATIC ${BeamSpectra_esources}
                             ${Spectra_SOURCES} 
)
target_include_directories(BeamSpectra_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(BeamSpectra_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS BeamSpectra_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(BeamSpectra_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME BeamSpectra)
set_target_properties(BeamSpectra_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(BeamSpectra SHARED ${BeamSpectra_esources}
                             ${Spectra_SOURCES} 
)
target_include_directories(BeamSpectra PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(BeamSpectra PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS BeamSpectra DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(BeamSpectra        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME BeamSpectra SOVERSION 0.0.0)
set_target_properties(BeamSpectra        PROPERTIES DEFINE_SYMBOL "")
endif()
