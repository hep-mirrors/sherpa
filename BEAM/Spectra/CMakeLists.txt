##### Automatically or semiautomaticaly generated, do not edit.#########
#Spectra/Makefile.am
 
 #pkglib_LTLIBRARIES = libBeamSpectra.la)
 set(GITTAG  BEAM/Spectra)
set(Spectra_SOURCES  Monochromatic.C Laser_Backscattering.C EPA.C DM_beam.C Spectrum_Reader.C)
set(Spectra_HEADERS Monochromatic.H Laser_Backscattering.H EPA.H DM_beam.H Spectrum_Reader.H)
sherpa_create_git_info(Spectra ../Spectra/)
list(TRANSFORM Spectra_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Spectra/")
list(APPEND Spectra_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(BeamSpectra SHARED ${Spectra_SOURCES} )
sherpa_handle_shared_library(BeamSpectra BeamSpectra)
