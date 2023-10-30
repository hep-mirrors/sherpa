set(Spectra_SOURCES Monochromatic.C
                    Laser_Backscattering.C
                    EPA.C
                    DM_beam.C
                    Fixed_Target.C)

set(Spectra_HEADERS Monochromatic.H
                    Laser_Backscattering.H
                    EPA.H
                    DM_beam.H
                    Fixed_Target.H)

sherpa_create_git_info(Spectra)
 
add_library(BeamSpectra SHARED ${Spectra_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(BeamSpectra BeamSpectra)
