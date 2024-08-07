set(Spectra_SOURCES Monochromatic.C
                    Laser_Backscattering.C
                    EPA.C
                    EPA_FF.C
                    EPA_Spectra_Plotter.C
                    Pomeron.C
                    Reggeon.C
                    DM_beam.C)

set(Spectra_HEADERS Monochromatic.H
                    Laser_Backscattering.H
                    EPA.H
                    EPA_FF.H
                    EPA_Spectra_Plotter.H
                    Pomeron.H
                    Reggeon.H
                    DM_beam.H)

sherpa_create_git_info(Spectra ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(BeamSpectra SHARED ${Spectra_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(BeamSpectra BeamSpectra)
add_dependencies(BeamSpectra Spectra_Git_Info)
