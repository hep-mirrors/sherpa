set(Spectra_SOURCES Monochromatic.C
                    Laser_Backscattering.C
                    EPA.C
                    Pomeron.C
                    Reggeon.C
                    Fixed_Target.C
                    DM_beam.C)

set(Spectra_HEADERS Monochromatic.H
                    Laser_Backscattering.H
                    EPA.H
                    Pomeron.H
                    Reggeon.H
                    Fixed_Target.H
                    DM_beam.H)

sherpa_create_git_info(Spectra)

add_library(BeamSpectra SHARED ${Spectra_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(BeamSpectra BeamSpectra)
