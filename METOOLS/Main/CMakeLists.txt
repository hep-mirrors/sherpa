set(Main_SOURCES Polarization_Tools.C
                 XYZFuncs.C
                 Polarization_Index.C
                 Spin_Structure.C
                 Partial_Amplitude_Base.C
                 Three_Particle_Amplitudes.C
                 Four_Particle_Amplitudes.C
                 SpinFuncs.C)

set(Main_HEADERS Polarization_Tools.H
                 XYZFuncs.H
                 Polarization_Index.H
                 Spin_Structure.H
                 Partial_Amplitude_Base.H
                 Three_Particle_Amplitudes.H
                 Four_Particle_Amplitudes.H
                 SpinFuncs.H)

sherpa_create_git_info(Main)
 
add_library(METoolsMain SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(METoolsMain METoolsMain)
