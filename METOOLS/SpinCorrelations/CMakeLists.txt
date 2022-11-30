##### Automatically or semiautomaticaly generated, do not edit.#########
#SpinCorrelations/Makefile.am
 
set(GITTAG  METOOLS/SpinCorrelations)

set(SpinCorrelations_SOURCES Amplitude2_Tensor.C
                             Amplitude2_Matrix.C
                             Decay_Matrix.C
                             Spin_Density.C)

set(SpinCorrelations_HEADERS Amplitude2_Tensor.H
                             Amplitude2_Matrix.H
                             Decay_Matrix.H
                             Spin_Density.H)

sherpa_create_git_info(SpinCorrelations)
 
add_library(METoolsSpinCorrelations SHARED ${SpinCorrelations_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(METoolsSpinCorrelations METoolsSpinCorrelations)
