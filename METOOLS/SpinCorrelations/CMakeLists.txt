set(SpinCorrelations_SOURCES Amplitude2_Tensor.C
                             Amplitude2_Matrix.C
                             Decay_Matrix.C
                             Spin_Density.C
                             PolWeights_Map.C
                             Polarized_CrossSections_Handler.C)

set(SpinCorrelations_HEADERS Amplitude2_Tensor.H
                             Amplitude2_Matrix.H
                             Decay_Matrix.H
                             Spin_Density.H
                             PolWeights_Map.H
                             Polarized_CrossSections_Handler.H)

sherpa_create_git_info(SpinCorrelations ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(METoolsSpinCorrelations SHARED ${SpinCorrelations_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(METoolsSpinCorrelations METoolsSpinCorrelations)
add_dependencies(METoolsSpinCorrelations SpinCorrelations_Git_Info)
