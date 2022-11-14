##### Automatically or semiautomaticaly generated, do not edit.#########
#SpinCorrelations/Makefile.am
 
 #pkglib_LTLIBRARIES = libMEToolsSpinCorrelations.la)
 set(GITTAG  METOOLS/SpinCorrelations)
set(SpinCorrelations_SOURCES  Amplitude2_Tensor.C Amplitude2_Matrix.C Decay_Matrix.C Spin_Density.C)
set(SpinCorrelations_HEADERS Amplitude2_Tensor.H Amplitude2_Matrix.H Decay_Matrix.H Spin_Density.H)
creategitinfo(SpinCorrelations ../SpinCorrelations/)
list(TRANSFORM SpinCorrelations_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../SpinCorrelations/")
list(APPEND SpinCorrelations_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(METoolsSpinCorrelations SHARED ${SpinCorrelations_SOURCES} )
target_include_directories(METoolsSpinCorrelations PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(METoolsSpinCorrelations)
install(TARGETS METoolsSpinCorrelations DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(METoolsSpinCorrelations        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME METoolsSpinCorrelations SOVERSION 0.0.0)
set_target_properties(METoolsSpinCorrelations        PROPERTIES DEFINE_SYMBOL "")
