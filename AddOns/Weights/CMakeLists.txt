add_library(SherpaWeightsOutput SHARED  Output_Weights.C)
install(TARGETS SherpaWeightsOutput DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaWeightsOutput       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaWeightsOutput SOVERSION 0.0.0)
